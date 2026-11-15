# SPDX-License-Identifier: GPL-3.0-or-later
"""Matched-input CPU/device PML + PEC dipole radiation qualification.

Fixture copied explicitly from the unchanged Cartesian_Subgrid_Radiation.py
control. FDTD uses the selected backend; DFT/NF2FF extraction remains on host.

Two solid rectangular PEC arms and one matched local-edge 50-ohm port.
The 0.9/1/1.1 GHz fixture is a public solver control, not a PCB antenna.
All six NF2FF faces and interpolation neighbours remain in homogeneous,
unrefined air, outside the local box and inside PML. No mirror faces.

The ordinary Yee comparison keeps the exact same tensor mesh and port cell.
The local/tensor comparison keeps the same fine cells around metal and port.
The accepted local/tensor grids use 5 mm interior Z intervals while retaining
the 10 mm gap edge. A preserved coarser-Z case must fail the mesh-sensitivity
gate. The finer-interior reference additionally refines transverse air, with
the same port cells and metal. This is a finite sensitivity control, not a
demonstration of asymptotic convergence for arbitrary antennas.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from _subgrid_device_support import run_native_device
from openEMS.physical_constants import Z0

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_Device_Radiation')
Path(Sim_Path).mkdir(parents=True, exist_ok=True)
FREQ = np.array([.9e9, 1e9, 1.1e9])
THETA = np.arange(0, 181, 3.)
PHI = np.arange(0, 361, 15.)
DT = 5e-12
METRICS = {}
RESULTS = {}
PHYSICAL = os.environ.get('OPENEMS_SUBGRID_DEVICE', '')
assert PHYSICAL in ('', 'metal', 'cuda')


def require(condition, message, value=None):
    assert condition, (message, value)


def rel_error(value, reference):
    peak = float(np.max(abs(reference)))
    require(peak > 1e-30, 'unexcited comparison')
    return float(np.max(abs(value-reference))/peak)


def integrate(values, coordinates, axis=-1):
    values = np.moveaxis(values, axis, -1)
    return np.sum(.5*(values[..., 1:]+values[..., :-1])*np.diff(coordinates), axis=-1)


def audit_geometry(path):
    """Read native XML independently of the objects used by the generator."""
    tree = ET.parse(path)
    boxes = []
    for metal in tree.findall('.//Properties/Metal'):
        for box in metal.findall('./Primitives/Box'):
            a = box.find('P1').attrib
            b = box.find('P2').attrib
            lo = np.array([float(a[n]) for n in 'XYZ'])
            hi = np.array([float(b[n]) for n in 'XYZ'])
            boxes.append((np.minimum(lo, hi), np.maximum(lo, hi)))
    require(len(boxes) == 2, 'expected two native PEC volumes', boxes)
    def span(element):
        box = element.find('./Primitives/Box')
        return [np.array([float(box.find(point).get(n)) for n in 'XYZ']) for point in ('P1', 'P2')]

    properties = tree.find('.//Properties')
    resistor = properties.find('LumpedElement')
    source = properties.find('Excitation')
    probes = {p.get('Name'): p for p in properties.findall('ProbeBox')}
    require(resistor is not None and source is not None, 'missing native source/resistor')
    a, b = span(source)
    require(np.array_equal(a, [0, 0, -5]) and np.array_equal(b, [0, 0, 5]), 'unexpected native source span')
    require(all(np.array_equal(x, y) for x, y in zip(span(resistor), (a, b))), 'resistor/source spans differ')
    require(all(np.array_equal(x, y) for x, y in zip(span(probes['port_ut_1']), (a, b))), 'voltage/source spans differ')
    require(all(np.array_equal(x, (a+b)/2) for x in span(probes['port_it_1'])), 'current probe is not centred')
    require(resistor.get('Direction') == probes['port_it_1'].get('NormDir') == '2' and
            float(resistor.get('R')) == 50 and float(probes['port_ut_1'].get('Weight')) == -1 and
            float(probes['port_it_1'].get('Weight')) == 1 and
            np.array_equal(np.fromstring(source.get('Excite'), sep=','), [0, 0, -1]),
            'native port direction, sign or resistance differs')
    for z in (a[2], b[2]):
        # The source terminal is in the interior of a finite metal end face.
        contacts = [lo[0] < 0 < hi[0] and lo[1] < 0 < hi[1] and
                    (lo[2] == z or hi[2] == z) for lo, hi in boxes]
        require(sum(contacts) == 1, 'missing/ambiguous finite terminal contact', z)
    require(not any(lo[2] < 0 < hi[2] for lo, hi in boxes), 'PEC shorts the port gap')
    return dict(pec_boxes_mm=[dict(start=lo.tolist(), stop=hi.tolist()) for lo, hi in boxes],
                source_start_mm=a.tolist(), source_stop_mm=b.tolist(),
                current_plane_mm=((a+b)/2).tolist(), finite_contacts=True,
                source_resistor_voltage_current_agree=True)


def setup(kind, pml=10, steps=4000, record_td=False):
    fdtd = openEMS(NrTS=steps, EndCriteria=1e-30, TimeStep=DT, OverSampling=4)
    fdtd.SetGaussExcite(1e9, .5e9, dc_free=True)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    end = 120+10*pml
    xy = np.arange(-end, end+1, 10.)
    z = np.arange(-145-10*pml, 146+10*pml, 10.)
    if kind in ('tensor', 'yee', 'tensor-zfine', 'yee-zfine'):
        xy = np.unique(np.r_[xy, np.arange(-30, 31, 5.)])
    elif kind == 'fine-interior':
        xy = np.unique(np.r_[xy, np.arange(-120, 121, 5.)])
    if kind.endswith('-zfine') or kind == 'fine-interior':
        extra = np.arange(-145, 146, 5.)
        z = np.unique(np.r_[z, extra[extra != 0]])  # keep the single gap edge
    for axis in 'xy':
        mesh.SetLines(axis, xy)
    mesh.SetLines('z', z)
    fdtd.SetBoundaryCond(['PML_'+str(pml)]*6)
    if kind not in ('yee', 'yee-base', 'yee-zfine'):
        fdtd.AddCartesianSubgrid(box=[-30, -30, -75, 30, 30, 75],
                                ratio=[2, 2, 1] if kind.startswith('local') else [1, 1, 1])
    metal = csx.AddMetal('dipole')
    metal.AddBox([-10, -10, -65], [10, 10, -5], priority=10)
    metal.AddBox([-10, -10, 5], [10, 10, 65], priority=10)
    port = fdtd.AddLumpedPort(1, 50, [0, 0, -5], [0, 0, 5], 'z', 1)
    # Passive arm contours on dual Z planes. In the fine metal/port region,
    # +/-17.5 are transverse dual lines and +/-20 are dual Z planes.
    # Sampling at base Z planes instead would choose one side of a dual chain.
    dual_z = .5*(z[:-1]+z[1:])
    plane_upper = min(dual_z[dual_z > 0], key=lambda value: (abs(value-20), value))
    for plane, label in [(-plane_upper, 'lower'), (plane_upper, 'upper')]:
        csx.AddProbe('arm_'+label+'_i', p_type=1, norm_dir=2).AddBox(
            [-17.5, -17.5, plane], [17.5, 17.5, plane])
    positions = {'near_e': [40, 20, 25], 'near_h': [45, 25, 30]}
    for name, position in positions.items():
        kind_number = 2 if name == 'near_e' else 3
        csx.AddProbe(name, p_type=kind_number).AddPoint(position)
        if record_td:
            csx.AddProbe(name+'_fd', p_type=kind_number, frequency=FREQ).AddPoint(position)
    # Node-interpolated E/H at the same physical point across mesh variants.
    # Raw vector probes above retain component-dependent Yee positions.
    for name, field_type in [('near_node_e', 10), ('near_node_h', 11)]:
        csx.AddDump(name, dump_type=field_type, dump_mode=1, file_type=1, frequency=FREQ).AddBox(
            [40, 20, 25], [40, 20, 25])
    boxes = [fdtd.CreateNF2FFBox(name, start=[-span, -span, -span-5],
             stop=[span, span, span+5], frequency=FREQ)
             for name, span in [('inner', 80), ('outer', 100)]]
    if record_td:
        boxes.append(fdtd.CreateNF2FFBox('inner_td', start=[-80, -80, -85], stop=[80, 80, 85]))
    metadata = dict(kind=kind, pml_cells=pml, steps=steps, timestep_seconds=DT,
        frequencies_hz=FREQ.tolist(), base_lines_mm=[xy.tolist(), xy.tolist(), z.tolist()],
        local_box_mm=[-30, -30, -75, 30, 30, 75],
        local_ratio=[2, 2, 1] if kind.startswith('local') else [1, 1, 1],
        port=dict(start_mm=[0, 0, -5], stop_mm=[0, 0, 5], resistance_ohm=50),
        arm_current=dict(normal='+Z on both arms',
                         contours_mm=[[-17.5, 17.5, -plane_upper], [-17.5, 17.5, plane_upper]],
                         note='Symmetric dual-Z planes; finer-interior planes differ and are not compared across meshes.'),
        near_probe_positions_mm=positions,
        scope='Solid PEC dipole in vacuum; not the project antenna or a conductor-loss model.')
    return fdtd, port, boxes, metadata


def audit_surfaces(path, name, metadata):
    """Check actual HDF5 coverage/coordinates before allowing NF2FF analysis."""
    planes = []
    for face in range(6):
        meshes = []
        normal = face//2
        for field in 'EH':
            source = path/('{}_{}_{}.h5'.format(name, field, face))
            require(source.is_file(), 'missing NF2FF face', str(source))
            with h5py.File(source) as data:
                axes = [np.array(data['Mesh/'+n]) for n in 'xyz']
                require(len(axes[normal]) == 1, 'NF2FF face is not planar', face)
                meshes.append(axes)
                if name != 'inner_td':
                    require(np.allclose(data['FieldData/FD'].attrs['frequency'], FREQ), 'wrong FD frequencies')
                    for index in range(len(FREQ)):
                        values = np.array(data['FieldData/FD/f'+str(index)])
                        require(np.isfinite(values).all() and np.max(abs(values)) > 1e-30,
                                'nonfinite or unexcited NF2FF face', (name, face, field))
        for a, b in zip(*meshes):
            require(np.array_equal(a, b), 'E/H NF2FF meshes differ', face)
        coordinate = float(meshes[0][normal][0])
        # Include interpolation neighbours: each face is outside the whole
        # refined region and metal, and at least one coarse cell from PML.
        region = .075 if normal == 2 else .030
        pml_inner = .145 if normal == 2 else .120
        require(abs(coordinate) >= region+.010-1e-12 and abs(coordinate) < pml_inner-.010,
                'NF2FF stencil approaches refinement/PML', (face, coordinate))
        planes.append(coordinate)
    require(all(planes[2*n] < 0 < planes[2*n+1] for n in range(3)), 'NF2FF faces do not enclose source')
    return planes


def check_probe_timing(path):
    metrics = {}
    for name in ('near_e', 'near_h'):
        td = np.loadtxt(path/name, comments='%')
        fd = np.loadtxt(path/(name+'_fd_FD'), comments='%')
        values = fd[:, 1::2]+1j*fd[:, 2::2]
        dt = td[1, 0]-td[0, 0]
        dft = 2*dt*(np.exp(-2j*np.pi*FREQ[:, None]*td[:, 0]) @ td[:, 1:])
        error = rel_error(dft, values)
        # FD and TD samplers need not use identical intervals. The complete
        # pulse is band limited; 0.1% is below the wrong-half-step control.
        require(error < .001, 'FD/TD probe disagreement', (name, error))
        shifted = 2*dt*(np.exp(-2j*np.pi*FREQ[:, None]*(td[:, 0]+DT/2)) @ td[:, 1:])
        negative = rel_error(shifted, values)
        require(negative > .005, 'wrong half-step was not detected', (name, negative))
        metrics[name] = dict(fd_td_relative_error=error, wrong_half_step_error=negative,
                             first_time_seconds=float(td[0, 0]), sample_interval_seconds=float(dt))
    return metrics


def run(label, kind, pml=10, steps=4000, record_td=False, engine='basic'):
    path = Path(Sim_Path)/label
    path.mkdir(exist_ok=False)
    fdtd, port, boxes, metadata = setup(kind, pml, steps, record_td)
    for prop in fdtd.GetCSX().GetAllProperties():
        prop.SetColor([96, 128, 160])
    fdtd.Write2XML(str(path/'input.xml'))
    metadata['native_contact_audit'] = audit_geometry(path/'input.xml')
    metadata['input_sha256'] = hashlib.sha256((path/'input.xml').read_bytes()).hexdigest()
    (path/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    previous = os.getcwd()
    try:
        run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1, exact_endcriteria=True)
    finally:
        os.chdir(previous)
    port.CalcPort(str(path), FREQ)
    require(np.isfinite(port.uf_tot).all() and np.isfinite(port.if_tot).all(), 'nonfinite port')
    require(np.min(port.P_acc) > 0 and np.min(abs(port.uf_inc)) > 1e-20, 'unexcited/nonpassive port')
    data = dict(frequencies=FREQ, theta=THETA, phi=PHI, z=port.uf_tot/port.if_tot,
                gamma=port.uf_ref/port.uf_inc, u_inc=port.uf_inc, pacc=port.P_acc,
                u=port.uf_tot, i=port.if_tot, pinc=port.P_inc)
    for name in ('near_node_e', 'near_node_h'):
        with h5py.File(path/(name+'.h5')) as field:
            coordinates = [float(field['Mesh/'+axis][0]) for axis in 'xyz']
            require(np.allclose(coordinates, [.040, .020, .025], atol=1e-12, rtol=0),
                    'near-field observation moved', (name, coordinates))
            data[name] = np.array([np.array(field['FieldData/FD/f'+str(i)]).reshape(3) for i in range(len(FREQ))])
    metrics = {}
    for box in boxes:
        planes = audit_surfaces(path, box.name, metadata)
        result = box.CalcNF2FF(str(path), FREQ, THETA, PHI, read_cached=False)
        data[box.name+'_et'] = np.array(result.E_theta)
        data[box.name+'_ep'] = np.array(result.E_phi)
        data[box.name+'_prad'] = np.array(result.Prad)
        data[box.name+'_dmax'] = np.array(result.Dmax)
        data[box.name+'_realized_gain'] = np.array(result.Dmax)*np.array(result.Prad)/port.P_inc
        norm = np.array(result.E_norm)
        normalized = norm/np.max(norm, axis=(1, 2))[:, None, None]
        theta_rad, phi_rad = np.deg2rad(THETA), np.deg2rad(PHI)
        sphere = np.array([integrate(integrate(e**2/(2*Z0), phi_rad, axis=1)*
                          np.sin(theta_rad), theta_rad) for e in norm])
        efficiency = np.array(result.Prad)/port.P_acc
        power_error = float(np.max(abs(efficiency-1)))
        sphere_error = float(np.max(abs(sphere/np.array(result.Prad)-1)))
        # Same lossless/surface-integration gates as NF2FF_Antennas.py.
        require(power_error < .05, 'lossless radiation/accepted-power balance', (label, power_error))
        require(sphere_error < .02, 'sphere integral/surface flux mismatch', (label, sphere_error))
        require(np.max(abs(np.array(result.Dmax)/1.64-1)) < .10,
                'dipole directivity physical check', result.Dmax)
        with np.errstate(divide='ignore', invalid='ignore'):
            analytic = np.where(abs(np.sin(theta_rad)) < 1e-9, 0,
                                np.cos(np.pi/2*np.cos(theta_rad))/np.sin(theta_rad))
        pattern_error = float(np.max(abs(normalized-analytic[None, :, None])))
        # Finite rectangular arms are not an analytical infinitesimal wire.
        require(pattern_error < .10, 'dipole pattern physical check', (label, pattern_error))
        require(np.max(normalized[:, [0, -1], :]) < .1, 'missing dipole axial null', label)
        require(np.max(np.ptp(normalized, axis=2)) < .02, 'excessive phi amplitude asymmetry', label)
        metrics[box.name] = dict(surface_planes_m=planes, directivity=result.Dmax.tolist(),
            radiation_accepted_ratio=efficiency.tolist(), sphere_surface_ratio=(sphere/np.array(result.Prad)).tolist(),
            maximum_pattern_error=pattern_error, maximum_phi_spread=float(np.max(np.ptp(normalized, axis=2))))
    if record_td:
        metrics['probe_timing'] = check_probe_timing(path)
        # Compare the complex vector against the total field peak. Ephi is
        # nearly zero for a dipole: its own tiny peak is not a useful scale.
        scale = np.max(np.sqrt(abs(data['inner_et'])**2+abs(data['inner_ep'])**2))
        delta = np.sqrt(abs(data['inner_td_et']-data['inner_et'])**2+
                        abs(data['inner_td_ep']-data['inner_ep'])**2)
        metrics['td_fd_vector_error'] = float(np.max(delta)/scale)
        metrics['td_fd_power_error'] = rel_error(data['inner_td_prad'], data['inner_prad'])
        require(metrics['td_fd_vector_error'] < .001 and metrics['td_fd_power_error'] < .001,
                'TD/FD NF2FF mismatch', metrics)
        # A complete closed surface is mandatory for this qualification.
        # Exercise the audit against separate incomplete links, without
        # changing the accepted inputs or raw field files.
        negative = path/'negative-missing-face'
        negative.mkdir(exist_ok=True)
        for source in path.glob('inner_[EH]_*.h5'):
            if source.name != 'inner_H_5.h5':
                link = negative/source.name
                if link.is_symlink():
                    link.unlink()
                link.symlink_to(source.resolve())
        try:
            audit_surfaces(negative, 'inner', metadata)
        except AssertionError as error:
            require('missing NF2FF face' in str(error), 'wrong negative-control failure', str(error))
        else:
            raise AssertionError('incomplete NF2FF surface was accepted')
        # A real break between source and metal must fail the independent
        # native terminal audit before interpreting RF results.
        tree = ET.parse(path/'input.xml')
        tree.find('.//Properties/Metal/Primitives/Box/P2').set('Z', '-15')
        broken = path/'negative-contact.xml'
        tree.write(broken)
        try:
            audit_geometry(broken)
        except AssertionError as error:
            require('terminal contact' in str(error), 'wrong contact-audit failure', str(error))
        else:
            raise AssertionError('broken native terminal was accepted')
    # Two closed surfaces must give the same complex far field and power.
    shell_errors = {name: rel_error(data['outer_'+name], data['inner_'+name]) for name in ('et', 'prad')}
    # Ephi is nearly zero for a dipole; use the total co-polar peak as its scale.
    shell_errors['ep'] = float(np.max(abs(data['outer_ep']-data['inner_ep']))/np.max(abs(data['inner_et'])))
    require(max(shell_errors.values()) < .01, 'NF2FF enclosure dependence', (label, shell_errors))
    metrics['shell_errors'] = shell_errors
    arms = [np.loadtxt(path/('arm_'+name+'_i'), comments='%') for name in ('lower', 'upper')]
    arm_symmetry = rel_error(arms[0][:, 1], arms[1][:, 1])
    require(arm_symmetry < 2e-4, 'symmetric dual-plane arm-current mismatch', (label, arm_symmetry))
    metrics['arm_current_symmetry'] = arm_symmetry
    np.savez(path/'results.npz', **data)
    (path/'metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
    RESULTS[label] = data
    METRICS[label] = metrics
    print('RADIATION', label, 'Z', data['z'], 'D', data['inner_dmax'],
          'Prad/Pacc', data['inner_prad']/data['pacc'], flush=True)
    return data


def compare(label, value, reference, impedance=.5, gamma=.005, field=.02):
    errors = dict(z_ohm=float(np.max(abs(value['z']-reference['z']))),
                  complex_gamma=float(np.max(abs(value['gamma']-reference['gamma']))))
    require(errors['z_ohm'] < impedance and errors['complex_gamma'] < gamma,
            'matched-port comparison', (label, errors))
    for name in ('near_node_e', 'near_node_h'):
        errors[name] = rel_error(value[name]/value['u_inc'][:, None], reference[name]/reference['u_inc'][:, None])
        require(errors[name] < field, 'common-point complex near-field comparison', (label, errors))
    for box in ('inner', 'outer'):
        scale = float(np.max(abs(reference[box+'_et']/reference['u_inc'][:, None, None])))
        for component in ('et', 'ep'):
            a = value[box+'_'+component]/value['u_inc'][:, None, None]
            b = reference[box+'_'+component]/reference['u_inc'][:, None, None]
            errors[box+'_'+component] = float(np.max(abs(a-b))/scale)
            require(errors[box+'_'+component] < field, 'complex far-field comparison', (label, errors))
        for quantity in ('dmax', 'realized_gain'):
            key = box+'_'+quantity
            errors[key] = float(np.max(abs(value[key]/reference[key]-1)))
            require(errors[key] < field, 'radiation quantity comparison', (label, errors))
    METRICS[label] = errors
    print('COMPARE', label, errors, flush=True)
    return errors



DEVICE_METRICS = {}
for kind in ['identity', 'local-zfine']:
    labels = {}
    reference = run(kind+'-basic', kind, record_td=kind=='local-zfine')
    reference_path = Path(Sim_Path)/(kind+'-basic')
    for engine in (['gpu-reference','gpu'] if PHYSICAL else ['gpu-reference']):
        label = kind+'-'+engine
        value = run(label, kind, record_td=kind=='local-zfine', engine=engine)
        candidate_path = Path(Sim_Path)/label
        require((candidate_path/'input.xml').read_bytes()==(reference_path/'input.xml').read_bytes(),
                'CPU/device radiation inputs differ', label)
        errors = compare(label+'-vs-basic', value, reference, impedance=.01, gamma=2e-4, field=2e-4)
        for quantity in ['pacc','u','i','pinc','inner_prad','outer_prad']:
            errors[quantity] = rel_error(value[quantity],reference[quantity])
            require(errors[quantity] < 2e-4,'CPU/device port/power difference',(label,quantity,errors[quantity]))
        DEVICE_METRICS[label] = errors
    Path(Sim_Path,'device-metrics.json').write_text(json.dumps(
        {'physical':PHYSICAL,'errors':DEVICE_METRICS},indent=2)+'\n')
print('PASS: native device dipole ports, near fields, two closed NF2FF surfaces, vector/power/TD-FD controls',flush=True)
