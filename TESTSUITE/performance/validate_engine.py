# SPDX-License-Identifier: GPL-3.0-or-later
"""Small transient regression fixtures for experimental openEMS engine builds.

Use the same script with each library build and compare the saved numerical
arrays. This is a solver regression, not an antenna design or convergence test.
"""
import argparse
from pathlib import Path
import json

import numpy as np


def simulate(args):
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    fdtd = openEMS(NrTS=2500, EndCriteria=0, OverSampling=4)
    fdtd.SetGaussExcite(8e9, 6e9)
    bc = ['PML_8'] * 6
    if args.case == 'mixed':
        bc = ['PML_8', 'PML_8', 'MUR', 'MUR', 'PEC', 'MUR']
    fdtd.SetBoundaryCond(bc)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    mesh.SetLines('x', np.linspace(-24, 24, 65))
    mesh.SetLines('y', np.linspace(-21, 21, 57))
    mesh.SetLines('z', np.unique(np.r_[np.linspace(-15, 18, 45), 0, 0.5, 1]))
    csx.AddMaterial('substrate', epsilon=3.66, kappa=0.002).AddBox(
        [-15, -12, 0], [15, 12, 1])
    if args.case == 'dispersive':
        from CSXCAD.CSProperties import CSPropLorentzMaterial
        material = CSPropLorentzMaterial(csx.GetParameterSet(), name='lorentz',
                                        epsilon=1.8, mue=1.2, order=2)
        csx.AddProperty(material)
        material.SetDispersiveMaterialProperty(0, eps_plasma=3e9, eps_pole_freq=12e9,
                                              eps_relax=50e-12, mue_plasma=1e9,
                                              mue_pole_freq=10e9, mue_relax=60e-12)
        material.SetDispersiveMaterialProperty(1, eps_plasma=2e9, eps_pole_freq=0,
                                              eps_relax=80e-12, mue_plasma=0.5e9,
                                              mue_pole_freq=0, mue_relax=80e-12)
        material.AddBox([-6, -6, 3], [6, 6, 6])
    metal = csx.AddMetal('copper')
    metal.AddBox([-15, -12, 0], [15, 12, 0], priority=10)
    csx.AddConductingSheet('patch', conductivity=5.8e7, thickness=35e-6).AddBox(
        [-4.5, -3.75, 1], [4.5, 3.75, 1], priority=10)
    ports = [fdtd.AddLumpedPort(1, 50, [-1.5, 0, 0], [-1.5, 0, 1],
                               'z', excite=1, priority=20),
             fdtd.AddLumpedPort(2, 75, [3, 0, 0], [3, 0, 1],
                               'z', excite=0, priority=20)]
    # Add a reactive extension as well as the ordinary lumped port resistance.
    csx.AddLumpedElement('load', ny=2, caps=True, R=120, L=1e-9, C=0.2e-12).AddBox(
        [0, 2.25, 0], [0, 2.25, 1], priority=20)
    csx.AddDump('electric', dump_type=10, file_type=1, frequency=[8e9]).AddBox(
        [-9, -9, 3], [9, 9, 3])
    fdtd.Write2XML(str(args.output / 'fixture.xml'))
    fdtd.Run(str(args.output), cleanup=False, numThreads=args.threads, verbose=1, engine=args.engine)
    maps = Path('/proc/self/maps')
    if maps.exists():
        (args.output / 'loaded-libraries.txt').write_text('\n'.join(
            line for line in maps.read_text().splitlines() if 'libopenEMS' in line))
    frequency = np.linspace(3e9, 13e9, 101)
    arrays = {'frequency': frequency}
    for number, port in enumerate(ports, 1):
        port.CalcPort(str(args.output), frequency, ref_impedance=50)
        arrays[f'port{number}_incident'] = port.uf_inc
        arrays[f'port{number}_reflected'] = port.uf_ref
    # Persist time traces without parsing header timestamps.
    for path in sorted(args.output.glob('port_*')):
        if path.is_file():
            arrays[path.name] = np.loadtxt(path, comments='%')
    np.savez_compressed(args.output / 'regression.npz', **arrays)
    (args.output / 'fixture.json').write_text(json.dumps(vars(args), default=str, indent=2))


def compare(args):
    import h5py
    report = []

    def check(name, a, b):
        if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise ValueError(f'Invalid shape or nonfinite values: {name}')
        scale = float(np.max(np.abs(a))) if a.size else 0.0
        error = float(np.max(np.abs(a - b))) if a.size else 0.0
        passed = bool(error <= 1e-6 * scale + 1e-20)
        report.append(dict(array=name, peak=scale, max_abs_error=error,
                           relative_peak_error=error / scale if scale else 0,
                           identical=bool(np.array_equal(a, b)), passed=passed))

    with np.load(args.reference / 'regression.npz') as ref, np.load(args.output / 'regression.npz') as new:
        if set(ref.files) != set(new.files):
            raise ValueError('Different regression arrays')
        for name in ref.files:
            check(name, ref[name], new[name])
        if (np.max(np.abs(ref['port1_incident'])) < 1e-15 or
                np.max(np.abs(ref['port_ut_1'][:, 1])) < 1e-4):
            raise ValueError('Trivial excitation cannot validate a solver')
    with h5py.File(args.reference / 'electric.h5') as ref, h5py.File(args.output / 'electric.h5') as new:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                check('electric.h5/' + name, np.asarray(obj[()]), np.asarray(new[name][()]))
        ref.visititems(visit)
    destination = args.output / 'comparison.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(f'{sum(r["passed"] for r in report)}/{len(report)} arrays passed; {destination}')
    if not all(r['passed'] for r in report):
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run', 'compare'])
    parser.add_argument('--case', choices=['pml', 'mixed', 'dispersive'], default='pml')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--engine', choices=['fastest', 'basic', 'sse', 'sse-compressed', 'multithreaded'],
                        default='fastest')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    if args.action == 'compare' and args.reference is None:
        parser.error('--reference is required for compare')
    if args.action == 'run':
        simulate(args)
    else:
        compare(args)
