"""Native SI conductivity mapping and scoped CPU/device parity.

Closed-boundary isotropic constant boxes only. Conductivity+absorber histories
are deliberately rejected; surface/conductor/ADE models remain unsupported.
"""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
from _subgrid_device_support import run_native_device

HERE = Path(__file__).resolve().parent
ROOT = Path(tempfile.gettempdir()) / 'Cartesian_Subgrid_Conductivity'
ROOT.mkdir(exist_ok=True)
PHYSICAL = os.environ.get('OPENEMS_SUBGRID_DEVICE', '')
assert PHYSICAL in ('', 'metal', 'cuda')
source = HERE / 'Cartesian_Subgrid_Materials.py'
ns = {'__file__': str(source), '__name__': 'native_conductivity_material_fixture'}
exec(compile(source.read_text().split('\nfor label, options', 1)[0], str(source), 'exec'), ns)


def fixture(electric, magnetic):
    fdtd, csx = ns['fixture']()
    materials = csx.GetPropertiesByName('slab')
    assert len(materials) == 1, 'expected exactly one slab material'
    materials[0].SetMaterialProperty(kappa=electric, sigma=magnetic)
    # The higher-priority pocket remains lossless. The C++ controls verify its
    # actual SI map independently, including reversed insertion and units.
    return fdtd, csx


def run(electric, magnetic, name, engine):
    fdtd, csx = fixture(electric, magnetic)
    path = ROOT / (name+'-'+engine)
    previous = os.getcwd()
    try:
        run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1,
                          exact_endcriteria=True)
    finally:
        os.chdir(previous)
    values = {key: np.loadtxt(path/key, comments='%') for key in ['received', 'magnetic']}
    assert all(np.isfinite(value).all() and np.max(abs(value[:, 1:])) > 0 for value in values.values())
    return values


baseline = run(0, 0, 'zero', 'basic')
results = {}
for name, electric, magnetic in [('electric', .04, 0), ('magnetic', 0, 120), ('both', .04, 120)]:
    expected = run(electric, magnetic, name, 'basic')
    change = max(float(np.max(abs(expected[key][:, 1:]-baseline[key][:, 1:])) /
                       np.max(abs(baseline[key][:, 1:]))) for key in baseline)
    assert change > 1e-4, (name, 'material loss did not affect the native trace', change)
    for engine in ['gpu-reference'] + (['gpu'] if PHYSICAL else []):
        actual = run(electric, magnetic, name, engine)
        errors = {}
        for key in expected:
            assert expected[key].shape == actual[key].shape
            np.testing.assert_allclose(expected[key][:, 0], actual[key][:, 0], rtol=1e-12, atol=0)
            errors[key] = float(np.max(abs(expected[key][:, 1:]-actual[key][:, 1:])) /
                                np.max(abs(expected[key][:, 1:])))
        assert max(errors.values()) < 2e-4, (name, engine, errors)
        results[name+'-'+engine] = dict(errors=errors, loss_effect=change)

# A real, otherwise valid PML setup must reach the explicit combined-physics
# guard rather than failing from geometry, collar size or source placement.
source = HERE / 'Cartesian_Subgrid_Device.py'
other = {'__file__': str(source), '__name__': 'native_conductivity_pml_guard_fixture'}
exec(compile(source.read_text().split('\nresults = {}', 1)[0], str(source), 'exec'), other)
fdtd, csx = other['fixture'](steps=8, source_x=0)
fdtd.SetBoundaryCond(['PML_4']*6)
csx.AddMaterial('unsupported_loss_in_pml', epsilon=1, kappa=.01).AddBox([-16]*3, [16]*3)
previous = os.getcwd()
try:
    try:
        run_native_device(fdtd, ROOT/'reject-pml-loss', 'basic', '', numThreads=1, setup_only=True)
    except (ValueError, RuntimeError) as error:
        assert 'conductivity with absorber histories' in str(error), str(error)
    else:
        raise AssertionError('unsupported conductivity/PML histories were admitted')
finally:
    os.chdir(previous)
(ROOT/'metrics.json').write_text(json.dumps(dict(physical=PHYSICAL, cases=results,
    scope='Native constant SI kappa/sigma only; no surface/conductor/antenna qualification.'), indent=2)+'\n')
print('PASS: native constant conductivity mapping, physical field parity and absorber guard')
