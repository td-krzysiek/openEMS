"""Generate actual CSXCAD Python material/XML roundtrip controls."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import CSXCAD
from CSXCAD import ContinuousStructure


def material(index):
    return dict(epsilon=[2, 2, 2] if index == 1 else [2, 3, 5],
                mue=[1.2, 1.2, 1.2] if index == 1 else [1.2, 1.8, 2.5],
                kappa=[0, .3, .9] if index == 0 else [.4, 0, 0] if index in [2, 3] else [0, 0, 0],
                sigma=[0, 2500, 7000] if index in [1, 3] else [0, 0, 0])


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_property(csx, values):
    require(csx.GetQtyProperties() == 1, 'Material property inventory changed')
    prop = csx.GetProperty(0)
    require(prop.GetIsotropy() is False, 'Explicit directional property lost isotropy flag')
    for name, expected in values.items():
        require(np.array_equal(prop.GetMaterialProperty(name), expected), 'Python material component changed: ' + name)


def main(folder):
    folder.mkdir(exist_ok=False)
    cases = []
    for index in range(5):
        values = material(index)
        unit = 1e-6 if index % 2 else .001
        origin = [31, -23, 17]
        csx = ContinuousStructure()
        grid = csx.GetGrid()
        grid.SetDeltaUnit(unit)
        for n, axis in enumerate('xyz'):
            grid.AddLine(axis, [(origin[n] + i) * .001 / unit for i in range(5)])
        prop = csx.AddMaterial('diagonal-material-' + str(index))
        prop.SetIsotropy(False)
        prop.SetMaterialProperty(**values)
        prop.AddBox(start=[v * .001 / unit for v in origin],
                    stop=[(v + 4) * .001 / unit for v in origin], priority=1)
        require(csx.Update() == '', 'Python fixture CAD update failed')
        check_property(csx, values)
        path = folder / ('case-' + str(index) + '.xml')
        require(csx.Write2XML(path), 'Python native XML writer failed')
        reloaded = ContinuousStructure()
        require(reloaded.ReadFromXML(str(path)) == '', 'Actual Python XML reload failed')
        check_property(reloaded, values)
        require(reloaded.GetGrid().GetDeltaUnit() == unit, 'Python XML reload changed drawing unit')
        roundtrip = folder / ('case-' + str(index) + '-roundtrip.xml')
        require(reloaded.Write2XML(roundtrip), 'Python roundtrip XML writer failed')
        cases.append(dict(index=index, values=values, unit=unit, origin=origin,
                          files=[path.name, roundtrip.name], expectedRejection='none'))
    # These are serialized by the actual Python API. Native setup must reject
    # their invalid finite directional coefficients with the intended reason.
    for name, property_name, axis, value, reason in [
        ('invalid-epsilon', 'epsilon', 1, -3, 'positive finite diagonal'),
        ('invalid-kappa', 'kappa', 1, -.3, 'finite nonnegative diagonal'),
        ('invalid-sigma', 'sigma', 2, -100, 'finite nonnegative diagonal')]:
        csx = ContinuousStructure()
        require(csx.ReadFromXML(str(folder / 'case-0.xml')) == '', 'Negative-control source reload failed')
        values = material(0)
        values[property_name][axis] = value
        csx.GetProperty(0).SetMaterialProperty(**{property_name: values[property_name]})
        path = folder / (name + '.xml')
        require(csx.Write2XML(path), 'Negative-control native XML writer failed')
        cases.append(dict(index=0, values=values, unit=.001, origin=[31, -23, 17],
                          files=[path.name], expectedRejection=reason))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result = dict(status='ACTUAL_PYTHON_API_XML_ROUNDTRIP_FINISHED_NATIVE_FIELDS_PENDING',
                  python=sys.version, csxcad=CSXCAD.__version__, numpy=np.__version__,
                  bindingFile=str(Path(CSXCAD.CSXCAD.__file__).resolve()), cases=cases,
                  files={path.name: sha(path) for path in sorted(folder.glob('*.xml'))})
    with (folder / 'generated.json').open('x') as stream:
        json.dump(result, stream, indent=2)
    print(result['status'])


if __name__ == '__main__':
    require(len(sys.argv) == 2, 'expected a new output directory')
    main(Path(sys.argv[1]))
