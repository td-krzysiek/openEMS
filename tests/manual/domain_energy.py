# SPDX-License-Identifier: GPL-3.0-or-later
"""Check one/two-device energy-selected stop and accumulated complex fields."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess

spec = importlib.util.spec_from_file_location('fixture', Path(__file__).with_name('startup_gate.py'))
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--two-device-control', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()
root, control = Path(args.output).resolve(), Path(args.two_device_control).resolve()
root.mkdir(parents=True, exist_ok=False)
shutil.copy2(control/'model.xml', root/'model.xml')
with (root/'solver.log').open('w') as log:
    subprocess.run([str(Path(args.binary).resolve()), 'model.xml', '--engine=gpu',
        '--gpu-device=0', '--numThreads=2', '--exact-endcriteria', '-v'], cwd=root,
        stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
pattern = r'end-criteria of -50(?:\.0+)?dB reached after (\d+) timesteps'
one = int(re.search(pattern, (root/'solver.log').read_text())[1])
two = int(re.search(pattern, (control/'solver.log').read_text())[1])
assert one == two, (one, two)
result = dict(one_gpu_steps=one, two_gpu_steps=two, **fixture.compare(root, control))
(root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
