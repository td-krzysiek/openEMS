# SPDX-License-Identifier: GPL-3.0-or-later
"""Repeated/manual/startup gates preserve a live simulation's state.

Uses the existing small startup fixture; no private antenna geometry. Compare
each engine with its own uninterrupted control: rtol1e-6/peak-atol1e-7.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import time
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location('startup_fixture', Path(__file__).with_name('startup_gate.py'))
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--binary', required=True)
    p.add_argument('--engine', default='basic')
    p.add_argument('--output', required=True)
    p.add_argument('--solver-option', action='append', default=[])
    a = p.parse_args()
    root = Path(a.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    command = [str(Path(a.binary).resolve()), 'model.xml', '--engine=' + a.engine,
               '--numThreads=2', '--exact-endcriteria', '-v'] + a.solver_option
    results = {}
    control = root / 'control'
    control.mkdir()
    fixture.model(control / 'model.xml', decay=True)
    with (control / 'solver.log').open('w') as log:
        subprocess.run(command, cwd=control, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=90)
    for label in ['automatic', 'manual-startup', 'timeout']:
        folder = root / label
        folder.mkdir()
        fixture.model(folder / 'model.xml', decay=True)
        token = label
        options = ['--runtime-control', '--decay-gate-id=' + token, '--dft-snapshots',
                   '--decay-gate-timeout=' + ('1' if label == 'timeout' else '30')]
        if label in ['automatic', 'timeout']:
            options += ['--decay-gates=' + ('-10' if label == 'timeout' else '-10,-20,-30,-40,-50')]
        if label == 'manual-startup':
            options += ['--startup-gate=64', '--startup-gate-id=startup', '--startup-gate-timeout=30']
        log_path = folder / 'solver.log'
        seen = []
        with log_path.open('w') as log:
            process = subprocess.Popen(command + options, cwd=folder, stdout=log, stderr=subprocess.STDOUT)
            try:
                if label == 'manual-startup':
                    fixture.wait_for(folder / 'startup-gate.xml', process, log_path)
                    fixture.atomic(folder / 'RUNTIME_PAUSE', token + '-manual')
                    fixture.atomic(folder / 'STARTUP_CONTINUE', 'startup')
                deadline = time.monotonic() + 90
                while process.poll() is None:
                    assert time.monotonic() < deadline, log_path.read_text()[-3000:]
                    witness = folder / 'runtime-gate.xml'
                    if witness.exists():
                        gate = ET.parse(witness).getroot()
                        if gate.get('state') == 'paused' and gate.get('id') not in seen:
                            seen.append(gate.get('id'))
                            steps = gate.get('timestep')
                            assert gate.get('snapshotAvailable') == 'true'
                            snapshot = folder / 'dft-snapshots' / gate.get('id') / 'snapshot.xml'
                            assert ET.parse(snapshot).getroot().get('timestep') == steps
                            before = (folder / 'probe').read_bytes()
                            time.sleep(.15)
                            assert (folder / 'probe').read_bytes() == before, 'FDTD advanced while paused'
                            if label != 'timeout':
                                fixture.atomic(folder / 'RUNTIME_CONTINUE', gate.get('id'))
                    time.sleep(.01)
                assert process.returncode == 0, log_path.read_text()[-3000:]
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
        if label == 'automatic':
            assert seen == [token + '-d' + str(n) for n in [10, 20, 30, 40]], seen
            assert not (folder / 'dft-snapshots' / (token + '-d50')).exists(), 'Duplicated final export'
        elif label == 'timeout':
            assert ET.parse(folder / 'runtime-gate.xml').getroot().get('state') == 'timed-out'
        else:
            assert seen == [token + '-manual']
        assert log_path.read_text().count('Creation time for operator:') == 1
        pattern = r'end-criteria of -50(?:\.0+)?dB reached after (\d+) timesteps'
        assert re.search(pattern, log_path.read_text())[1] == re.search(pattern, (control / 'solver.log').read_text())[1]
        results[label] = dict(gates=seen, **fixture.compare(folder, control))
    # A manual stop calls normal PostProcess; malformed/stale IDs fail closed.
    for label in ['stop', 'abort', 'wrong-id']:
        folder = root / label
        folder.mkdir()
        fixture.model(folder / 'model.xml', decay=True)
        log_path = folder / 'solver.log'
        with log_path.open('w') as log:
            process = subprocess.Popen(command + ['--decay-gates=-10', '--decay-gate-id=test',
                '--decay-gate-timeout=30', '--dft-snapshots'], cwd=folder, stdout=log, stderr=subprocess.STDOUT)
            try:
                fixture.wait_for(folder / 'runtime-gate.xml', process, log_path)
                if label == 'abort':
                    (folder / 'ABORT').touch()
                else:
                    fixture.atomic(folder / ('RUNTIME_CONTINUE' if label == 'wrong-id' else 'RUNTIME_STOP'),
                                   'wrong' if label == 'wrong-id' else 'test-d10')
                code = process.wait(timeout=30)
                assert code != 0 if label == 'wrong-id' else code == 0
                if label != 'wrong-id':
                    assert 'simulation aborted after' in log_path.read_text()
                results[label] = dict(exit_code=code)
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
    (root / 'result.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
