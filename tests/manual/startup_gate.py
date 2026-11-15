# SPDX-License-Identifier: GPL-3.0-or-later
"""Native startup qualification: paused continuation versus uninterrupted FDTD.

Run with --binary PATH --engine basic|gpu-reference|gpu --output DIR.
Compare each engine with itself: rtol1e-6, absolute floor1e-7 times reference
peak; report measured errors. Test failure paths on the selected actual engine.
"""
import argparse
import importlib.util
import json
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


def model(path, steps=12000, decay=False):
    fdtd = openEMS(NrTS=steps, EndCriteria=1e-5 if decay else 0)
    fdtd.SetGaussExcite(3e9, 2e9)
    fdtd.SetBoundaryCond(["PML_8" if decay else "PEC"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis in "xyz":
        grid.AddLine(axis, np.linspace(-30, 30, 25))
    if decay:
        # Bulk conductivity damps the source's electrostatic residual as well.
        csx.AddMaterial("lossy", epsilon=1, kappa=.02).AddBox([-30]*3, [30]*3)
    csx.AddExcitation("source", exc_type=0, exc_val=[0, 0, 1]).AddBox([0, 0, -2.5], [0, 0, 2.5])
    csx.AddProbe("probe", p_type=2).AddPoint([5, 5, 0])
    for field, kind in [("E", 10), ("H", 11)]:
        dump = csx.AddDump("nf2ff_" + field, dump_type=kind, dump_mode=1,
                           file_type=1, frequency=[3e9, 5e9])
        for n in range(6):
            a, b = [-7.5] * 3, [7.5] * 3
            a[n // 2] = b[n // 2] = [-7.5, 7.5][n % 2]
            dump.AddBox(a, b)
    fdtd.Write2XML(str(path))


def fields(folder):
    result = {}
    for path in sorted(folder.glob("nf2ff_*.h5")):
        with h5py.File(path, "r") as h:
            for name, dataset in h["FieldData/FD"].items():
                result[path.name + "/" + name] = dataset[()]
    assert len(result) > 0
    return result


def compare(left, right):
    a, b = fields(left), fields(right)
    assert a.keys() == b.keys()
    errors = []
    nonzero = 0
    for key in a:
        assert np.isfinite(a[key]).all() and np.isfinite(b[key]).all()
        peak = float(np.max(abs(b[key])))
        np.testing.assert_allclose(a[key], b[key], rtol=1e-6, atol=peak * 1e-7)
        errors.append(float(np.max(abs(a[key] - b[key]))))
        nonzero += int(peak > 0)
    assert nonzero > 0
    x = np.loadtxt(left / "probe", comments="%")
    y = np.loadtxt(right / "probe", comments="%")
    assert x.shape == y.shape
    np.testing.assert_allclose(x, y, rtol=1e-6, atol=np.max(abs(y)) * 1e-7)
    return dict(max_field_absolute_error=max(errors), nonzero_arrays=nonzero,
                max_probe_absolute_error=float(np.max(abs(x-y))))


def atomic(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(value + "\n")
    temporary.replace(path)


def wait_for(path, process, log, timeout=30):
    deadline = time.monotonic() + timeout
    while not path.exists():
        assert process.poll() is None, log.read_text()[-3000:]
        assert time.monotonic() < deadline, path
        time.sleep(.01)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--engine", default="basic")
    parser.add_argument("--output", required=True)
    parser.add_argument("--legacy", action="store_true")
    parser.add_argument("--library-path", help="Child executable LD_LIBRARY_PATH; keep model bindings isolated")
    args = parser.parse_args()
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    command = [str(Path(args.binary).resolve()), "model.xml", "--engine=" + args.engine,
               "--numThreads=2", "--exact-endcriteria", "-v"]
    if args.legacy:
        command.append("--legacyHDF5Dumps")
    if args.library_path:
        command = ["env", "LD_LIBRARY_PATH=" + args.library_path, *command]
    results = {}
    for label, decay, gate_step in [("fixed", False, 4096), ("energy", True, 128)]:
        control, live = root / (label + "-control"), root / (label + "-gated")
        for folder in [control, live]:
            folder.mkdir()
            model(folder / "model.xml", decay=decay)
        with (control / "solver.log").open("w") as log:
            subprocess.run(command, cwd=control, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=90)
        token = "qualification-" + label
        log_path = live / "solver.log"
        with log_path.open("w") as log:
            process = subprocess.Popen(command + ["--dft-snapshots", "--startup-gate=" + str(gate_step),
                "--startup-gate-id=" + token, "--startup-gate-timeout=30"], cwd=live,
                stdout=log, stderr=subprocess.STDOUT)
            try:
                wait_for(live / "startup-gate.xml", process, log_path)
                witness = ET.parse(live / "startup-gate.xml").getroot()
                assert witness.get("id") == token and witness.get("converged") == "false"
                actual = int(witness.get("timestep"))
                assert actual >= gate_step and actual < 12000
                probe = (live / "probe").read_bytes()
                assert len(np.loadtxt(live / "probe", comments="%")) > 2
                time.sleep(.25)
                assert process.poll() is None and (live / "probe").read_bytes() == probe
                # DFT export during the pause must not advance or reset the solver.
                atomic(live / "DFT_SNAPSHOT", token)
                snapshot = live / "dft-snapshots" / token
                wait_for(snapshot / "snapshot.xml", process, log_path)
                assert int(ET.parse(snapshot / "snapshot.xml").getroot().get("timestep")) == actual
                assert len(list(snapshot.glob("*.h5"))) == 12
                atomic(live / "STARTUP_CONTINUE", token)
                assert process.wait(timeout=90) == 0, log_path.read_text()[-3000:]
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
        text = log_path.read_text()
        assert text.count("Creation time for operator:") == 1
        assert text.count("Startup gate: accepted;") == 1
        if decay:
            import re
            pattern = r"end-criteria of -50(?:\.0+)?dB reached after (\d+) timesteps"
            assert re.search(pattern, text)[1] == re.search(pattern, (control / "solver.log").read_text())[1]
        results[label] = dict(actual_gate_step=actual, **compare(live, control))
    for label in ["reject", "timeout", "wrong-id", "abort"]:
        folder = root / label
        folder.mkdir()
        model(folder / "model.xml", steps=1000)
        log_path = folder / "solver.log"
        with log_path.open("w") as log:
            process = subprocess.Popen(command + ["--startup-gate=128", "--startup-gate-id=test",
                "--startup-gate-timeout=1"], cwd=folder, stdout=log, stderr=subprocess.STDOUT)
            try:
                wait_for(folder / "startup-gate.xml", process, log_path)
                if label == "reject":
                    atomic(folder / "STARTUP_REJECT", "test")
                elif label == "wrong-id":
                    atomic(folder / "STARTUP_CONTINUE", "stale-token")
                elif label == "abort":
                    (folder / "ABORT").touch()
                code = process.wait(timeout=10)
                text = log_path.read_text()
                assert "continuing without reinitialization" not in text
                if label == "abort":
                    assert code == 0 and "simulation aborted after" in text
                else:
                    assert code != 0 and "Startup gate:" in text
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
        results[label] = dict(exit_code=code, rejected_or_aborted=True)
    for label, flags, stale in [
        ("missing-id", ["--startup-gate=128"], None),
        ("invalid-id", ["--startup-gate=128", "--startup-gate-id=../bad"], None),
        ("past-cap", ["--startup-gate=1000", "--startup-gate-id=test"], None),
        ("stale", ["--startup-gate=128", "--startup-gate-id=test"], "STARTUP_CONTINUE"),
    ]:
        folder = root / label
        folder.mkdir()
        model(folder / "model.xml", steps=1000)
        if stale:
            atomic(folder / stale, "test")
        with (folder / "solver.log").open("w") as log:
            result = subprocess.run(command + flags, cwd=folder, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=10)
        text = (folder / "solver.log").read_text()
        assert result.returncode != 0 and "Creation time for operator:" not in text
        results[label] = dict(rejected_before_operator=True)
    # Preserve experimental conductor histories as well as ordinary Yee fields.
    fixture_path = Path(__file__).resolve().parents[2] / 'TESTSUITE/opaquesurface/test_gpu.py'
    spec = importlib.util.spec_from_file_location('opaque_fixture', fixture_path)
    fixture_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture_module)
    opaque = fixture_module.fixture(2, uniform=True)
    opaque_probes = [p.get('Name') for p in opaque.findall('ContinuousStructure/Properties/ProbeBox')]
    left, right = root / 'opaque-control', root / 'opaque-gated'
    for folder in [left, right]:
        folder.mkdir()
        ET.ElementTree(opaque).write(folder / 'model.xml')
    with (left / 'solver.log').open('w') as log:
        subprocess.run(command, cwd=left, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=90)
    log_path = right / 'solver.log'
    with log_path.open('w') as log:
        process = subprocess.Popen(command + ['--startup-gate=4096', '--startup-gate-id=opaque',
            '--startup-gate-timeout=30'], cwd=right, stdout=log, stderr=subprocess.STDOUT)
        try:
            wait_for(right / 'startup-gate.xml', process, log_path)
            atomic(right / 'STARTUP_CONTINUE', 'opaque')
            assert process.wait(timeout=90) == 0, log_path.read_text()[-3000:]
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    maximum = 0.
    nonzero = 0
    for name in opaque_probes:
        a, b = np.loadtxt(left / name, comments='%'), np.loadtxt(right / name, comments='%')
        assert a.shape == b.shape and np.isfinite(b).all()
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=np.max(abs(a))*1e-7)
        maximum = max(maximum, float(np.max(abs(a-b))))
        nonzero += int(np.any(abs(a[:, 1:]) > 0))
    assert nonzero > 0 and 'EXPERIMENTAL_OPAQUE_SURFACE_V1 active' in log_path.read_text()
    results['opaque_history'] = dict(max_probe_absolute_error=maximum, nonzero_probes=nonzero)
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
