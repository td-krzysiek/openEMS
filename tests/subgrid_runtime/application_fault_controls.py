"""Bounded application failure/timeout semantics on the frozen R61 executable.

Arguments: solver, guarded public protocol XML, new output directory,
cpu|reference-two|cuda-two|metal.
This driver owns only the children it starts. It never launches containers or
changes engine source, physical device selection, or production runs.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import h5py
import numpy as np


SOLVER = Path(sys.argv[1]).resolve()
XML = Path(sys.argv[2]).resolve()
OUT = Path(sys.argv[3]).resolve()
MODE = sys.argv[4]
assert MODE in ["cpu", "reference-two", "cuda-two", "metal"]
OUT.mkdir(exist_ok=False)
SOURCE = Path(os.environ["SUBGRID_TEST_SOURCE"])


def load(name):
    spec = importlib.util.spec_from_file_location(name, SOURCE / "python/openEMS" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


control = load("run_control")
startup = load("startup_gate")
snapshots = load("snapshots")
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
workers = 1 if MODE in ["cpu", "metal"] else 2
flags = ["--engine=basic" if MODE == "cpu" else "--engine=gpu-reference" if MODE == "reference-two" else "--engine=gpu", "--numThreads=1", "--exact-endcriteria"]
if workers == 2:
    flags += ["--subgrid-domain-axis=x", "--subgrid-domain-cuts=0", "--subgrid-domain-workers=" + ("reference,reference" if MODE == "reference-two" else "auto@0,auto@1")]
gate_flags = ["--startup-gate=132", "--startup-gate-id=fault-start", "--startup-gate-timeout=30", "--runtime-control", "--decay-gate-id=fault-run", "--decay-gate-timeout=30", "--dft-snapshots"]
receipt = dict(status="RUNNING", mode=MODE, workers=workers, solverSHA256=sha(SOLVER), inputSHA256=sha(XML), cases=[])


def save():
    (OUT / "fault-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


def wait(predicate, child, label, timeout=60):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        result = predicate()
        if result:
            return result
        assert child.poll() is None, (label, child.returncode)
        time.sleep(0.02)
    raise TimeoutError(label)


def paused(path, name, token):
    file = path / name
    if not file.exists():
        return None
    root = ET.parse(file).getroot()
    if root.get("state") != "paused" or root.get("id") != token:
        return None
    assert int(root.get("nativeCompletedWorkers")) == workers
    assert root.get("nativeCompletedEpoch") == root.get("timestep") == "132"
    assert root.get("converged") == "false"
    return dict(root.attrib)


def arrays(path):
    result = {}
    for file in sorted(path.glob("*.h5")):
        with h5py.File(file) as h:
            def visit(name, value):
                if isinstance(value, h5py.Dataset):
                    result[file.name + "/" + name] = value[...]
            h.visititems(visit)
    for name in ["port_ut_1", "port_it_1", "arm_lower_i", "arm_upper_i"]:
        result[name] = np.loadtxt(path / name, comments="%")
    assert result
    return result


def exact(left, right):
    assert left.keys() == right.keys()
    for name in left:
        np.testing.assert_array_equal(left[name], right[name], err_msg=name)
    return len(left)


def run(label, action=None, overrides=None, stale=None, expected=0):
    path = OUT / label
    path.mkdir()
    (path / "input.xml").write_bytes(XML.read_bytes())
    if stale:
        (path / stale).write_text("stale-token\n")
    args = [str(SOLVER), "input.xml", *flags]
    if action or stale:
        options = list(gate_flags)
        for prefix, value in (overrides or {}).items():
            options = [x for x in options if not x.startswith(prefix + "=")]
            options.append(prefix + "=" + str(value))
        args += options
    started = time.monotonic()
    detail = {}
    with (path / "solver.log").open("x") as log:
        child = subprocess.Popen(args, cwd=path, stdout=log, stderr=subprocess.STDOUT)
        try:
            if action:
                detail = action(path, child)
            code = child.wait(timeout=300)
            assert code == expected, (label, code, expected)
        finally:
            if child.poll() is None:
                child.terminate()
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
    assert sha(path / "input.xml") == receipt["inputSHA256"]
    receipt["cases"].append(dict(label=label, exit=code, seconds=time.monotonic() - started, **detail))
    save()
    return path


def startup_ack(path, child):
    return wait(lambda: paused(path, "startup-gate.xml", "fault-start"), child, "startup fence")


def runtime_ack(path, child):
    first = startup_ack(path, child)
    control.control_run(path, "pause", "fault-run-manual")
    startup.decide_startup_gate(path, "fault-start", True)
    second = wait(lambda: paused(path, "runtime-gate.xml", "fault-run-manual"), child, "runtime fence")
    return first, second


def timeout_action(path, child):
    first, second = runtime_ack(path, child)
    return dict(acknowledgements=[first, second])


def stop_action(path, child):
    first, second = runtime_ack(path, child)
    control.control_run(path, "stop", "fault-run-manual")
    return dict(acknowledgements=[first, second])


def reject_action(path, child):
    ack = startup_ack(path, child)
    startup.decide_startup_gate(path, "fault-start", False)
    return dict(acknowledgements=[ack])


def startup_timeout_action(path, child):
    return dict(acknowledgements=[startup_ack(path, child)])


def write_failure_action(path, child):
    ack = startup_ack(path, child)
    target = snapshots.request_dft_snapshot(path, "preserved")
    wait(lambda: (target / "snapshot.xml").exists(), child, "initial snapshot")
    original = {str(p.relative_to(target)): sha(p) for p in target.rglob("*") if p.is_file()}
    assert original
    # Bypass the producer's intentional used-ID guard to exercise native rejection.
    (path / "DFT_SNAPSHOT").write_text("preserved\n")
    error = target.with_name("preserved.error")
    wait(lambda: error.exists() and not (path / "DFT_SNAPSHOT.active").exists(), child, "native reused-ID error")
    assert original == {str(p.relative_to(target)): sha(p) for p in target.rglob("*") if p.is_file()}
    stage = path / "dft-snapshots/.blocked.tmp"
    stage.write_text("intentional staging obstruction\n")
    blocked = snapshots.request_dft_snapshot(path, "blocked")
    wait(lambda: blocked.with_name("blocked.error").exists() and not (path / "DFT_SNAPSHOT.active").exists(), child, "native staging error")
    assert not blocked.exists()
    assert stage.read_text() == "intentional staging obstruction\n"
    assert paused(path, "startup-gate.xml", "fault-start") == ack
    assert child.poll() is None
    # A fresh request must still work after both I/O failures, at the same epoch.
    recovered = snapshots.request_dft_snapshot(path, "recovered")
    wait(lambda: (recovered / "snapshot.xml").exists(), child, "snapshot after errors")
    left, right = {}, {}
    for root, values in [(target, left), (recovered, right)]:
        for f in root.glob("*.h5"):
            with h5py.File(f) as h:
                def visit(name, value):
                    if isinstance(value, h5py.Dataset):
                        values[f.name + "/" + name] = value[...]
                h.visititems(visit)
                assert h.attrs["sample_count"] > 0
                assert int(h.attrs["last_sample_timestep"]) <= 132
    assert left and exact(left, right) > 0
    for f in target.glob("*.h5"):
        with h5py.File(f) as x, h5py.File(recovered / f.name) as y:
            for key in ["sample_count", "last_sample_timestep", "last_sample_time", "sample_interval"]:
                np.testing.assert_array_equal(x.attrs[key], y.attrs[key])
    startup.decide_startup_gate(path, "fault-start", True)
    return dict(acknowledgements=[ack], preservedSnapshotFiles=original, errors=["reused-id", "blocked-staging"], recoveredDFTArrays=len(left))


try:
    baseline = run("baseline")
    reference = arrays(baseline)
    timed = run("runtime-timeout", timeout_action, {"--decay-gate-timeout": 1})
    witness = ET.parse(timed / "runtime-gate.xml").getroot()
    assert witness.get("state") == "timed-out" and witness.get("timestep") == "132"
    assert "controller timed out; resuming FDTD" in (timed / "solver.log").read_text()
    receipt["runtimeTimeoutExactArrays"] = exact(reference, arrays(timed))
    stopped = run("runtime-stop", stop_action)
    witness = ET.parse(stopped / "runtime-gate.xml").getroot()
    assert witness.get("state") == "stopped" and witness.get("timestep") == "132" and witness.get("converged") == "false"
    stop_log = (stopped / "solver.log").read_text()
    assert "132" in stop_log and "aborted" in stop_log.lower()
    snapshot_root = stopped / "dft-snapshots/fault-run-manual"
    assert ET.parse(snapshot_root / "snapshot.xml").getroot().get("converged") == "false"
    for hdf in snapshot_root.glob("*.h5"):
        with h5py.File(hdf) as h:
            assert not bool(h.attrs["simulation_converged"])
    for name in ["port_ut_1", "port_it_1"]:
        trace = np.loadtxt(stopped / name, comments="%")
        assert 0 < trace[-1, 0] <= 132 * 5e-12
    rejected = run("startup-reject", reject_action, expected=1)
    assert "qualification rejected" in (rejected / "solver.log").read_text()
    expired = run("startup-timeout", startup_timeout_action, {"--startup-gate-timeout": 1}, expected=1)
    assert "qualification timed out" in (expired / "solver.log").read_text()
    # Existing startup semantics leave the paused witness after the child exits.
    # Liveness is independently checked from wait()/exit above, not the XML state.
    assert ET.parse(expired / "startup-gate.xml").getroot().get("state") == "paused"
    recovered = run("snapshot-write-errors", write_failure_action)
    receipt["writeErrorResumeExactArrays"] = exact(reference, arrays(recovered))
    for stale in ["startup-gate.xml", "STARTUP_CONTINUE", "runtime-gate.xml", "RUNTIME_PAUSE"]:
        result = run("stale-" + stale.replace(".", "-"), stale=stale, expected=1)
        raw = (result / "solver.log").read_text()
        assert "stale control file" in raw and "Operator creation" not in raw
        assert not (result / "port_ut_1").exists()
    assert sha(SOLVER) == receipt["solverSHA256"]
    receipt.update(status="PASS_APPLICATION_TIMEOUT_STOP_AND_SNAPSHOT_IO_FAILURES", limits=["Small public protocol fixture, not RF qualification", "Backend failure fencing remains a separate R60 core test", "Startup witness is not a liveness certificate"])
except BaseException as error:
    receipt.update(status="FAILED_NO_RESTART", error=repr(error))
    raise
finally:
    save()
print(receipt["status"])
