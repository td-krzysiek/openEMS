# Startup qualification without rebuilding the operator

Generated with assistance from OpenAI Codex.

Opt in with `--startup-gate=4096 --startup-gate-id=UNIQUE_ID`. Run the final
XML from the beginning, including its final timestep cap and decay criterion.
The solver pauses once at the first existing processing boundary at or after
4096 steps. It does not change the stepping/sampling schedule or reset fields,
material histories, maximum energy or DFT sums. The operator is built once.

`startup-gate.xml` is published atomically after asynchronous processing has
finished and existing probe data are available. It records the ID, actual
timestep, timestep duration and `converged="false"`. A controller must inspect
startup logs, resource usage, probes and any required snapshots independently;
the witness alone does not qualify a run.

Publish a `STARTUP_CONTINUE` file atomically, containing the matching ID plus
an optional newline, to continue the same process. `STARTUP_REJECT` with that
ID fails the run. Missing acceptance times out after300seconds; set
`--startup-gate-timeout=SECONDS` between1and3600. Existing control files,
invalid IDs, IDs from another request and an unreachable configured cap are
errors. A fresh run directory and unique ID are required. This is a local
single-controller protocol, not a network control interface.

`ABORT` and SIGINT remain responsive during the pause and retain the normal
graceful-stop behavior. With `--dft-snapshots`, the controller can request and
read an immutable DFT snapshot while paused, then accept or reject. No ordinary
final field dump is labelled complete by the gate. A paused run retains its
RAM/VRAM; resource ownership must remain held until continuation or exit.

The qualification wait is excluded from reported FDTD throughput/time and the
wall-clock end-criterion mode switch. It remains part of controller wall time.
For reproducible end steps use `--exact-endcriteria`. No gate flags preserves
the ordinary run path. The protocol is shared by CPU, GPU reference and device
backends; it is not a disk checkpoint or restart facility.

Python forwards these command-line options through
`FDTD.Run(path, startup_gate=4096, startup_gate_id='ID')`. Matlab users can pass
native flags via `RunOpenEMS(path, xml, options)`. Rebuild bindings against the
updated library; do not replace a library beneath old compiled bindings.
No XML material/geometry attributes are added by this feature.

Native regression control:

```sh
python tests/manual/startup_gate.py --binary /path/to/openEMS \
    --engine basic --output /new/result/directory
```

Repeat with GPU reference and supported actual GPU backends. The control checks
final nonzero complex DFT data/probes and exact decay-stop steps, paused output
stability, live snapshots, rejection, timeout, invalid/stale IDs and graceful
abort. Numerical comparisons declare tolerances and report actual errors.
