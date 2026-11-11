# In-memory run gates

Opt-in run control pauses at a processing boundary, drains asynchronous writes
and flushes probes. Fields, dispersive histories and DFT sums remain resident;
continuing does not rebuild the operator or repeat a processing sample. This
is not a restart checkpoint and does not release GPU memory.

Startup qualification and runtime gates use the same native wait primitive.
Existing `--startup-gate` flags and `STARTUP_CONTINUE`/`STARTUP_REJECT` protocol
remain supported. Startup rejection/timeout still fails qualification. Runtime
timeout instead resumes the original calculation so a failed external analysis
cannot leave a long run waiting indefinitely.

## Manual control

Start with `--runtime-control --decay-gate-id=RUN_TOKEN`. Atomically publish
`RUNTIME_PAUSE` containing a unique `RUN_TOKEN-REQUEST_ID` (ASCII letters,
digits, `-` and `_`, at most 64 characters including the token). Wait for
`runtime-gate.xml` with that ID and `state="paused"`. It records the actual
step, time, last energy estimate and deadline. Publish `RUNTIME_CONTINUE` or
`RUNTIME_STOP` with the same ID. Stop uses the normal graceful finalization,
but is not energy convergence. `ABORT` remains effective while paused.

Python: `openEMS.run_control.control_run(path, action, request_id)`.
Matlab/Octave: `ControlRun(path, action, request_id)` (one controller per run).
Command-line Python arguments to `Run()` are the corresponding underscore
names, for example `runtime_control=True, decay_gate_id='run1'`.

## Automatic energy gates

```
--dft-snapshots --exact-endcriteria --runtime-control
--decay-gates=-10,-20,-30,-40,-50 --decay-gate-id=RUN_TOKEN
--decay-gate-timeout=1800
```

The list contains strictly decreasing negative integer dB values. Each level
triggers once at an energy-check boundary after all enabled Gaussian sources,
including their delays, have finished. A checkpoint ID is `RUN_TOKEN-d10`,
etc. The native gate publishes an immutable DFT snapshot before reporting
`paused`. The external controller runs NF2FF over that snapshot and continues
the same solver. NF2FF orchestration and its GPU-memory policy remain outside
openEMS. A snapshot export failure is recorded and resumes without waiting.

The ordinary end criterion takes precedence. If −50 dB ends FDTD, its normal
final field export is used instead of a duplicate −50 gate. Actual decay may
overshoot a requested level; if multiple levels were passed before the pulse
finished, their snapshots share the actual step. Metadata reports that actual
step/energy rather than inventing an earlier recording. This does not change
the numerical end criterion or establish mesh/radiation convergence.

Use a fresh run directory. Decisions must match the active ID; mismatched IDs
fail rather than resume another gate. Gate timeout is 1–3600 seconds. Runtime
wait/snapshot time is excluded from native stepping speed and iteration time.
An external job must have its own shorter hard timeout: pausing alone does not
terminate an orphaned NF2FF process. The controller should retain checkpoint
results independently and keep its GPU lease while admitting a child process.

Run `tests/manual/runtime_gates.py` with `--binary`, `--engine` and `--output`
to compare repeated gates, manual/startup continuation and timeout with an
uninterrupted run; it also checks stop, abort and mismatched-token behavior.
