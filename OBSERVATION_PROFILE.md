# Optional field-observation timing

Set `OPENEMS_OBSERVATION_PROFILE` to a unique per-process JSONL path before
starting the solver. It is off by default. Existing XML, DFT sampling, fields,
accumulation order, kernels and checkpoint schedules are unchanged. This does
not enable device DFT, surface-only gather or UPML fusion for grouped GPUs.

```sh
OPENEMS_OBSERVATION_PROFILE="$PWD/observation-profile.jsonl" openEMS model.xml --engine=gpu
```

Python: set `os.environ['OPENEMS_OBSERVATION_PROFILE']` before `Run()`.
Matlab/Octave: use `setenv('OPENEMS_OBSERVATION_PROFILE', filename)` before
`RunOpenEMS`. For containers, pass the variable and a writable mounted path.
The environment is read once, when the first processing object is constructed;
it cannot instrument an already running executable. CPU and Metal can use the
same option. A shared-memory Metal mirror does not generate transfer records.

The schema-1 header identifies host steady-clock wall durations, followed by
scope rows with `stage`, `owner`, `step`, `host_wall_seconds`, `logical_bytes`,
`sample_interval` and `frequencies`. Records are buffered, mutex-protected for
background dumps and flushed at processing boundaries/background-sample end.
Use a local filesystem. Appended headers separate successive process runs;
do not let multiple processes write the same file. Preserve XML/runtime hashes,
engine/device identities and run timestamps alongside the profile.

| Stage | What is timed |
| --- | --- |
| `gpu_full_download_E/H` | Existing full host-mirror download call, including existing completion waits; bytes are the logical copied field size. |
| `dft_extract` | Surface interpolation/gather on the host; can include the full downloads above. |
| `dft_cpu_accumulate` | CPU complex multiply/add over one extracted surface and its configured frequencies; bytes describe the input, not PCIe traffic. |
| `dft_snapshot_wait` | The existing wait for an asynchronous field snapshot; no extra wait is inserted. |
| `dft_device_submit` | Host submission of a device DFT update, not GPU kernel completion time. |

Durations are **not additive** across nested scopes or asynchronous workers.
In particular, do not add full-download time again to `dft_extract`. Separate
E/H mirrors are cached: twelve surface dumps do not imply twelve full copies.
Full mirror downloads triggered by other observations are recorded too; do not
attribute every transfer to DFT. Range reads for scalar probes and halo DMA are
not measured here. Host-fallback extension transfers and backend-internal
snapshot DMA are outside the full-download counters. This is not a PCIe
bandwidth meter or GPU idle/device-clock measurement. Standalone NF2FF surface
integration runs in another process and needs its own existing timing log.

There are no added device waits, events or field copies. Enabled monitoring
does read host clocks and format records; measure its overhead on a suitable
public observation-heavy fixture before using it for performance claims. The
disabled scopes read no clock, copy no names and write no profile. An output
open failure warns and disables monitoring; later I/O failure warns without
aborting the solver, and can leave incomplete records. Ignore a partial final
line after a crash.
