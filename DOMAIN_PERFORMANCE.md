# Sampled domain performance records

This opt-in observer records a fixed multi-GPU run for an **external planner of
the next run**. It never moves fields, changes cuts, or chooses resources.

```sh
OPENEMS_DOMAIN_PROFILE=/new/run/domain-performance.jsonl \
  openEMS model.xml --engine=gpu --gpu-devices=0,1 --gpu-domain-cuts=0,100,300
```

Python can set `os.environ['OPENEMS_DOMAIN_PROFILE']` before `FDTD.Run(...)`;
Matlab/Octave can use `setenv('OPENEMS_DOMAIN_PROFILE', path)` before
`RunOpenEMS(...)`. Clear the variable afterwards when it should not apply to
another run. The file must not already exist. A failed profile open/write is
reported and disables recording; it is not a reason to discard a valid solver
result. A writer-close receipt does not certify solver success or convergence.

Defaults: skip the first 4096 full E/H steps, then sample one full step every
4096 steps. Optional `OPENEMS_DOMAIN_PROFILE_INTERVAL` (minimum 128) and
`OPENEMS_DOMAIN_PROFILE_WARMUP` override that cadence. There is no full-grid
copy. Timers are reused, file output is buffered, and the buffer flushes every
eight samples and at destruction. A crash can lose the last buffered records;
complete records already written remain useful. A short run may have no samples.

The HIP adapter records device-stream events around main updates and individual
native extension hooks. Results are queried after the **existing** step
synchronizations, with no profiling event synchronization. An unready timer is
reported as unavailable, never waited for. Event intervals may contain launch
gaps within a hook: they are not hardware utilization counters. Profiling adds
some event/CPU work; benchmark its overhead before selecting a production cadence.

Backends without a nonblocking device timer, including the current Metal and
CPU-reference adapters, still report host durations but have null device times.
Do not treat missing device timing as a fast or idle device. The older
`OPENEMS_HIP_KERNEL_TIMES` diagnostic waits after every kernel and must remain
disabled when measuring normal multi-device throughput.

## JSONL contract (`openems.domain-performance.v1`)

- `header`: solver version, exact grid, ordered domain IDs/backend names, owned
  half-open X intervals/node counts, warmup and cadence. Index order is the
  order of `--gpu-devices`, not necessarily physical host GPU numbering.
- `sample`: zero-based full-step index and wall-clock timestamp; each domain's summed main
  and native-extension device intervals, interval counts/completeness, separate
  host submission/blocking durations, pack/unpack wall durations and sent/received
  halo byte counts. `exchange_wall_seconds` includes the coordinator's transfer
  work and waits; it overlaps per-domain pack/unpack metrics, so do not add them.
- `closed`: count of serialized samples, written on normal backend destruction
  including error unwinding; it says only that the writer was closed.

Host blocking in the sequential completion loop is **not GPU idle time** and is
order-dependent. Use device times to assess compute cost; inspect communication
and host waits separately. Device scopes exclude setup, host-fallback extensions,
processing/DFT, energy checks, NF2FF and most physical-boundary handling. That last
host interval is reported for the first domain separately. These metrics are not
an additive accounting of complete wall time.

The launch controller must bind this stream to input/library/plugin hashes,
ordered physical UUIDs, device memory budgets, topology/driver, material model
and mesh identity in its execution metadata. Backend names alone are not stable
device identity. An external reader can normalize sampled device seconds by
owned nodes, but different PML/material densities, surface histories and halo
costs make that an approximation, not an automatic optimal split.

The next-run planner owns compatibility checks, memory constraints and any new
partition choice. Record a proposed split as a new run's execution configuration
and qualify it independently. The solver never learns or repartitions silently.

## Per-extension attribution

Set `OPENEMS_DOMAIN_PROFILE_EXTENSIONS=1` alongside `OPENEMS_DOMAIN_PROFILE`
to label the **existing** sampled event pairs by extension instance and hook.
This adds no timers, GPU events, device waits, field transfers or kernel changes.
Without this opt-in the aggregate records remain available as before. Names
are registered once at setup; IDs distinguish repeated instances such as PML
regions. The sample adds `extensions` per domain, containing instance totals
and the seven pre/post/apply E/H and synchronize hooks. A hook may be a no-op
on a domain, so small event intervals include dispatch/launch gaps; they are
not an isolated kernel benchmark. No allocation is added in unsampled steps.

Each entry reports calls, missing event intervals, device seconds (null if
unavailable) and separate host submission seconds. Missing measurements must
not rank as zero/fast. Existing event pairs, warmup/cadence and buffered writing
are reused. Group totals must reproduce the aggregate extension time. Added
CPU aggregation/JSON output still has a cost: measure aggregate versus detailed
logging at the same cadence, not profiling-off versus many more GPU events.
Family totals should be formed within each sample before taking medians, and
reported independently for each ordered physical GPU. These are not percentages
of full application wall time. The synchronized per-kernel diagnostic remains
unnecessary for this low-overhead attribution.

This observer cannot be attached to an existing executable; enable it on a
subsequent qualified launch. It does not enable UPML fusion or any other solver
optimization. The current Metal/reference adapters still expose host times only.

This feature adds an optional virtual backend method: rebuild the core, plugins
and any C++ dependants together. Do not load an old plugin into this build.
