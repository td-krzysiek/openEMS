# Experimental CPU LossyMetal setup

`OPENEMS_LOSSY_SETUP` selects surface scanning during operator preparation:

| Value | Scan workers | Reuse of metal extents |
|---|---|---|
| `serial` (default) | 1 | off |
| `cache` | 1 | on |
| `parallel-no-cache` | setup thread count | off |
| `parallel` | setup thread count | on |

Use the existing `--numThreads=N` to bound setup parallelism. Invalid modes
fail explicitly. This is a host setup optimization, independent of the later
CPU/GPU engine. It does not change the mesh, surface-impedance equations,
roughness parameters, fitting accuracy or time-stepping kernels.

Each worker scans consecutive X columns with private candidate lists, metal
boundary caches, contributions, fit maps and counters. After joining, results
are merged in original node/first-fit order. Primitive-used flags and operator
coefficient mutations remain serial; the compressed coefficient override map
is not a concurrent container. A fit may be calculated once per worker rather
than once globally. In particular, nonzero-roughness fitting is not free.

An extent entry stores a known interval of metal grid nodes, the last inside
node and the **unchanged** boundary-bisection tail. Reuse evaluates the original
coordinate difference plus that tail; it does not approximate a contour or
assume a uniform grid. Keys include the edge direction and both fixed line
coordinates. Positive/negative searches and disconnected bodies stay separate.
The cache only knows intervals actually traversed by the original search.

Cached modes retain at most 4096 candidate lines per worker and 64 extent
records per sign/line before clearing and recomputing. These are entry limits,
not a byte budget: candidate lists have geometry-dependent lengths. Uncached
modes retain the original unbounded candidate-line behavior. Contributions,
positions and duplicate worker fits also consume host memory. Surface scanning
is the parallel part; ADE assembly and the interior magnetic-edge mask remain
serial. No speedup proportional to thread count is promised.

`OPENEMS_PROFILE_LOSSY_SETUP=1` still emits one combined record. Phase times
are wall time. Query counts and sampled durations are summed across private
workers: inclusive estimates represent worker time, **not** parallel wall time,
and nested query families must not be added together.

## Regression

Build with `OPENEMS_BUILD_TESTS=ON`, run CTest, then:

```sh
python3 tests/check_lossymetal_setup.py build/tests/inspect_lossymetal_setup /tmp/lossy-state-check
```

The output directory must be new. An optional third argument selects an
inspector built against the unchanged parent. The helper reads protected state
only from its own extension instance. It compares all LossyMetal ADE positions,
internal/external coefficients, flags, surface counters and fitted G/a values,
plus native base VV/VI/II/IV and PEC dumps, for 1/2/4-worker cached and uncached
modes. Synthetic controls cover nonuniform grids, disconnected/boundary metal,
curves, priority air voids, subcell and zero-height sheets, missing-thickness
fallback, two metal properties, Rq>0 and a PEC-only control. No FDTD steps run.
Byte equality is a strict optimization regression on the same build/host, not
a requirement of binary equality between different numerical backends.

Long-run RF, CUDA/HIP execution, sanitizer race checks and full-antenna setup
scaling remain separate qualification. Do not replace a frozen production
runtime solely because these local setup controls pass.

## Local measurements

Apple M5 Pro, four setup threads, three interleaved repeats, Release build.
Values are median native operator-only setup wall time; no FDTD steps run in
these approximately one-million-location cases.

| Fixture | Parent | Serial | Cache | Parallel, no cache | Parallel + cache |
|---|---:|---:|---:|---:|---:|
| Public 1,044,480-location cavity |0.535s|0.527s|0.500s|0.245s|0.254s|
| Local 1,000,500-location diagnostic |12.863s|12.902s|10.754s|5.096s|4.297s|

The local diagnostic's LossyMetal stage decreases from11.236s to2.834s;
separately measured process peak RSS is213.3MiB versus189.3MiB. Neither the
private geometry nor its outputs are distributed here. The default serial
refactor differs by0.31% on that fixture. Cache bookkeeping is not always a
benefit: the simple cavity is slightly faster with parallel scanning alone.
These are local fixture results, not a full-model scaling forecast.

The standard benchmark adapter can compare two builds and their time stepping:

```sh
OPENEMS_LOSSY_SETUP=parallel python3 benchmark/lossymetal_cpu.py run \
  --openems /path/to/parent/openEMS --label parent \
  --openems /path/to/candidate/openEMS --label candidate \
  --cases ideal,lossymetal --size XS --engine basic,cpu,gpu \
  --threads 4 --duration 0.4 --repeat 3 --workdir /tmp/lossy-speed
python3 benchmark/lossymetal_cpu.py compare /tmp/lossy-speed/result.json
```

The unchanged parent ignores the new environment selector. All six local
basic/CPU/Metal probe comparisons were identical; throughput changes were
within measured noise. Separate4096-step nonzero-probe controls also matched
across setup modes within each engine. This does not qualify HIP/CUDA or a
long-running antenna calculation.
