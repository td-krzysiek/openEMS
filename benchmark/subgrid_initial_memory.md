# Experimental subgrid initial-state host memory

`subgrid_initial_memory.cpp` is the exact driver compiled for the R55 matched
host controls. It constructs homogeneous zero fields, consumes four streamed
worker packages, records setup peak RSS and then performs a CPU step that
materializes dense fields. It samples64 curl rows and12 field positions and
checks ownership/counts and zero energy. No GPU is used.

Compile the same driver separately against each matching source/install pair
(CSXCAD and other normal build dependencies must be installed). On a POSIX
host, for example:

```sh
c++ -std=c++17 -O3 -ffp-contract=off -I"$SUBGRID_SOLVER_SRC" \
  benchmark/subgrid_initial_memory.cpp -L"$SUBGRID_SOLVER_PREFIX/lib" \
  -Wl,-rpath,"$SUBGRID_SOLVER_PREFIX/lib" -lopenEMS \
  -o /tmp/subgrid-initial-memory
/tmp/subgrid-initial-memory patterns 128
/tmp/subgrid-initial-memory blocks 128
```

Use a fresh process for each case: `ru_maxrss` is a process high-water mark.
The driver explicitly opts into a16M-leaf count budget; the solver's public
default remains1M. The128 case has3,932,160 leaves.160/200 are larger optional
characterizations and need separately checked RAM/time limits. The driver does
not enforce an external process/container budget. Do not infer antenna/device
capacity from the count setting or from this homogeneous geometry.

Matched Linux controls used4CPU/4GiB containers, unchanged dependencies and
actual parent `d2d6f20870ed84a31b64a1a644e95b2545d034e5` versus the implicit-zero
source promoted in `4b6a31a0346c2416d8fc8c84eef74879cb015da4`. Each case was a
separate process. Both layouts save281,743,856 accounted initial bytes
(268.69MiB). After the CPU write/step, allocated field storage is comparable.

| Layout | Parent plan peak MiB | Implicit-zero plan peak MiB | Reduction |
|---|---:|---:|---:|
| patterns | 1409.14 | 1140.34 | 19.08% |
| blocks | 1329.76 | 1060.47 | 20.25% |

These are single host observations, not repeated statistical estimates or
FDTD-throughput/GPU/RF measurements. Compare exact counts, ownership, sampled
curls/fields and energy before interpreting memory. Native sampled PCB, device
observations and complete antenna resource qualification remain separate.
