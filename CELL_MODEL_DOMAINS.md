# CellModel / multi-device integration

Experimental integration of upstream `fc002359254cda9dfbb5c503cc5933306083eb57`
with the fork runtime `2c3246099763d631903a2400166e883b67350455`.
This is a separate runtime. Do not replace a frozen simulation's libraries or
plugins with this build. Rebuild CSXCAD, openEMS and GPU plugins together.

## Implemented scope

- Keep upstream's CellModel, shared point-operation kernels and single-device
  temporal blocking. Do not restore the superseded per-backend physics kernels.
- Slice the expanded cell model along X. Each device stores its owned nodes
  plus two halo planes. ADE state slots are densely renumbered for owned nodes
  only; halo histories are not allocated or advanced. Immutable coefficient,
  pole and Z-profile tables remain unchanged.
- Preserve arbitrary ordered unequal X cuts and the backend-independent
  coordinator, half-step halo exchange, pinned/borrowed host buffers, sampled
  device timers and HIP current-interior overlap.
- After **all** native extension hooks, exchange corrected tangential fields.
  A group uses one E/H step at a time. Temporal blocking and whole-step fusion
  across domains are disabled until a dependency-cone schedule is qualified.
- The eight-branch LossyMetal representation is retained. Upstream's <=4-pole
  cell fusion does not consume it. The generic fallback runs all eight poles;
  optional `OPENEMS_GPU_LOSSY_FUSION=off|apply|all` combines point operations
  only after proving identical unique positions and voltage-only R-L branches.
  This preserves the branch ordering and arithmetic, not a reduced pole fit.
- Domain point lists upload only owned excitation, lumped-RLC and dispersive
  points. Input field reads and energy reductions gather/reduce ownership once,
  excluding halo copies. Runtime gates and DFT snapshots stay at completed
  timestep/batch boundaries.

## Invocation

```
OPENEMS_DOMAIN_HALOS=overlap OPENEMS_DOMAIN_HOST_BUFFERS=borrowed \
OPENEMS_GPU_LOSSY_FUSION=all openEMS public.xml --engine=gpu \
  --gpu-devices=0,1 --gpu-domain-cuts=0,32,96
```

HIP supports `serial`, `async` and `overlap` halos. Metal logical-domain tests
use `serial`/`copy`; multiple queues on one Apple GPU are correctness tests,
not evidence of physical multi-GPU scaling. `--engine=gpu-reference` exercises
partitioning without a GPU SDK. Device indices belong to this build's device
inventory; record the actual ordered UUIDs and topology when qualifying a run.

## Explicit limits

- Cartesian, fixed X slabs; no cylindrical subgrid or repartitioning in flight.
- Mur/local absorbers, TFSF, steady-state device extensions and the experimental
  opaque-surface extension are not qualified for this domain adapter. Reject
  them instead of silently running host physics on stale device fields.
- CellModel replaces the legacy raw/auto main coefficient tables. The old
  `OPENEMS_GPU_MAIN_COEFF=raw` request is rejected. `auto` or unset uses the
  current upstream representation, which is not the previous auto-table policy.
- Group device DFT/gather remains a host-observation fallback, as before.
- A successful build or logical-Metal comparison is not physical mixed-device,
  full-array RF, memory or speed qualification. Preserve the existing production
  runtime until those controls and benchmarks pass.

## Reproducible controls

`ctest` includes cell-model slab/slot/point ownership checks and the existing
transport, field and sampled-timer controls. Public full-field comparisons:

```
python3 tests/manual/cell_model_domains.py --binary build/openEMS \
  --fixtures tests/manual/fixtures/cell-domains --output /tmp/cell-domains \
  --engine gpu-reference --basic
```

The script compares the complex E/H fields and the time probe, normalised by
the corresponding nonzero control peak (tolerance 2e-5). It covers vacuum,
PEC, lumped RLC, two-pole Lorentz and eight-branch LossyMetal with UPML, equal
and unequal cuts, and a one-owned-line slab. Repeat on actual HIP devices with
`--engine gpu --two-only`, all fusion modes, serial/async/overlap, and runtime
pause/DFT controls. Benchmark separately with timers disabled, identical
inputs/stopping steps, interleaved repetitions and explicit baseline hashes.

## Bounded Metal timing evidence

Public graded grid 193x257x97 (4,811,297 nodes), two copper bodies,
259,924 surface edges, 1536 fixed steps, one Apple GPU. Three interleaved
repetitions per mode; setup is excluded, observations are included. Temporal
blocking, whole-step fusion and timing instrumentation are disabled on both
builds. Baseline: the previous optimized fork with full eight-branch fusion.

| Variant | Median FDTD time |
| --- | ---: |
| Previous optimized runtime, full copper fusion | 3.65081 s |
| CellModel, separate copper branches | 5.22000 s |
| CellModel, apply-only copper fusion | 3.89041 s |
| CellModel, full copper fusion | 3.27232 s |

The last row uses 10.37% less time than the optimized old control (1.116x
throughput). Maximum complex-field/probe differences from that control are
6.01e-7 / 8.99e-7 of the respective control peak. This is not a prediction for
an antenna or mixed physical GPUs. Automatic temporal-block selection is a
separate comparison, not part of this table. Numerical source: `852c26e`.

```
python3 tests/manual/gpu_lossy_benchmark.py --binary NEW/bin/openEMS \
  --baseline OLD/bin/openEMS --fusion-only --baseline-fusion all --one-step \
  --shape 193,257,97 --steps 1536 --repeats 3 --output /tmp/cellmodel-timing
```

Raw per-repeat times, input/binary hashes and tolerances are preserved in
`benchmark/cellmodel-domains-metal-20261115.json`. Keep rounded native timing
precision and the synthetic workload in mind when interpreting the ratios.
