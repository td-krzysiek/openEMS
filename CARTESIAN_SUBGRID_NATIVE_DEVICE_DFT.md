# R58 — native device FieldDFT

Status: staged, unbuilt successor of immutable R57 R02. Build admission waits
for the actual R57 CPU73 receipt; physical qualification also waits for its
exact accepted physical68 terminal. No production or antenna runtime changes.

## Native contract

Use existing `Engine_Interface_Base::{Create,Accumulate,Read}FieldDFT` and
`ProcessFieldsFD` source weights, dual-time phase, frequencies and intervals.
Only matching startup `Field_Gather_Subgrid` instances opt in. Unsupported
CPU/interface/dump cases or whole-budget rejection retain host gathered DFT.

Prepare owning SI recipes and unique normalized-state support per dump. The
core R57 group stores owned selected E or H sums; the coordinator maps each
state to exactly one worker. After a successful upload mark only that gather's
dedicated recipe IDs device-only. Host observation sealing uses empty
placeholders for those IDs, so no DFT support is also read every batch. Genuine
U/I/TD supports remain active. Explicit later TD use of a device-only gather
must invoke the existing live complete-field fallback, never return zero/stale
placeholders. No field or material history is migrated.

The default shared DFT payload budget is256MiB; a C++ startup setter controls
it before group registration. Charged payload is4*N+16*N*F bytes, identical
across single/domain ownership (no per-worker weight buffers). Recipes, maps,
backend overhead, temporary host reconstruction/readback and total VRAM/RSS
are separate. Table preparation keeps the R56 budget; no global guard increase.

Register only before sealing/advancing. A sample requires sealed observations
and a complete epoch, finite correctly sized weights, and one sample/group/
epoch. Invalid requests reject before submission; partial allocation/device
failures invalidate future native observation/advance. Preallocate host registry
bookkeeping before device mutation. Shape/frequency/interface/engine/mode
identity checks precede reads. Empty regions use a nonnull complex sentinel.

At snapshot, read unique normalized sums, reconstruct SI complex raw samples,
apply physical edge lengths and round to complex float; reuse the existing
`EvaluateSums` interpolator. Check finite raw and output fields. This commutes
linear reconstruction with DFT but changes per-sample raw-float rounding, so
it is a separately qualified numerical derivative, not bit-identical legacy FD.

## Predeclared qualification

Native whole-complex-array relative maximum gate2e-5 versus direct same-build
host float accumulation with exact supplied timestep-scaled weights. No
relaxation of R56's original2e-6 gathered TD/host-DFT controls. Twelve new cases
cover all cuts, single/two workers, PML/bulk loss, nonuniform child coordinates,
drawing transforms, all interpolation modes, both fields, three frequencies,
eight batches, repeated non-resetting snapshots and two empty groups. Preserve
source/ports and existing numerical gates. Eleven cases have zero host E/H DFT
support; one retains a real over-budget E fallback. Check no full reads before
explicit fallback, no sample sum downloads and current explicit TD fallback.

Also run the unchanged Python native TD/FD/U/I controls, exact actual-parent
operator/field arrays, full CPU74, then physical72 on the authorized4090/3080
pair in both orders under existing guards. Recent Metal, matched native RF/
NF2FF, collective pause and representative PCB/antenna remain separate.


## Scoped qualification — 17 November 2026

Exact R58 passes Linux CPU74/74,16 exact actual-R57 operator/field arrays and
72 physical CUDA invocations on each RTX4090/RTX3080 and both pair orders.
Forty-eight physical native complex-array cases have maximum relative error
2.85387e-7 under the unchanged2e-5 gate. Eleven of twelve cases per lane have
no host E/H selection; empty support and live whole fallback remain covered.
Independent review rehashes794 source/439 archived files. This qualifies these
small array/lifecycle controls, not recent Metal, native RF/NF2FF, collective
pause or a representative antenna. Production and frozen prefixes are unchanged.
Project evidence: `native-device-dft-r58/independent-qualification-r01.json`.


The separate R59 R06 public1GHz PEC dipole comparison now passes with the
same native XML, one RTX4090 and a mixed RTX4090/RTX3080 pair, host versus
device DFT. Port U/I/Z/Gamma are exact; largest relative complex far-field
error is1.80653e-7, near-field2.15114e-7 and derived radiation quantities below
8.067e-8. Both NF2FF enclosures and retained public dipole/power/contact gates
pass. Single versus pair exports are exact. Independent review covers172 raw
artifacts. This is small-fixture radiation evidence, not Metal, PCB/copper RF
convergence, a speed benchmark or representative antenna acceptance.
Evidence: `native-dft-radiation-r59-r06/independent-qualification-r01.json`.
