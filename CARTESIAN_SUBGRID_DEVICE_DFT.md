# R57 — owned-state device DFT

Status: the exact R57 R02 core is qualified by CPU73/73,16 exact actual-R56
parent arrays and68 physical4090/3080 controls, then promoted on PR29. Native
FieldDFT integration is separately staged as R58. Recent Metal, native RF,
collective pause and representative antenna acceptance remain open. No
production or frozen antenna runtime changes.

## Implementation contract

- Register immutable E or H selections before advancing or sampling. Preserve
  order and duplicates; each distributed state is assigned to its unique owner,
  never a halo. Each dump owns a separate group, allowing different E/H dual-time
  weights, frequencies and sampling schedules.
- Submit finite complex FP32 weights only at complete E/H epochs, once per group
  per epoch. Registration/invalid/budget errors reject before allocations or
  dispatch. Partial worker/backend submission poisons the DFT path/coordinator;
  it must not silently retry a partially sampled epoch.
- Shared point-operation kernel, compiled through the existing reference,
  HIP/CUDA and Metal paths. No new engine, dependency or FP64 device requirement.
  Compensated FP32 real/imaginary sums include FMA product residuals. Inline
  launch weights avoid the existing blocking `DeviceUpload` contract. One kernel
  per frequency is an explicit initial implementation cost to measure later.
- No sample downloads. Only a blocking, non-resetting snapshot downloads the
  selected frequency's compensated pairs (16 bytes per selected state). It is
  taken at a complete collective epoch; all sums reflect that group's recorded
  samples, which need not include the latest completed field epoch.
- Exact device buffer payload: `4*N + 16*N*F` bytes per nonempty group. Empty
  groups allocate/launch/transfer nothing. The coordinator preflights the sum
  over all workers and groups against the supplied aggregate budget. This is
  neither total VRAM nor host RSS: backend allocation metadata, routing maps,
  host vectors and temporary snapshot arrays are separate.
- Source/fields, material histories, curl/ports, halo transport and fixed domain
  ownership remain unchanged. No device wait is added solely for accumulation.

## Required controls

Predeclared core complex-array gate: relative whole-vector maximum error2e-5
against independent double accumulation of the same-build eager FP32 fields;
exact unchanged single-device E/H arrays, original domain gates retained.
Test all three cut axes, nonuniform base/refinement, E/H, duplicate/empty IDs,
three frequencies, separately phased E/H weights, nonzero initial fields/source,
50-ohm resistor, multiple snapshots and intervening nonsampled steps. Reject
invalid/late registration, half-step sample/read, duplicate epoch, nonfinite or
wrong-size weights, invalid frequency/index and insufficient/overflow budget.
Long cancellation has an independent 99,999-sample point-kernel oracle and a
1,023-sample actual-backend control with tiny/large finite fields.

Full GPU-free suite and unchanged parent arrays precede fresh physical
4090/3080 single/pair/reversed controls under the existing resource gates.
Recent Metal qualification is separate; a Linux pass does not establish Metal.

## Remaining integration after core qualification

The native `Engine_Interface_FDTD::CreateFieldDFT` adapter must reuse these
sums and the registered R49 reconstruction recipes. It must mark device-only
recipes before sealing so R53 does not also download them every batch. Explicit
TD/U/I/full fallback must continue to return current fields. Reconstruction
and interpolation commute with the DFT algebraically, but native per-sample
raw-float rounding does not: retain this as an explicit bounded numerical
change and qualify native whole complex arrays and matched NF2FF/RF output.

Collective external pause/checkpoint integration and representative PCB/antenna
memory, radiation and mesh controls remain required. Core C++ DFT tests do not
complete native NF2FF or qualify a full antenna.

## R02 correction

R01's exact point-kernel negative reproducer fails the unchanged long
cancellation gate. Renormalize the high/low pair after each weighted sample
using a general TwoSum residual; keep the same tests/gates and all physics.
The original R01 source/build/reproducer remain frozen.

## Qualification evidence

Independent review rehashes792 source files and434 compact archive files before
exact promotion. Twenty-four physical domain/port DFT cases have maximum
complex-array error1.35356e-14, exact unchanged single-device fields and no
sample downloads. The99,999-sample cancellation and1,023-sample backend
controls pass without gate relaxation; the original R01 defect is preserved.
All587 frozen code/build files match the promoted branch. Small controls do
not establish whole-NF2FF/antenna memory, speed or RF convergence.
