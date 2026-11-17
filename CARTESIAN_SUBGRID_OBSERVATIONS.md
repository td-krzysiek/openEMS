# Experimental selective normalized-state observations

`GPU_CartesianSubgrid::ConfigureObservations(electricIDs, magneticIDs)` defines
one ordered selection at unused startup. IDs address normalized compact E/H
states, not physical point-field values. Duplicates and caller order are
retained. Invalid IDs reject before changing the configuration. Reconfiguration,
configuration after stepping and observation between electric/magnetic half
steps reject. Output vectors must be distinct. An empty selection is valid.

`ReadObservations(electric, magnetic)` returns only the selected FP32 states at
one completed epoch. The existing backend indexed-copy kernel performs the
gather; no new backend or CPU stepping fallback is introduced. Metadata/device
buffers are prepared once. Local workers read only owned IDs; the coordinator
resolves global IDs to the unique owner and restores requested output order.
Unsupported process/transport adapters explicitly reject this operation.

Both single-device and coordinator complete-field host mirrors are lazy:
`ReadFields` preserves the prior complete observation contract and allocates its
legacy buffers on first use. `HostObservationStorageBytes` accounts only for
owned host observation buffers/selection destinations. It excludes global-grid
storage, device allocations, caller outputs and general process overhead; it is
not an estimate of complete simulation RAM.

This C++ API is a building block for planned native probe/DFT gathering. The
native engine currently still calls complete `ReadFields` at processing batch
boundaries. It therefore does **not** yet avoid complete native field mirrors
or establish production transfer/RAM savings. Native field interpolation,
weighted U/I functionals, energy reduction and DFT integration remain separate.
Normalized-state samples must not be presented as physical E/H without applying
the correct reconstruction weights.

No time-step arithmetic, material histories, partition, halo exchange or
observation-time definition is changed. Matching callers/backend plugins must
be rebuilt. This candidate requires its own CPU, Metal and physical-device
qualification; parent results alone do not qualify the new interface.


## Scoped qualification

The separately rebuilt R48-r04 passes66/66 Linux CPU tests and16 bit-identical
actual-R46 operator/field files. All44 physical CUDA gates pass on RTX4090 and
RTX3080 individually and in both domain orders, including24 selective-state
controls, retained native TD/FD/ports/PML and sampled bulk losses. Exact
same-epoch selective/complete downloads agree. Device/source/runtime identities
and401 archived files were independently rechecked. Qualified core SHA256:
`40d3b111e04ae7208c682aa0a867c424b6f130840f653791f1fc8bd2e8ead59f`.

The first test's independent FP32 helper lacked resistor support and failed
also against unchanged R46. The retained correction uses independent double
reference for port cases under the established2e-4 bound; source-only and
single/domain FP32 gates remain2e-6. Original failure/diagnosis remain evidence,
not a numerical code fix or relaxed solver threshold. Metal is still open;
no native DFT/antenna/performance acceptance follows from these small controls.


## Physical E/H/U/I observations — R49

`ElectricSample` and `MagneticSample` build owning linear SI reconstruction
stencils, preserving local interpolation, hanging traces, electric mass factors
and leaf permeability. `CartesianSubgridObservations` accepts these and existing
voltage-path/current-contour functionals, deduplicates state IDs at startup,
and reconstructs ordered outputs from selected complete-epoch readbacks.
The device must originate from the same unchanged grid and outlive the reader;
count checks alone do not identify topology. No native Processing/energy/DFT
readback behavior changes in this increment.

Exact R49 source passes67/67 CPU tests,16 bit-identical frozen-R48 files and48
physical CUDA gates on4090/3080 and both pair orders. The24 new physical controls
include105 functionals/case, analytic SI normalization, opaque/material
interfaces, direct reconstruction and independent double stepping. Maximum
relative double-reference error is3.23037e-7; selected/complete and single/domain
outputs agree exactly. Qualified core SHA256:
`4f0d52ceb9c75d15be8d9ef4c404d77f7b8b16429af9a611ab00156f0f4df948`.
Source and406 archived files were independently rehashed. Metal remains open.

The external physical-test owner checked an older completed job instead of its
intended RF predecessor. This harness defect allowed overlap with our own RF
series; resource limits still passed and numerical cases used independent
state. Preserve the original launch/audit and use no timing-isolation or
performance claim from these controls. The next owner requires exact container
name AND ID terminal checks twice, with positive/negative regression coverage.

## R53 staged native observation integration

Unqualified native processing registration uses `RegisterObservation` before
stepping and `SealObservations` after processing setup. `EvaluateObservation`
returns registered U/I (or explicit C++ SI E/H functionals) at the complete
processing epoch. The native application seals automatically after setup;
manual C++ engine users remain on the existing eager path until explicit seal.
Registration, sealing and sample IDs are checked. Native CPU remains direct.

Sealed device execution gathers only registered states and scalar work/loss at
batch boundaries. `CalcFastEnergy` reads optional device field energy once per
epoch and retains separately mirrored surface-history energy, preserving the
existing vacuum scaling of the engine interface. GPU reduction rounding may
change a threshold-edge termination batch; field stepping is unchanged.

`ObservationGrid()` explicitly refreshes a complete host mirror on demand.
Generic E/H field queries/dumps still use that full fallback, serialized once
per epoch before concurrent read-only access. Raw `Grid()` is setup/reference
storage, not a live device mirror after sealing; previously held raw references
must not be used as live fields. Sparse diagnostic calls do not refresh it.
Full device mirror removal, planned E/H/DFT gathers and device DFT accumulation
are subsequent work, not claimed by this increment. Complete initial double
host fields are still allocated. Counters describe this component only and
exclude setup, halos, other extension transfers and total RSS/VRAM.

The new control compares eager/sealed single/two workers, three cut axes,
PML with/without bulk loss, native U/I and SI functionals, energy/epoch caching,
unchanged complete fields, invalid operations and concurrent full fallback.
Representative RF/antenna, Metal and actual physical CUDA remain separate gates.


R53 qualification: CPU69/69,16 bit-identical actual-R51 arrays and56 physical
CUDA gates pass on RTX4090/3080 and both pair orders. All780 source and418
archived files were independently rehashed. The48 new physical native cases
retain exact selected/full functional values, maximum relative energy error
7.60313e-9, zero full reads over four batches and one full refresh per observed
epoch under concurrent readers. Exact source is promoted; Metal, prepared
E/H/DFT, PCB scaling and representative RF/antenna remain open. Evidence:
`native-sparse-observations-r53/independent-qualification-r01.json`.


## Qualified material/observation composition — R54

R54 composes the exact bounded material palette with native selected probes
and scalar energy. Exact repeated epsilon/mu retain their original bits;
homogeneous media need no per-leaf IDs, up to4096 pairs use32-bit IDs, and
fully varying media fall back to dense storage. Setup dictionary memory is
released before field construction. CPU70/70,16 actual-R53 operator/field
files and56 physical CUDA gates pass. Independently rehashed782 source and419
archive files. Exact numerical/API/test source is promoted.

R52's3.54–3.68% peak-RSS measurements remain separate homogeneous host evidence;
they are not a native PCB or full-antenna capacity measurement for R54. Initial
double fields and generic full fallback remain in this increment. R55 stages
implicit zero fields separately; Metal, native PCB scaling, E/H/DFT and
representative antenna acceptance remain open.
