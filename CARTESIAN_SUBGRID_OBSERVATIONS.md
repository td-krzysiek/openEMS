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
