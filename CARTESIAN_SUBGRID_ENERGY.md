# Experimental device field-energy observations

Exact R51-r02 passes scoped CPU/CUDA qualification below. Native Processing/
DFT/end-criterion behavior is unchanged; no full-field download has been removed
from native execution by this increment.

`ConfigureEnergy(grid)` is optional and startup-only. The supplied grid must be
the same unchanged topology/material normalization used to create the device
or fixed domain plan. Count checks do not establish identity with an unrelated
same-sized grid. A configuration attempt is single-use; partial backend-owned
allocations must not be retried as though freed. Matching builds are required.
Unsupported adapters explicitly reject. No fields or time steps are modified.

At a complete epoch `ReadEnergy()` returns the same positive field energies as
`CartesianSubgrid::ElectricEnergy()` and `MagneticEnergy()`:

- E: half the squared normalized electric state plus each off-diagonal
  normalized Galerkin mass product, stored once.
- H: half the squared normalized magnetic state, counting each owned face once.

This excludes the timestep-modified curl cross term, surface stored history,
resistor/bulk losses, source work and PML work. Those are separate diagnostics.
The result is not the complete modified energy/passivity ledger by itself.

Domain setup assigns every complete mass component to one owner. Energy setup
verifies both ends of every mass pair have the same unique owner; ghosts carry
no energy ownership. Workers return scalar reductions at the matching epoch,
which the coordinator sums in double. The worker-local coupling interface is
for already verified owned metadata; arbitrary incomplete supplied couplings
cannot be recognized as the actual grid's matrix.

Two shared point operations perform maximum-magnitude and quadratic reductions
on reference/HIP/Metal implementations through the existing backend dispatch.
Each field is scaled by its current maximum before squaring; the host restores
SI energy in double, avoiding square underflow/overflow for normal FP32 states.
Paired compensated sums retain low product/coefficient residuals. Nonfinite
fields, invalid pairs, incomplete epochs and negative/nonfinite energy reject.
The existing initial-state FTZ policy remains; no arbitrary subnormal-state
support or physical accuracy is inferred from this observation operation.

Only partial reduction arrays and interface mass pairs are allocated, on opt-in.
`EnergyDownloadBytes()` counts actual scalar downloads by this component:
up to24 bytes per nonzero E/H read per worker, excluding transport, processing,
other diagnostics and allocation/setup traffic. It is not a measured wall-time
speedup or a count of all transfers. Default unconfigured execution allocates
none of these buffers. Metal and physical GPU qualification remain separate.


## Scoped qualification

Fresh CPU68/68 and52 physical CUDA gates pass, with16 bit-identical actual-R49
operator/field arrays. Actual4090/3080 and both domain orders include36 targeted
energy controls plus retained loss/PML energy compositions. Source779 files and
413 archived files are independently rehashed; same-epoch energy errors retain
2e-6 FP32 and2e-4 independent-double gates. Core SHA256:
`f5b7d10c42ffbcf321a6392f039c933dbbf0f22f2d50bf55bafd4a915734c170`.

The original test configuration failure and test-only correction remain frozen.
The original external collector confused52 invocations with48 inherited native
box rows; the corrected collector changes only this row cardinality. Actual
numerical commands were not repeated and every numerical bound is unchanged.
Exact predecessor Name+ID guards passed at owner and remote admission. Metal,
native sparse processing and representative antenna acceptance remain open.
