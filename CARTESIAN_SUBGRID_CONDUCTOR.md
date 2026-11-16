# Cartesian subgrid native smooth conductors and resistive paths

Experimental extension of the passive surface component in
[CARTESIAN_SUBGRID_SURFACE.md](CARTESIAN_SUBGRID_SURFACE.md). It does not qualify
production PCB antennas, arbitrary CAD or distributed lossy execution.

## Native conductor mapping

With a compatible CSXCAD `CSPropLossyMetal` API, untransformed positive-volume
boxes with finite positive conductivity and Rq=0 can be mapped to opaque local
cells. Physical height remains geometric: 18/35 µm is not substituted by a
zero-height sheet. Interior E/H, curl and mass contributions are excluded;
exposed top/bottom and sidewall faces supply the surface quadrature.

Each resolved normal column is walked through actual occupied cells. Contiguous
same-conductivity runs determine the fit thickness. Splitting, overlapping or
reordering boxes therefore preserves the physical union and its constitutive
fit. Material priorities are applied before occupancy/thickness inference.
Touching different-conductivity bodies explicitly reject; they require a
separate interface model. Box boundaries and surface traces must resolve on
the supported mesh. Unresolved coarse/fine faces reject rather than disappear.

The public GPL smooth-slab fitter supplies a positive G+RL admittance. Local
edge descriptors avoid a global voltage-path search for every surface trace.
The paired conductivity extension owns the midpoint advance, surface histories,
Ohmic losses and stored energy. The fit band is captured before excitation
preparation resets its temporary frequency metadata. This is a smooth exposed-
face approximation, not a general conductor-diffusion or roughness model.

## Native resistive ports

The native adapter accepts complete axis-aligned paths of local edges. A CSX
resistor spanning N series cells maps to N resistors of R/N, matching the
ordinary solver's line-box convention. The explicit C++ whole-functional
`Grid::SetResistors` API remains a separate circuit description.

Partial/off-edge paths, finite transverse boxes, overlapping separate native
resistors and reactive native elements remain unsupported. The existing
64-circuit guard applies to resolved segments. Physical voltage-path/current-
contour observations must be matched before comparing impedances.

## Dependency compatibility

CSXCAD is independently versioned. Builds without the optional volume-
LossyMetal header retain supported material/PEC/port behavior. Only its two
native volume-metal feature tests return 77 and are marked skipped. A header
without the required conductivity/Rq/roughness metadata is not qualification.
Native conductor use requires the matching API; there is no PEC fallback.

## Evidence and candidate status

Candidate r02 passes the complete 34-stage CPU/reference/physical-Metal matrix,
CTest 37/37 and Python 14/14. Separate physical RTX4090 and RTX3080 lanes pass all 27
CUDA gates and CTest 37/37, with actual ordered-device/runtime identities retained.
11,255 prior numerical arrays are bit-identical. Two-device controls cover the
retained lossless field/port lane; they do not distribute conductor histories.

An independent per-cell circuit detects the initial unpublished whole-path
native-R mismatch. A separate negative control reproduces the earlier
first-primitive-thickness error for an unchanged same-metal union. Both are
experimental-adapter defects, not ordinary upstream-solver defects.

Short physical-height 3/4 mm microstrip fixtures, excited from both ends, compare
S=B A^-1 with an independently qualified ordinary smooth-conductor solver.
For the scoped ratio2 mesh, maximum band complex-S difference is 0.002025;
transmission phase difference is 0.0745–0.1143°. The incremental Cu/PEC correction
agrees within 1.92e-5 in complex ratio. This is not characteristic-Z, measured
loss, full mesh convergence or antenna qualification. Current-plane alignment
and the native R/N correction are necessary parts of this comparison.

The separate closed parallel-plate TEM pilot compares conductor-induced decay
and frequency shift against an analytic smooth-slab root: sampled decay error
is below 0.39% and shift error below 0.22% for 24/48 longitudinal cells. These are
scoped analytic checks, not proof of general interface or port accuracy.

A real older-CSXCAD build passes its 3/3 tests and 35 solver tests, with exactly
two explicit volume-metal skips. The ordinary material/PEC/native path-port
runtime links the actual older library. Python ABI compatibility and volume
metal on that absent-feature dependency are not claimed.

The current r03 source adds the optional-dependency guard to that numerical
candidate. Its modern-dependency CTest 37/37 and focused physical Metal controls
pass; the complete 34-stage integration replay and separate physical CUDA replay
are still running at this publication. Do not treat the r02 receipt as their
terminal result. Local r02 core SHA256 is
`7a7d8b5e94c2624f829e2f46592ce44ebcf0722d553c79013ba494b310a9f167`;
its CUDA core is
`a8d2815b012cf64e9e66c0a103c7ce6207d791d6d256570632e94d24755f0240`.

A separate fixed-current-plane two-resolution diagnostic retains within-level
native/ordinary S differences below 0.00223, but the mesh change is approximately
0.099 in BOTH solvers. The point-port fixture is not a converged RF standard.
It cannot establish antenna accuracy or characteristic impedance. A two-length
attempt also produces negative apparent attenuation in both lanes and is not
admitted as a physical conductor-loss reference.

Further acceptance requires spatial/time convergence, broad conductor CAD and
ports, scalable field/interface representation, distributed loss/PML/DFT state
and representative antenna comparisons. F02 remains incomplete.
