# Cartesian subgrid native smooth conductors and resistive paths

Experimental extension of the passive surface component in
[CARTESIAN_SUBGRID_SURFACE.md](CARTESIAN_SUBGRID_SURFACE.md). It does not qualify
production PCB antennas, arbitrary CAD or distributed lossy execution.

## Native conductor mapping

With a compatible CSXCAD `CSPropLossyMetal` API, untransformed positive-volume
boxes with finite positive conductivity and Rq=0 can be mapped to opaque local
cells. The separately qualified shaped increment below extends this to simple
linearly extruded polygons and axial cylinders/hollow shells. Physical height remains geometric: 18/35 µm is not substituted by a
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
LossyMetal header retain supported material/PEC/port behavior. Its native volume-metal feature tests return77 and are marked skipped
(three tests in the current inventory). A header
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
candidate. All 34 local stages pass: CTest 37/37, Python 14/14 and 11,255 prior
arrays bit-identical. All 27 physical CUDA gates and CTest 37/37 pass, including
native conductor/path controls on each RTX4090/RTX3080 and the retained lossless
pair controls in both orders. The source inventory is identical to fork commit
`c006caa9e1da257175a8798b0fc80f0fb1f5dbab` for every numerical/test source.

Actual r03 local core SHA256:
`8c306301e648183b29f76486cf718f8ab09d087723e27bb78bf9e17da5b8f120`.
Actual r03 CUDA core SHA256:
`584f586e463c87911fb628adb7a94356836f062afb25b064325d3c121b5d9fa4`.

A separate fixed-current-plane two-resolution diagnostic retains within-level
native/ordinary S differences below 0.00223, but the mesh change is approximately
0.099 in BOTH solvers. The point-port fixture is not a converged RF standard.
It cannot establish antenna accuracy or characteristic impedance. A two-length
attempt also produces negative apparent attenuation in both lanes and is not
admitted as a physical conductor-loss reference.

A third-level fixed-interior-voltage diagnostic completes local/ordinary
ratio2/4/8 cases with PEC and smooth copper. The maximum band phase change
falls from approximately12.9rad/m (2→4) to5.7rad/m (4→8) in both solvers;
within-level phase differences remain below0.637rad/m. The incremental conductor
alpha/beta differences are below0.825%/0.562%. Ratio8 uses5fs, versus10fs at
ratio2/4, so this is combined spatial/time evidence. Five fixed voltage paths
avoid the changing point-current normalization, but the finite-ground,
point-line launch and single-mode fit remain diagnostic, not modal-port,
physical attenuation or antenna acceptance. A separately frozen r04 enlarged-air-collar diagnostic completes both
PEC/copper cases; changing only that collar shifts fitted beta by less than
0.001rad/m. It does not resolve the finite-fixture/modal-port qualification gap.

## Native shaped conductors

The r08 adapter accepts untransformed simple linearly extruded polygons and
axis-aligned cylinders/hollow shells with physical volume. It reuses native
CSXCAD membership and priorities, centre-classified opaque leaves, exposed
Cartesian faces and the existing passive surface advance. The supplied CAD
coordinates are unchanged. Invalid/self-crossing/backtracking contours,
zero-height smooth sheets, transforms and unsupported axial contacts reject.
There is no new fit, numerical backend, reference-plane shift or roughness model.
The current r08 guard requires omission of a repeated final closing vertex; a
reproduced equivalent closed-ring rejection is an experimental-adapter
restriction under follow-up, not an ordinary upstream solver defect.
Direct C++ callers must update CSXCAD before operator setup, as native XML/
Python setup already does. Holes need supported shell/priority geometry or
simple decompositions, not self-touching polygon bridges.

`tests/test_cartesian_subgrid_shaped_surface.cpp` compares concave and rounded
contours, solid/hollow Z cylinders, polygon normals X/Y/Z, Cu18/Cu35, reversed
extrusion/winding, substrate loss and alternative units/origin. Every native
mask matches an independent analytic centre classifier; native coordinates
remain exact. CPU fields, loss and occupied volume match independent occupied-
box decompositions exactly. All hollow-bore, passive-energy and five invalid-
input controls pass. Maximum normalized shaped-field errors are
9.71489e-07 on physical Metal and 9.71489e-07 across physical RTX4090/RTX3080,
within the unchanged 2e-4 gate.

The composed r08 replay passes 15 local stages, CTest40/40, Python14/14 and
11,255 unchanged prior numerical arrays. Physical CUDA passes 35 gates and
CTest40/40, including each device and retained lossless pair controls in both
orders; terminal state is exit0/no OOM with source/runtime rechecked. A separate
build against the actual older CSXCAD library passes three scoped native
controls and explicitly skips the three missing-volume-feature tests. This
six-test control is not a full old-dependency or old-Python-ABI replay.

This is Cartesian sampling, **not conformal curved-surface qualification**.
A local isolated rounded-cap diagnostic improves footprint area error from
about 8.5% at 50–100µm to 0.68% at 12.5µm. Its stair sidewall perimeter nevertheless
remains about 24% above the smooth contour. Purely geometric surface area is
not RF loss: field distribution and quadrature also matter. Do not infer
conductor-loss/phase accuracy from footprint convergence; a matched curved-
line RF comparison is required before antenna adoption. Native conductor/PML
composition remains scoped by the retained disjoint controls; distributed
conductor histories are not implemented by this increment.

Actual r08 local core SHA256:
`5e9b2a6d385ecb7d98c8ba306dd63e2fe139475c1f79685716dd2c2dbe7e2203`.
Actual r08 CUDA core SHA256:
`69b8fa61dead03bf134f3aaf815dffcd997084bbabdd8529d42faafade2386a5`.

Further acceptance requires spatial/time convergence, broad conductor CAD and
ports, scalable field/interface representation, distributed loss/PML/DFT state
and representative antenna comparisons. F02 remains incomplete.
