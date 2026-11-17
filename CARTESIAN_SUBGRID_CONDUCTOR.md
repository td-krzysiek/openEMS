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
The original r08 guard rejected a conventional repeated closing vertex.
The r09 follow-up ignores one such terminal sentinel only during validation,
preserving every supplied coordinate. Extra repetitions still reject. The
unchanged rejected r08 input passes; this fixes an experimental admission
restriction, not an ordinary upstream solver defect.
Direct C++ callers must update CSXCAD before operator setup, as native XML/
Python setup already does. Holes need supported shell/priority geometry or
simple decompositions, not self-touching polygon bridges.

`tests/test_cartesian_subgrid_shaped_surface.cpp` compares concave and rounded
contours, solid/hollow Z cylinders, polygon normals X/Y/Z, Cu18/Cu35, reversed
extrusion/winding, substrate loss and alternative units/origin. Every native
mask matches an independent analytic centre classifier; native coordinates
remain exact. CPU fields, loss and occupied volume match independent occupied-
box decompositions exactly. All hollow-bore, passive-energy and malformed-
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

The r09 closed-ring follow-up passes eight positive and six negative shape
controls, the composed 15-stage CPU/Metal replay, CTest40/40, Python14/14,
35 physical CUDA gates/CTest40 and 11,255 unchanged prior arrays. Maximum
Metal/CUDA shaped-field error remains 9.715e-7 within the unchanged 2e-4 gate.
Three older-CSXCAD scoped tests pass with three explicit missing-feature skips;
this is not a full old-dependency replay. R09 local core SHA256:
`e7f9501ef239d13f657e869fd9a04e6dade05cb9b40a5f47a03b936a3d995955`; CUDA core:
`5f273237e245d42dcde5340fc6be837b6b4cb1ad36fa4c43c1dd09358f4f60bd`. Curved RF is a separately tracked diagnostic,
not qualification supplied by this admission fix. The first completed U-line
ratio4 pair has band S21 phase differences of 0.095° for PEC and 12.7° for
smooth copper, with complex-S11 differences of 0.00024 and 0.149 respectively. Matched
native CAD/probes are checked. This unresolved experimental discrepancy is
retained. The boundary follow-up below explains its dominant contribution;
finer-grid and antenna acceptance remain open. It is not an established
ordinary upstream solver defect or proof that sidewall area alone determines
the RF error.

Further acceptance requires spatial/time convergence, broad conductor CAD and
ports, scalable field/interface representation, distributed loss/PML/DFT state
and representative antenna comparisons. F02 remains incomplete.

## Curved-line boundary diagnosis

The first coarse U-line RF pilot is not accepted as a curved-conductor
reference. Its matched continuous-CAD PEC pair agrees within 0.00024 in
complex S11 and 0.095° in S21 phase, but nominal smooth copper differs by
0.149 and 12.7°. A follow-up separates constitutive response from geometry
classification without changing any frozen runtime or RF input.

A synthetic PEC row-box control represents the exact opaque-cell mask.
Independent native preparation from the original XML checks all 126,390 leaves:
zero mismatches, including 10,008 opaque cells and 792 strip cells. Increasing
conductivity from 5.8e7 to 5.8e11 S/m approaches this matching PEC limit within
1.52e-5 in complex S11 and 0.00335° in S21 phase. Changing only the discrete PEC
boundary reproduces a 12.627° shift. With consistent PEC boundaries, the local
and ordinary incremental copper corrections differ by at most 0.02953° and
0.00258 dB in transmission.

For this fixture the dominant discrepancy is therefore geometric sampling:
native PEC classifies edges, while smooth-volume metal classifies occupied
cells and their exposed faces. Identical continuous CAD and probe planes do
not imply identical discrete conductor boundaries. This does not establish
an ordinary upstream defect, a general conductor-loss guarantee or converged
curved RF. Finer and selectively refined boundary/port controls remain open.

The row boxes are numerical diagnostics, never fabrication geometry or a
replacement for the continuous copper master. No production antenna runtime,
geometry or qualification is changed. The representative-array and scaling
acceptance requirements remain incomplete.

## Internal opaque refinement faces

A flat continuous opaque slab reproduced a separate surface-extraction
admission defect: coarse/uniform controls passed, but a local region rejected
an internal opaque/opaque face before checking its complete material coverage.
The r14 correction inspects every finest-lattice tile when the centre
neighbour does not cover the whole face. It omits that face only if every
intersecting neighbour is opaque. A centre-only shortcut would lose an exposed
corner; that negative case retains the original rejection.

`test_cartesian_subgrid_opaque_faces` covers all three normals, unequal 2/4
tangential ratios, coarse/uniform/local grids and two separate regions with
different ratios. Twelve positive cases check physical surface area, unique
faces and the 256-step invariant (maximum relative drift 2.670e-15). The new
regression fails on both the frozen R09 core and the actual R13 parent after
its coarse/uniform controls pass.
The partial-face negative remains rejected. No CAD, occupancy, PEC trace guard,
surface fit or update equation changes.

The composed local qualification passes 19 stages, CTest42/42 and Python14/14,
including 15 Metal constitutive cases and 11,255 identical prior arrays.
Physical CUDA replay passes all39 gates and CTest42/42, including58
constitutive comparisons in total:15 per device and14 per pair order. Separately, a contained native
Cu18 U-line with 961,584 leaves, 153,504 exposed faces and 614,016 traces now
prepares on Metal and completes eight actual CUDA steps on an RTX4090, without
host half-step fallback. Eight steps qualify allocation/dispatch only, not RF
accuracy or throughput. The CLI `--no-simulation` control intentionally returns
1 before engine creation and is retained as host preprocessing evidence only.

Local core: `d2de3973bfdbd81d616c7d211606054ec9fd09de807c9cadd522ef178624c791`;
CUDA core: `80f62648b4667c956e75aa0a34cab3d25c64edfa91d274190f15adad0cfed83a`.
General exposed partial-face subdivision and unresolved PEC interface traces
remain unsupported. Finer curved-line RF, scaling and representative antenna
qualification remain open; the new RF comparison uses a separately identified
graded ordinary reference with an unchanged outer/PML collar.


## Contained curved-line graded RF control

The r15 diagnostic retains continuous Cu18 U geometry, a contained16/16/1
region, common physical probes and5fs timestep. Its ordinary Cartesian
reference inserts the local fine axes while retaining the same coarse outer
collar/PML. This explicitly differs from R10's globally refined far field.
Four400,000-step/2ns Metal cases complete using frozen R14 and separately
pinned ordinary PR242 runtimes. No engine change follows from the comparison.

| Model | Maximum band complex S11 difference | S21 phase | S21 magnitude | Apparent Z difference |
|---|---:|---:|---:|---:|
| PEC | 0.002923 | 0.35744 deg | 0.00314 dB | 0.22067 ohm |
| Smooth copper | 0.012696 | 0.96489 deg | 0.04453 dB | 0.94782 ohm |

Values are maxima over24.000–24.250GHz. PEC passes the predeclared diagnostic
limits; copper misses complex-S11<0.01. Independent interior-voltage phase
error reaches0.94166deg, so the residual is not only current-port extraction.
Raw-TD prefix extraction agrees with the official full port extraction and
changes complex S11 by less than7.5e-9 from1ns to2ns. This bounds the observed
recorded tail, not arbitrary later evolution or spatial error.

The identical local-copper input also completes on a physical RTX4090: CUDA
versus Metal complex-S11 difference is9.62e-9 and phase difference2.67e-7deg.
That backend agreement is distinct from local-versus-ordinary RF accuracy.
Point-line launches are not calibrated characteristic-Z ports; conductor
boundary/fit and spatial convergence remain open. No antenna qualification,
production speedup or relaxed acceptance threshold is claimed.

## Compact curved-line two-resolution diagnosis

R23 retains a separate public terminal-scale U: physical Cu18,139.105um
width,250um centreline radius,0.5mm legs, fixed Z and400,000steps/5fs.
All eight frozen native inputs complete on one physical RTX4090. Local
execution uses qualified frozen R22 core SHA256
`778c6ee7478fc284889a0c26acc5a4f233cc06a0486b8d37c39c10f25db7fb14`;
the ordinary control uses its own pinned PR242 CLI/CSXCAD/plugin, retaining
LossyMetal. Input, runtime,152 raw-trace hashes and terminal exit0/no OOM are
independently checked. DFT reconstructs each trace at its own E/H timestamps;
the two half-weight current files are summed once.

| XY refinement | Conductor | Max band complex S11 difference | S21 phase | S21 magnitude | Original band gates |
|---|---|---:|---:|---:|---|
|8/8/1|PEC|0.001069911|0.006005deg|0.004534dB|Pass|
|8/8/1|Smooth copper|0.010262262|0.037780deg|0.051395dB|Fail S11|
|16/16/1|PEC|0.001264785|0.005358deg|0.005079dB|Pass|
|16/16/1|Smooth copper|0.006944176|0.008171deg|0.034716dB|Pass|

The original gates remain complex S11<0.01, phase<1deg and magnitude<0.1dB
over24.000–24.250GHz. Although the finer copper comparison passes, changing
XY resolution8→16 changes phase3.37–3.43deg and complex S11 by0.0308–0.0331
in both local and ordinary engines. All four two-level spatial-sensitivity
gates fail. The1.5ns-to-full2ns complex-S11 change is at most2.44e-8; every
late-prefix gate passes. The observed tail does not explain the spatial
change. This is unresolved curved-boundary/mesh/point-port sensitivity,
not an established backend defect or asymptotic convergence result.

These are capped diagnostic runs with fixed Z, not energy-converged antenna
qualification. The larger R15 failure is retained. No conductor/antenna
milestone is completed or threshold weakened by this additional evidence.
