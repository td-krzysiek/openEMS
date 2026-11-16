# Cartesian subgridding — implementation and qualification checklist

Last reviewed: **17 November 2026**.

This is the canonical detailed checklist for extending the Cartesian subgridding prototype
to the functionality of the regular Cartesian solver. Check an item only after
its stated evidence exists. Implementation and numerical qualification are
separate tasks; a passing test on one engine does not qualify another engine.

## Baseline and scope

- Development: [fork draft PR #29](https://github.com/td-krzysiek/openEMS/pull/29),
  branch `feature/cartesian-local-subgrid-20261112`. Keep subsequent increments
  in this PR rather than opening one PR per checklist item.
- Initial multi-region prototype commit: `907391792b55bee2a557c498eb6774812306e28d`.
  Subsequent qualified increments are recorded in the completion log below.
- Pinned master baseline: `92cd92a8d924c78a93453c25fa541384c2ac1d0e` (explicitly merged and requalified after the original `1d710c2` reference).
  Reconcile this inventory explicitly when adopting newer upstream commits.
- Public numerical summary and reproduction commands: [CARTESIAN_SUBGRID.md](CARTESIAN_SUBGRID.md).
- Executable controls: [C++ tests](tests/test_cartesian_subgrid.cpp),
  [graded-mesh controls](tests/test_cartesian_subgrid_mesh.cpp),
  [Python API tests](python/Tests/test_CartesianSubgrid.py),
  [native integration tests](python/Tests/Cartesian_Subgrid.py), and
  [Matlab/Octave test source](TESTSUITE/unittests/cartesian_subgrids.m).
  Raw qualification logs/build hashes remain in the contributor's local
  evidence archive; they are not attached to this PR. Checked items refer to
  the recorded runs summarized in the numerical documentation, not CI coverage
  of every platform or engine.
- Current capability includes a **serial CPU numerical reference**, with homogeneous
  lossless background, cell-aligned isotropic lossless material boxes and PEC/PMC
  outer faces, scoped parallel first-order Mur faces, and classical convolution
  PML on a checked homogeneous unrefined collar. Internal PEC boxes/sheets and axis-aligned cylinders/shells have
  scoped geometry, contact and numerical controls. One-edge positive resistive
  ports and local U/I integration have short two-port CPU controls. A public
  vacuum PEC dipole has scoped port, near-field, power and NF2FF CPU evidence;
  this does not qualify the project PCB antenna or general antenna workflows.
  A narrower native device field/source path has reference-backend and physical
  Metal controls. One-edge device RF port controls are separately qualified;
  the reproduced DC stress is corrected in scoped controls below. Scoped device PML and the public vacuum PEC dipole now have separate
  reference/physical Metal evidence; a separately pinned RTX3080 lane is recorded
  below. Losses, PCB antenna and scalability remain open.
- The new feature concerns Cartesian coordinates. Preserve the existing
  cylindrical modes and test them for regressions; implementing a new
  cylindrical subgrid is outside this scope.

### Agreed spatial interface

A region is specified by two opposite XYZ corners in mesh drawing units and
independent subdivision factors in X/Y/Z. Multiple regions may have different
factors. Corners must lie on base-grid lines, with at least one full base cell
between regions and between each region and the outer boundary. Reject invalid
coordinates rather than silently snapping them. No overlap, touching or nested
refinement is requested.

```python
fdtd.AddCartesianSubgrid(box=[1, 2, 2, 4, 6, 6], ratio=[2, 2, 4])
fdtd.AddCartesianSubgrid(box=[7, 2, 2, 10, 6, 6], ratio=[4, 1, 2])
```

`AddCartesianSubgrid` appends; `SetCartesianSubgrids([...])` replaces the list.
Prototype guards currently allow 64 regions, factors 1/2/4/8/16 per axis,
16,384 base intervals per axis, one million base cells, one million total leaves,
and 16,384 coupled degrees of freedom per interface factor. Above 1,024 edges,
half-bandwidth must be at most four and the one-RHS inverse bound must succeed. These are implementation guards, not
universal limits of subgridding. A shared timestep is retained.

## 0. Reference scheme and configuration — implemented

- [x] **R01** Local refinement in all three Cartesian dimensions, including
  independent per-axis factors and unchanged axes (`ratio=1`).
- [x] **R02** Several independent regions, each with its own subdivision;
  remove covered coarse cells without allocating global fine parent fields.
- [x] **R03** Constrained tangential electric traces and consistent Galerkin
  interface mass, shared curl/transpose, and a coupled timestep bound.
- [x] **R04** C++/Python list replacement and append API; optional XML region
  container; legacy single-box XML and default 2/2/2 compatibility.
- [x] **R05** Reject overlap, nesting, touching, invalid ratios, malformed XML,
  misaligned geometry and combined allocation-budget overflow.
- [x] **R06** Verify counts/local dimensions, all box faces, div(curl), region
  insertion-order independence, and a valid 64-region control.
- [x] **R07** Numerical controls: native Yee equivalence on uniform grids,
  cavity and packet convergence, million-step single-region and 100k-step
  multi-region invariant tests. See the evidence above for tolerances/results.
- [x] **R08** Native soft-E source, E/H probes and E/H field dumps on the CPU
  reference; existing cylindrical-cavity and regular CPU-engine controls pass.
- [x] **R09** Matlab/Octave configuration helper and test source are present.
- [ ] **R10** Execute the Matlab/Octave helper/test and check generated XML
  with the native solver. The local interpreter was unavailable.
- [ ] **R11** Extend interface-accuracy tests to oblique incidence, broader
  spectra, different ratios/aspect ratios and polarization. Existing coarse
  packet cases remain inaccurate; energy conservation alone is insufficient.

## 1. Geometry and material mapping

- [ ] **G01** Map CSXCAD geometry and material priorities onto actual local
  cells, edges and faces; preserve coordinate units and physical dimensions.
- [x] **G02** Implement spatially varying isotropic dielectric/permeability
  coefficients with a consistent interface mass and stability bound.
  Qualified for lossless, piecewise-constant leaf coefficients on serial CPU;
  native geometry currently accepts only cell-aligned volume boxes.
- [x] **G03** Implement internal PEC objects, including finite physical metal
  height, sheets where supported, vias, gaps and finite-width contacts.
  Qualified on serial CPU for aligned boxes/sheets and sampled axis-aligned
  cylinders/hollow barrels with aligned ends; unresolved interface traces fail.
  Later CPU/reference/Metal controls also qualify untransformed planar/linearly
  extruded PEC polygons, Cu18 port contacts and strict axial-alignment guards;
  see the [native polygon scope](CARTESIAN_SUBGRID.md#native-planar-and-extruded-pec-polygons).
  This is not general CAD, arbitrary interface/contact or antenna qualification.
- [ ] **G04** Qualify dielectric and conductor crossings of coarse/fine
  boundaries. Any interim placement restriction must fail explicitly.
- [ ] **G05** Support the anisotropic material coefficients available in the
  pinned regular Cartesian solver, or document a remaining parity exclusion.
- [ ] **G06** Add layered-medium, dielectric-interface, cavity and PEC-object
  controls against analytic references and the regular fine-grid solver.

Acceptance: independently inspect material occupancy and contacts, verify
geometry priorities, and demonstrate convergence and bounded evolution.
G01/G04/G06 are **partial**: native volume-box priorities, units and leaf
occupancy, dielectric crossings of refinement faces, analytic layered-cavity
convergence and a graded dielectric/native-Yee control pass. PEC adds volume,
sheet and barrel occupancy, resolved/broken Cu18 contacts, supported conductor
crossings, invariant/closed-cavity tests and four native-Yee shape controls.
General primitives, arbitrary contact/interface patterns and broader scattering
qualification remain open; unsupported interface patterns are rejected.
The polygon increment extends partial G01/G04/G06 through native membership,
axis/unit/origin controls and finite-height contacts, retaining the earlier
material-priority regressions. It preserves
the supplied contours and qualifies neither curved material interfaces nor
the spatial accuracy of a production PCB/array.

## 2. Nonuniform base mesh

- [x] **M01** Replace constant per-axis spacing assumptions with the actual
  base-grid line coordinates and local physical edge/face/volume measures.
- [x] **M02** Subdivide each selected base cell using its own dimensions,
  while retaining the agreed corner-alignment and separation rules.
- [x] **M03** Update point lookup, reconstruction, interface coefficients,
  source placement and timestep calculation for nonuniform base intervals.
- [x] **M04** Test graded meshes, small spacings, stretched cells and region
  boundaries against matched regular-grid references.

Acceptance: no hidden coordinate snapping, unexpected geometry change or
unaccounted minimum-spacing/timestep change.

## 3. Physical sources, ports and circuit elements

- [ ] **P01** Support the source types, directions, delays and waveforms used
  by the regular Cartesian solver; preserve the intended source normalization.
- [x] **P02** Implement lumped-port source/resistor coupling on local fields,
  plus voltage-path and current-contour integration across local cells.
  Scoped CPU qualification: native positive R on one independent local edge,
  zero-delay soft E source, exact local U paths and discrete dual H contours.
  Distributed boxes, multi-edge/hanging native resistors, overlapping native
  resistors and LC/modal ports remain explicitly unsupported.
- [ ] **P03** Support the remaining circuit/lumped RLC elements and qualify
  energy injection, stored energy and dissipation.
- [ ] **P04** Implement microstrip/coaxial/waveguide port workflows as
  applicable, including modal processing and reference-plane conventions.
- [ ] **P05** Validate single- and multiport power balance, complex impedance,
  reflection/transmission and passive terminations against regular fine grids.

Acceptance: qualify a short 50-ohm transmission line before an antenna.
Sampling reconstructed E/H on base observation coordinates is not, by itself,
a qualified physical port. P01/P03/P05 are **partial**: source normalization,
50-ohm resistors, complex Z/S11/S21, phase, polarity and two-port power balance
pass the scoped [port controls](CARTESIAN_SUBGRID.md#local-resistors-and-ui-integration).
The outside/inside-refinement controls preserve point-port cells in their
conventional Yee references. Energy/source-work/loss agrees within 1.46e-12
through 100k steps; broader sources, distributed R/LC/modal workflows and
material/interface placements remain open.

The frozen native-adapter Step/sample0 discrepancy is corrected by the
[source-tail increment](CARTESIAN_SUBGRID.md#native-source-waveform-tails).
Scoped CPU/reference/Metal and each physical RTX4090/RTX3080 waveform controls,
plus actual native Step scripting, pass. The full local 21-stage suite and
physical CUDA 26/26 CTest pass;17 older prefixes and11,255 arrays are retained.
Direct C++ periodic controls do not admit normal native periodic/steady-state
scripting. H-source/delay/TFSF and broader P01 coverage remain open.

## 4. Boundaries and radiation postprocessing

- [x] **B01** Implement and test PMC and supported symmetry boundaries.
  Scoped serial CPU reference: physical-face weak PMC, mixed closed-cavity
  convergence and all six mirrored-domain controls, including graded meshes.
  Exact-face normal E/tangential H are unilateral leaf reconstructions, not
  imposed point traces; staggered native-Yee PMC planes differ in position.
  See [boundary contract](CARTESIAN_SUBGRID.md#outer-pecpmc-boundaries).
- [ ] **B02** Integrate UPML and Mur boundaries. An initial requirement to keep
  refinement outside absorbing layers must be explicit and checked.
- [ ] **B03** Measure absorbing-boundary reflection over relevant frequencies,
  angles and polarizations, including nonuniform meshes.
- [x] **B04** Qualify DFT accumulation, field timing, spatial interpolation and
  integration weights on NF2FF surfaces sampling the local grid.
  Scoped serial CPU evidence: closed, node-interpolated E/H surfaces and all
  interpolation neighbours in homogeneous unrefined air, outside local
  refinement and PML; no mirror faces or material-crossing surfaces.
- [x] **B05** Compare small-antenna input response, near fields, radiation cuts,
  power and efficiency with regular fine-grid and analytic controls.
  Scoped to the public 1 GHz solid rectangular PEC dipole in vacuum, one
  local-edge 50-ohm port, finite mesh/time/PML sensitivity and thin-wire
  pattern sanity checks. This is not a PCB or lossy-efficiency qualification.

B02/B03 are **partial**: paired operator/engine Mur extensions support one
face or two opposite parallel faces, with a checked homogeneous, source-free,
unrefined two-cell collar. All six normals, a local packet, a 100k-step signed
boundary-work ledger and native Yee/normal/oblique-TE/graded controls pass.
The measured first-order reflection and temporal-window sensitivity are
[documented separately](CARTESIAN_SUBGRID.md#first-order-mur-absorbing-faces).
The separate classical convolution PML extension now supports any subset of
six faces, including edge/corner overlaps, with an unrefined background-matched
collar plus one inner guard cell. All six normals, a 100k-step local packet,
a 10k-step six-face corner and separate anisotropic-region controls pass on
serial CPU. Native 9.5/10/10.5 GHz normal, TE 30/45, TM 30, graded and local
controls pass; the largest sampled amplitude return is 1.71e-3. Native scope
rejections and the near-unit PEC reflection negative control pass.
See [PML contract, build identity and reproduction](CARTESIAN_SUBGRID.md#classical-convolution-pml).
This is not ordinary-Yee UPML or CFS parity. Adjacent Mur corners, mixed
Mur/PML, broader/grazing/near-cutoff spectra and general material traversals
remain open. The signed-work ledger is not a general passivity proof.

Acceptance: a radiating PEC antenna must work before treating the prototype as
an antenna solver. Correct E/H file creation alone does not qualify NF2FF.
The [PEC dipole controls](CARTESIAN_SUBGRID.md#scoped-pec-dipole-radiation)
now supply this evidence within the stated scope. B04/B05 do not cover
refinement-crossing or nonhomogeneous NF2FF surfaces, arbitrary CAD/ports,
the project array. The separately qualified device PML/dipole increment below
adds reference/physical Metal coverage within the same restricted surface scope. A failed coarser-Z candidate is
preserved and remains an expected-negative mesh control.

## 5. Scalable representation and CPU execution

- [ ] **C01** Use structured storage/update kernels for regular block interiors;
  reserve compact constraint metadata and solves for the interfaces.
- [ ] **C02** Replace or scale the current dense interface-factor approach;
  measure setup cost, factor memory and per-step cost before raising guards.
  Partial: exact banded Cholesky in the existing breadth-first order, local
  inverse-bound columns and one reusable block-sized solve workspace. The
  initial 1024-edge guard was retained at that milestone. A subsequent one-RHS signed
  inverse norm applies to positive bipartite graphs, with the old all-column
  fallback and scoped CPU/ordinary-Metal regression. General large-interface scaling
  remains open. The [bounded larger-interface increment](CARTESIAN_SUBGRID.md#bounded-larger-interface-factors)
  now admits at most 16,384 edges with strict bandwidth/bound-method guards,
  CPU/reference/Metal controls and unchanged small-model results. This does not
  qualify general large matrices or structured production fields. See also
  [factor representation and measurements](CARTESIAN_SUBGRID.md#banded-interface-factors).
- [ ] **C03** Reduce reference-only storage overhead, with regression coverage
  for any change to magnetic-face histories or numerical quadrature.
  Partial: banded interface factors/workspace and exact translated curl
  patterns now reduce setup/steady storage. Independent-edge addressing now has the separately qualified
  [structured slot index](CARTESIAN_SUBGRID.md#structured-independent-edge-address-index);
  other per-leaf topology, duplicated magnetic quadrature histories and field
  precision remain unchanged. The pattern
  representation is not a structured-interior kernel; CPU updates are slower
  in the scoped memory benchmark. See [exact stencil evidence](CARTESIAN_SUBGRID.md#exact-translated-curl-stencils).
- [ ] **C04** Integrate the regular CPU execution options, threading and SIMD
  where applicable; document remaining engine differences explicitly.
- [ ] **C05** Benchmark peak setup/steady-state memory, setup time and time to
  the same simulated physical duration and accuracy, against regular grids.
  Partial: three interleaved repetitions of seven standalone core cases against
  the previous PR29 source, with actual peak RSS, accounted storage, unchanged
  timestep and E/H/invariant checks. This is not a resource comparison against
  an optimized ordinary Yee/GPU engine or a production antenna.

Acceptance: do not infer production memory or speed gains from leaf counts.
Record minimum timestep, interface overhead, field precision and build identity.

Further partial C01/C03/U01–U03/U05: the
[spatial-curl increment](CARTESIAN_SUBGRID_BLOCKS.md) adds exact host block
admission, canonical E-slot maps and sparse exceptions. Default CPU retains its
old pattern loop. Actual shared GPU interiors use exact prepared coefficient
pairs shared by curl/transpose, with a byte gate retaining smaller old layouts.
Focused reference/Metal and physical RTX4090/RTX3080 geometry controls pass;
the composed local replay passes14 stages/CTest39/39/Python14/14 with11,255
prior arrays unchanged. Full CUDA replay passes33 gates and CTest39/39 on the actual RTX4090/RTX3080 pair. Interleaved Metal
layout timing is1.28–1.40x slower despite metadata savings; no speedup is claimed.
Expanded domain upload,
reference topology/field layout and production scaling remain open, so C01/C03
remain unchecked. The host-only r02 matrix passes13 stages/CTest38/38/Python14/14
with11,255 prior arrays unchanged; it does not qualify later device code.

Further partial C01/C03/C05/U01–U03/U05: exact transpose slots and stored
four-term order retain the r06 physical CUDA33-gate/CTest39 controls.
The composed r07 single-device row-stream preparation avoids both expanded
CSR operators, with local14-stage/CTest39/Python14 qualification and11,255
unchanged prior arrays. At491,520 leaves, three interleaved Metal repetitions
reduce peak RSS23.60% (1.615→1.234GB), with3.34% longer setup and effectively
unchanged timestep time. Full r07 CUDA replay passes33 gates/CTest39 on each physical RTX4090/RTX3080 and the retained lossless pair in both orders. See the
[exact-order/streamed evidence](CARTESIAN_SUBGRID_BLOCKS.md#exact-ordering-and-streamed-preparation).
CPU reference performance, general field/interface scaling, distributed
preparation and antenna qualification remain open; no item is completed by
this bounded memory result.

## 6. One GPU through the existing backends

- [ ] **U01** Define device-resident block fields, interface metadata and update
  scheduling within the existing GPU execution architecture.
- [ ] **U02** Implement regular interior updates and the constrained interface
  operations in the HIP backend, including the CUDA/NVIDIA path.
- [ ] **U03** Implement the equivalent operations in the Metal backend.
- [ ] **U04** Keep supported material/source/boundary histories and DFT state
  consistent on device; avoid per-half-step full-field transfers to the host.
- [ ] **U05** Validate CPU/GPU fields, ports, energy and radiation with stated
  numerical tolerances. Validate each physical backend independently.
- [ ] **U06** Measure actual VRAM, setup time, update/interface costs and total
  time to an equivalent result on one GPU. Check allocation cleanup.

Partial U01/U02/U03/U05/U06: shared portable point operations now run compact
lossless E/H, banded interface solves and ordered sources through the existing
backend. Independent FP32 full-field and 20k-step controls plus native E/H/FD
reconstruction pass on the reference backend and physical Metal (Apple M5 Pro).
A separately pinned [single RTX3080 CUDA lane](CARTESIAN_SUBGRID_CUDA.md#qualified-single-rtx3080-lane)
now passes all eleven original gates. It predates the later CAD/larger-interface
increments and does not qualify those source changes, other devices or multiple GPUs. The subsequent [device resistor scope](CARTESIAN_SUBGRID.md#scoped-device-resistor-ports)
adds native RF ports, midpoint damping and compensated source/loss ledgers.
The retained strong-DC failure of `0b255a8` is corrected by the separately
qualified [curl/timestep precision increment](CARTESIAN_SUBGRID.md#device-curl-and-timestep-precision-correction),
with the original field gate unchanged and permanent DC/cancellation regressions.
The later [exact device patterns](CARTESIAN_SUBGRID.md#exact-device-curl-patterns)
reduce device metadata with unchanged fields and scoped cost measurements.
Expanded host upload, host batch mirrors
and guards remain; this is not structured scalability or device DFT.
The subsequent [PML/dipole increment](CARTESIAN_SUBGRID.md#scoped-device-pml-and-dipole-radiation)
adds device absorber and small-antenna evidence only within its explicit scope. [Scope, source/build and tests](CARTESIAN_SUBGRID.md#scoped-native-device-execution).

Acceptance: one supported antenna case runs on one GPU without hidden CPU
execution. The ordinary master GPU path still working does not complete this
stage for subgridding.

## 7. Losses and remaining master functionality

- [ ] **F01** Implement and qualify finite electric/magnetic conductivity and
  dielectric loss, with passive energy dissipation and revised stability checks.
- [ ] **F02** Integrate the conducting-sheet/surface-conductor functionality
  present in the pinned master; check physical thickness and surface conventions.
- [ ] **F03** Support the pinned solver's dispersive material extensions,
  including their material histories and CPU/GPU implementations.
- [ ] **F04** Implement TFSF and remaining excitation/extension behavior.
- [ ] **F05** Complete the remaining probe/dump quantities, modal processing,
  SAR and processing features from an explicit pinned-master inventory.
- [ ] **F06** Run the complete applicable Python and Matlab/Octave suites;
  preserve existing Cartesian/cylindrical and non-subgridded backend behavior.
- [ ] **F07** Complete a feature/engine compatibility matrix. Every unsupported
  combination must remain explicit; do not call partial coverage full parity.

Acceptance: each extension needs numerical controls on every claimed backend.
Additional downstream conductor, runtime-control and multi-GPU extensions need
their own integration checks; master parity does not automatically qualify them.

Partial F01: the [paired constant-conductivity implementation](CARTESIAN_SUBGRID_CONDUCTIVITY.md)
passes independent CPU/reference/Metal and physical RTX4090/RTX3080 recurrence
controls, including full XYZ refinement, pure electric/magnetic losses, ports
and real sources. Native SI box mapping passes analytic decay/time convergence,
priority/unit and scripting controls. The combined local 25-stage qualification,
CTest 29/29, Python 14/14 and physical CUDA 13-gate/CTest 29/29 suite pass;
21 older prefixes and 11,255 prior arrays are unchanged. A separate native closed-mode
pilot shows second-order frequency convergence, with residual spatial waveform
error explicitly recorded. The later prepared-support increment qualifies separated
conductivity/PML on CPU/Metal and each RTX4090/RTX3080: all26 local stages,
all15 CUDA gates and CTest30/30 pass;22 older prefixes and11,255 arrays remain
unchanged. Overlapping/Mur loss, broader propagation, surface/ADE conductors
and distributed histories remain open, so F01 is not checked off.

Partial F02: the [passive physical-path surface-state component](CARTESIAN_SUBGRID_SURFACE.md)
passes a full monolithic oracle, seven CPU/reference/Metal and physical
RTX4090/RTX3080 cases, two disjoint surface/PML controls and fourteen negative
guards. All30 local stages, all19 CUDA gates, CTest33/33 and Python14/14 pass;
11,255 prior arrays are bit-identical. Native physical conductor faces/interiors,
short RF lines, roughness and distributed surface histories are not admitted.
F02 stays unchecked. A reproduced reconstruction-roundoff defect in PR29 is
corrected without moving mesh/CAD coordinates or weakening mass-graph guards.

Further partial F02/G01/P02/P05: the
[native smooth-conductor/path-port increment](CARTESIAN_SUBGRID_CONDUCTOR.md)
implements physical opaque volume boxes, same-metal unions, priority-clipped
normal thickness and native complete-edge R/N paths. Candidate r02 passes all
34 local stages, CTest37/37/Python14/14, all27 physical CUDA gates and CTest37/37;
11,255 older arrays remain bit-identical. This is single-device conductor
coverage plus retained lossless two-device controls, not distributed histories.
The current r03 adds old-CSXCAD compatibility: its real old-dependency build
passes 35 tests with exactly two optional-feature skips. All 34 modern local
stages, CTest 37/37/Python 14/14, all 27 physical CUDA gates and CTest 37/37 pass;
11,255 prior arrays remain bit-identical. Short RF controls agree within each mesh, but the two-level
point-port fixture is not yet spatially converged. Broad CAD/ports, scaling and
antenna acceptance remain open; F02 is unchecked.

Further partial F02/G01/G04/G06/U04/U05: the
[native shaped increment](CARTESIAN_SUBGRID_CONDUCTOR.md#native-shaped-conductors)
adds simple extruded smooth-volume contours and axial cylinders/shells with
independent occupancy, unchanged coordinates, exact occupied-box decomposition,
hollow-bore, fields/passive energy and malformed-shape controls. Composed r08:
15 local stages, CTest40/40/Python14/14,35 physical CUDA gates/CTest40 and11,255
unchanged arrays. A separate actual older-CSXCAD build passes3 scoped native
controls with3 explicit missing-feature skips. This is centre-classified CAD,
not conformal or curved-line RF qualification. Sidewall-area bias and required
RF convergence remain explicit; F02/G01 and representative-antenna items remain
incomplete. No production adoption or allocation-guard increase follows.

The r09 follow-up accepts one native closing sentinel without editing CAD and
retains malformed-contour rejection. Eight shape/six negative controls and the
composed CPU/Metal 15-stage/CTest40/Python14 plus physical CUDA 35-gate/CTest40
replay pass; 11,255 prior arrays remain identical. This completes the scoped
closed-ring admission fix only, not curved RF, distributed histories or F02.

Further F02/A02/A03 diagnosis: the coarse U-line copper discrepancy is largely
reproduced by the different edge/opaque-cell boundary classifications. An
actual native mask check and high-conductivity PEC-limit control separate this
from the much smaller constitutive correction. See [boundary diagnosis](CARTESIAN_SUBGRID_CONDUCTOR.md#curved-line-boundary-diagnosis).
This is a retained accuracy gap requiring finer boundaries/ports, not curved
RF, representative-antenna or production qualification.

Further partial G01/G04/F02: [internal opaque refinement faces](CARTESIAN_SUBGRID_CONDUCTOR.md#internal-opaque-refinement-faces)
are omitted only after complete tile coverage is established. The reproduced
slab admission failure is corrected; partially exposed faces and unresolved PEC
traces still reject. Twelve positive face controls and a negative corner case,
local19-stage/CTest42/Python14 qualification and11,255 exact prior arrays pass.
The contained native Cu18 U prepares on Metal and completes eight CUDA steps;
this is allocation/dispatch evidence, not curved RF acceptance. Physical CUDA
compatibility replay passes39 gates/CTest42, including58 constitutive cases. No general geometry, conductor or antenna
item is completed by this admission fix.

## 8. Representative antenna qualification

- [ ] **A01** Freeze a regular-grid PEC antenna reference, geometry/material/
  port hashes, comparison frequencies and numerical acceptance tolerances.
- [ ] **A02** Select refinement regions and verify actual metal, dielectric,
  via, contact and source occupancy against the frozen reference.
- [ ] **A03** Compare complex Gamma/Z, full-band S11, power balance and near
  fields before interpreting changes in radiation.
- [ ] **A04** Compare beam directions, widths, both first and worst sidelobes,
  polarization quantities and gain/efficiency with matched definitions.
- [ ] **A05** Establish spatial and temporal convergence; angular sampling is
  not a substitute for spatial/temporal qualification.
- [ ] **A06** Repeat on the selected lossy-conductor model and demonstrate
  acceptable errors and actual end-to-end resource savings.
- [ ] **A07** Record qualified runtime/input hashes before adopting the new
  scheme for subsequent antenna decisions. Preserve frozen prior results.

## 9. Multiple GPUs — follow-up after single-device qualification

- [ ] **D01** Define ownership of blocks, constrained interface states and
  material histories using a backend-independent domain contract. This is follow-up work,
  not functionality already present in the pinned master.
- [ ] **D02** Implement exchanges and synchronization through a separately qualified
  transport layer, including host staging when P2P is unavailable.
- [ ] **D03** Qualify heterogeneous device sizes/speeds, uneven partitions,
  memory limits, and fixed ownership throughout a run.
- [ ] **D04** Preserve full-step observation, pause/resume, DFT snapshots and
  NF2FF consistency across domains, with separately qualified runtime controls where applicable.
- [ ] **D05** Compare fields/RF results with one device and benchmark exchange
  overhead, peak host/device memory and time to an equivalent result.

Partial D01/D02/D03/D04/D05: the [C++ domain lane](CARTESIAN_SUBGRID_DOMAINS.md)
implements fixed whole-component ownership, owned/ghost state, backend-neutral
coordinator/worker contracts and preallocated host-staged boundary exchange.
One through four logical domains, XYZ cuts and mixed CPU-reference/GPU workers
have complete-field/epoch/failure controls on reference, physical Metal and
one physical RTX3080. Local 19-stage qualification passes, with 13 older runtime
prefixes and 11,255 arrays unchanged; CUDA also passes CTest 23/23 and single-device
fields/ports/DC. These are logical-domain controls on one physical GPU per lane,
not physical multi-GPU or distributed RF/processing qualification. Subsequent
scoped coupled-resistor/source-ledger controls pass on three logical workers
on one Metal GPU and one RTX3080, retaining complete small port solves on one
owner. All 20 local stages and CUDA CTest 24/24 pass; 14 older prefixes and 11,255
arrays remain unchanged. Domain absorber/bulk-loss histories, DFT/NF2FF,
collective pause/resume, native selectors,
asynchronous transport and scalability remain open. The D items are incomplete.

Further partial D02/D03/D05: [physical mixed-device context controls](CARTESIAN_SUBGRID_DOMAINS.md#physical-mixed-device-contexts)
reproduce and correct missing device selection in the experimental local
adapter. A physical RTX4090/RTX3080 pair passes both device orders, uneven
ownership, field/port/source-ledger and failure guards with identical
single-device fields; CUDA CTest 24/24 passes. This is bounded C++ field/port
qualification, not native distributed PML/loss/DFT or antenna support. The
host-only activation microbenchmark is not an end-to-end speedup result.
All 20 local compatibility stages pass, with CTest 24/24, Python 14/14,
15 older prefixes and 11,255 numerical arrays retained unchanged. The D items
remain partial because histories, native distributed processing and scaling
are still open.

Further partial D01/D02/F01/F02: [prepared conductivity domains](CARTESIAN_SUBGRID_DOMAINS.md#prepared-conductivity-domains)
retain owned bulk/surface histories and complete-epoch loss/stored-energy
observations. Fourteen reference and fifteen Metal constitutive controls pass;
the composed local19-stage/CTest41/Python14 replay retains11,255 exact arrays.
Physical CUDA39-gate/CTest41 qualification passes, including58 constitutive
comparisons across each GPU and both physical pair orders. FP32 fields equal
one-device results exactly; maximum double field error5.830e-6 and balance
2.638e-7 retain the2e-4 gates. Native distributed
selectors, absorber/DFT/pause histories, RF accuracy and production scaling
remain open; no D or antenna milestone is completed by this increment.

Further partial D01/D02/D03/D05/B02/U05: cold PML-domain composition passes
20 local stages, CTest43/Python14 and11,255 exact arrays. R16's physical pair
exposes a compact-PML device-context defect after35 retained gates pass. The
[r17 owning-device correction](CARTESIAN_SUBGRID_DOMAINS.md#compact-pml-owning-device-correction)
passes the same local replay and five physical4090/3080 PML controls in each
order. Full physical replay passes43 gates/CTest43, including24 PML and58
constitutive comparisons, with checked terminal source/runtime identities. Native
selectors/DFT/pause, general scaling and antenna acceptance remain open.

The [contained curved-line RF control](CARTESIAN_SUBGRID_CONDUCTOR.md#contained-curved-line-graded-rf-control)
passes the PEC diagnostic targets. Copper has0.012696 complex-S11 difference,
0.96489deg phase and0.04453dB magnitude error against the matched graded
ordinary reference; the0.01 S11 gate is missed. The recorded1-to-2ns tail is
negligible, but spatial/conductor accuracy is not yet accepted. Physical
CUDA/Metal agreement on the same local copper input is separately verified.
No F02 or representative-antenna item is completed by these controls.

Further partial D01/D04/D05/P05/B04/U05: the [native local-domain adapter](CARTESIAN_SUBGRID_RUNTIME.md)
connects normal CLI/Python runs to cold loss/PML domains and complete-batch
host TD/FD/physical U/I processing. R18-r01 passes20 local stages, CTest44,
Python14 and11,255 exact arrays. R02's explicit one-worker `none` sentinel
passes its own CTest44 and native reference/Metal/CLI/Python/mixed controls.
A scoped public dipole has native logical-Metal-domain NF2FF evidence. R18-r01
also passes48 physical CUDA gates/CTest44, including both pair orders and16
actual native Python cases, with terminal source/runtime rechecks. R02's
parser-only physical follow-up passes seven gates/CTest44, native Python and
the actual one-worker CLI control. No representative array,
device DFT, external collective pause or production scaling is qualified;
the D/A/U items remain incomplete.

Further partial C01/C03/C05/D01/D05: [streamed global domain preparation](CARTESIAN_SUBGRID_DOMAINS.md#streamed-global-domain-preparation)
removes both temporary global expanded curl operators, preserving local
uploads, metadata/order and update kernels. Six exact metadata/256-step field
cases per reference/Metal backend, native8/16-case TD/FD/port controls and
CTest45/45 pass. Against the actual frozen R18-r02 runtime, three interleaved
491,520-leaf setup measurements reduce peak RSS23.68%, with10.07% longer whole
setup. The original R19 CUDA fixture fails before backend creation; the corrected
test separately passes reference/Metal on the unchanged core. Its local
source/binary scope correction is documented with the domain evidence.
Physical qualification remains open. Local expanded uploads, hash
maps, host fields and general interfaces remain scaling work; no C/D/U/A item
or allocation guard is completed/raised by this scoped result.


Further partial C01/C03/C05/D01/D05: [streamed local domain preparation](CARTESIAN_SUBGRID_DOMAINS.md#streamed-local-domain-preparation)
retains exact local row patterns directly and borrows them only for backend
upload, without changing kernels or ownership. The complete R20-r03 source
was fixed/hashed before compilation. Six reference/six Metal exact pattern,
transfer and256-step port-field cases,13 pre-allocation malformed-package
guards, native8/16-case TD/FD/ports and fresh CTest45/45 pass. An actual frozen
R19 comparison at491,520 leaves reduces peak RSS30.37%, with7.93% longer whole
setup; a graded two-material control saves only2.88% RAM. A960,000-leaf setup
passes below1.5GiB, without FDTD/RF acceptance at that size. Physical CUDA
passes14 gates/CTest45,24 exact metadata/field cases, both pair orders,
13 malformed-package guards and16 actual native Python cases, with terminal
source/runtime rechecks. Global host fields/topology and large-interface
scaling remain open; guards and production runtimes are unchanged.

Further partial C03/C05/D01/D05: [shared owned-state domain addresses](CARTESIAN_SUBGRID_DOMAINS.md#shared-owned-state-domain-addresses)
replace owned hash entries with shared direct slots; private ghost hashes and
the independent expanded oracle remain. Sanitized map/negative controls,
12 reference/Metal exact metadata/field cases, native8/16 TD/FD/U/I and fresh
CTest46/46 pass. Interleaved comparison against frozen R20-r03 reduces
491,520-leaf peak RSS24.72% and whole setup9.78%; the graded61,440-leaf
reduction is9.78% RAM/7.47% setup. A960k-leaf setup-only control remains below
1.5GiB. Physical CUDA passes14 gates/CTest46,24 exact metadata/field cases,
13 package guards, both physical pair orders and16 native Python cases,
with terminal source/runtime rechecks. This is preparation-only memory
work; reference fields/topology, interfaces and actual antenna
convergence/scaling remain open. No guard increase or production adoption.

Further partial C03/C05/U04/D05: [lazy serial reference workspace](CARTESIAN_SUBGRID_DOMAINS.md#lazy-serial-reference-workspace)
keeps CPU sweep/constitutive scratch empty during device-only preparation;
host mirror fields and physical histories remain. First-use/capacity/negative
controls,12 reference/Metal exact cases, native8/16 and CTest47/47 pass.
Interleaved comparison against actual frozen R21 reduces491,520-leaf peak
RSS7.04%, with exactly35,047,680B lower accounted cold grid storage. The graded
control reduces RSS8.17%; no update speedup is claimed. GPU-free CUDA build
passes; physical R22 passes14 gates/CTest47,24 exact cases,13 package guards,
both4090/3080 orders and16 native Python cases, with terminal source/runtime
rechecks. No allocation-guard increase, antenna adoption or large-scale
acceptance follows.

## Milestone decisions

| Milestone | Completion condition | Status |
|---|---|---|---|
| Reference prototype | R01–R08 demonstrated on small CPU controls | Achieved, limited scope |
| Useful PEC antenna on CPU | Required geometry/mesh/port/PML/NF2FF items, with small-antenna evidence | Achieved for the public vacuum PEC dipole only; broader antenna workflows open |
| Useful PEC antenna on one GPU | Previous milestone plus scalable representation and a qualified device backend | Open |
| Representative lossy antenna | Required conductor support and A01–A07 | Open |
| Pinned-master Cartesian parity | Full applicable feature/engine matrix and regression suites | Open |
| Production multi-GPU subgridding | Relevant single-device physics qualified, then D01–D05 | Open |

Next implementation milestones: further interface scaling and required PCB
CAD/material mapping, while structured interior field kernels remain open. Exact banded interface factors
and the one-solve inverse bound have scoped CPU evidence. Exact translated
curl patterns reduce memory but increase measured CPU update time; current
allocation guards remain. The public PEC dipole
provides a port/radiation regression for those changes. Dielectric PCB geometry,
broader scattering/ports/NF2FF placements and representative antenna spatial
qualification remain open, followed by lossy antenna validation and parity work.
Nonuniform base geometry is implemented and qualified on small homogeneous
controls; application compatibility still requires the material, port and
boundary gates above. These are dependency gates, not a
commitment that every master feature must precede the first useful PEC model.

Further partial D01/D02/B02/U04: the [cold PML domain candidate](CARTESIAN_SUBGRID_DOMAINS.md#prepared-pml-domains)
implements target-owned convolution histories with explicit disjoint loss
composition, ordered field-only halos and complete-epoch signed boundary work.
Its original R16 local controls passed, while the physical pair exposed the
context defect corrected and qualified in R17 above. Later native selection
and processing evidence belongs to R18, not this initial component milestone.
Collective runtime control and representative-array/scaling acceptance remain
open; no full D/B/F/U/A item is completed here.

## Completion log

Append one row per qualified increment. Link actual evidence and record the
source/build identity and applicable engines. Reopen affected items when a
numerical change invalidates their prior qualification.

| Date | Items | Source | Evidence / limits |
|---|---|---|---|
| 2026-11-12 | R01 (XYZ portion), R03/R07/R08, initial reference | `7bade0c` | [Numerical summary](CARTESIAN_SUBGRID.md#tests-and-observed-accuracy); homogeneous CPU/PEC-wall tests only |
| 2026-11-12 | R01/R04/R05, independent axis ratios | `bd214a3` | [Axis-ratio controls](tests/test_cartesian_subgrid.cpp) and [API regression](python/Tests/test_CartesianSubgrid.py) |
| 2026-11-12 | R02/R04–R09, multiple separate regions | `9073917` | [Multi-region controls](tests/test_cartesian_subgrid.cpp) and [native integration](python/Tests/Cartesian_Subgrid.py); R09 is source availability, R10 remains unexecuted |
| 2026-11-12 | M01–M04, nonuniform base lines | This PR, graded-mesh increment | [Physical geometry and numerical tests](tests/test_cartesian_subgrid_mesh.cpp), [native source/translation controls](python/Tests/Cartesian_Subgrid.py), [results](CARTESIAN_SUBGRID.md#tests-and-observed-accuracy); serial homogeneous CPU only |
| 2026-11-12 | G02; partial G01/G04/G06 | This PR, lossless-material increment | [Material/occupancy/invariant/cavity controls](tests/test_cartesian_subgrid_material.cpp), [graded native Yee](tests/test_cartesian_subgrid_mesh.cpp), [native priorities/units/guards](python/Tests/Cartesian_Subgrid_Materials.py); serial CPU, cell-aligned isotropic lossless boxes only; internal PEC and broad scattering still open |
| 2026-11-12 | G03; extended partial G01/G04/G06 | This PR, internal-PEC increment | [Cu18 contact/PEC invariant/cavity/native-Yee controls](tests/test_cartesian_subgrid_pec.cpp), [native integration and guards](python/Tests/Cartesian_Subgrid_PEC.py); serial CPU, aligned boxes/sheets and sampled axial cylinders/shells; not physical-port, broad scattering or GPU qualification |
| 2026-11-12 | P02 scoped; partial P01/P03/P05 | This PR, physical-port increment (parent `4186af8`) | [Core power/integral controls](tests/test_cartesian_subgrid_ports.cpp), [native two-port comparisons](python/Tests/Cartesian_Subgrid_Ports.py), [contract and numerical results](CARTESIAN_SUBGRID.md#local-resistors-and-ui-integration); one independent native R edge, serial CPU; broader ports remain open |
| 2026-11-13 | B01 scoped PMC/symmetry | This PR, PMC boundary increment (parent `19dd0f1`) | [Core cavity/energy/face controls](tests/test_cartesian_subgrid_boundaries.cpp), [native mirror-domain controls](python/Tests/Cartesian_Subgrid_Boundaries.py), [boundary planes and reconstruction limits](CARTESIAN_SUBGRID.md#outer-pecpmc-boundaries); serial CPU only. PML/NF2FF/GPU remain open. |
| 2026-11-13 | Partial B02/B03, parallel Mur faces | This PR, paired Mur extension increment (parent `9de9020`) | [Pulse/energy/half-step controls](tests/test_cartesian_subgrid_mur.cpp), [native reflection and scope guards](python/Tests/Cartesian_Subgrid_Mur.py), [contract and limitations](CARTESIAN_SUBGRID.md#first-order-mur-absorbing-faces); serial basic CPU only. Adjacent corners, PML/NF2FF/GPU remain open. |
| 2026-11-15 | Extended partial B02/B03, classical convolution PML | This PR, paired PML extension increment (parent `eaa1a31`) | [Six-normal/100k-packet/10k-corner/multiple-region controls](tests/test_cartesian_subgrid_pml.cpp), [native normal/TE/TM/graded/local reflection and rejection controls](python/Tests/Cartesian_Subgrid_PML.py), [contract and build evidence](CARTESIAN_SUBGRID.md#classical-convolution-pml). Scoped serial basic CPU only; not UPML/CFS parity, NF2FF, antenna, losses or GPU qualification. |
| 2026-11-15 | B04/B05 scoped, public PEC CPU antenna milestone; extended partial P05 | This PR, radiation qualification increment (numerical source `55fe3e2`) | [Ten-case native dipole regression](python/Tests/Cartesian_Subgrid_Radiation.py), [results, actual build hashes and limitations](CARTESIAN_SUBGRID.md#scoped-pec-dipole-radiation). Existing numerical engine/NF2FF algorithms unchanged. Homogeneous unrefined closed surfaces, vacuum solid PEC dipole and one-edge port only; coarser-Z case retained as an expected negative. Not project-array, losses, GPU, speedup or production qualification. |
| 2026-11-15 | Partial C02/C03/C05, banded interface representation | This PR, factor increment (parent `fb98953`) | [Frozen dense-oracle/guard regression](tests/test_cartesian_subgrid_factor.cpp), [interleaved core benchmark](benchmark/subgrid_compare.py), [factor contract and qualified build](CARTESIAN_SUBGRID.md#banded-interface-factors). Same SPD mass/curl/timestep scheme; 1024-edge guard unchanged. Serial CPU only, no structured-interior, GPU or production performance qualification. |
| 2026-11-15 | Further partial C02/C05, signed inverse-norm setup | This PR, signed-bound increment (parent `571766f`) | [Independent dense-inverse controls](tests/test_cartesian_subgrid_bound.cpp), frozen factor oracle, [derivation/build and regression results](CARTESIAN_SUBGRID.md#signed-inverse-norm-setup). One solve for positive bipartite graphs; old general fallback and guards retained. CTest13/13, Python13/13, long9/9, native8/8 and ordinary4/4; ten frozen dipole arrays unchanged. Ordinary Metal is a regression only, not subgrid GPU qualification. |
| 2026-11-15 | Further partial C01/C03/C05, exact translated curl patterns | This PR, stencil increment (parent `cb36dde`) | [Independent CSR oracle](tests/test_cartesian_subgrid_curl.cpp), [memory/performance and loaded build evidence](CARTESIAN_SUBGRID.md#exact-translated-curl-stencils). CTest 14/14, Python13/13, long9/9, native8/8, ordinary4/4; ten frozen dipoles unchanged. Memory reduced about40–50%, CPU updates about9–12% slower. Structured interiors, topology scaling and GPU remain open. |
| 2026-11-15 | Reconcile pinned master; preserve scoped CPU/ordinary GPU coverage | Merge upstream `92cd92a` into parent `f5d1b16` | [Integration and rounding evidence](CARTESIAN_SUBGRID.md#upstream-92cd92a-integration). CTest 14/14, Python13/13, long9/9, native8/8, ordinary4/4. Ten frozen dipoles within 6.45e-7, not bit-identical after upstream disables extension FMA. Isolated old-contraction identity control recovers exact frozen arrays. This does not qualify subgrid GPU/lossy/PCB operation. |
| 2026-11-15 | Partial U01/U02/U03/U05/U06, shared native device field/source path | This PR, native-device increment (parent `039b128`) | [Shared operations](FDTD/kernels/subgrid_update.h), [independent FP32 oracle and fields](tests/test_cartesian_subgrid_device.cpp), [native probes/dumps and guards](python/Tests/Cartesian_Subgrid_Device.py), [scope/build evidence](CARTESIAN_SUBGRID.md#scoped-native-device-execution). CTest 15/15; physical fields 4/4; Python 13/13; native device 4/4, CPU base plus six groups, two frozen CPU dipoles unchanged, ordinary three groups. Metal on Apple M5 Pro and explicit reference backend; no physical HIP/CUDA, ports/absorbers/losses, device antenna or production scaling qualification. |
| 2026-11-15 | Extended partial P01/P03/P05 and U04/U05/U06: scoped device RF resistors | This PR, device-port increment (parent `4a7cf69`) | [Core passive/work controls](tests/test_cartesian_subgrid_device_ports.cpp), [native RF ports](python/Tests/Cartesian_Subgrid_Device_Ports.py), [factor range guards](tests/test_cartesian_subgrid_device_precision.cpp), [scope/build and open DC failure](CARTESIAN_SUBGRID.md#scoped-device-resistor-ports). CTest17/17, Python13/13, physical8/8, native fields/ports4+4, CPU base+six groups, two frozen dipoles unchanged, ordinary3/3. Strong-DC field target remains failed; not full source accuracy, absorber/lossy/antenna or scaling qualification. |
| 2026-11-15 | Extended partial U05/U06: device curl/timestep precision | This PR, precision increment (parent `0b255a8`) | [Correction, unchanged numerical gates, exact runtime and measured costs](CARTESIAN_SUBGRID.md#device-curl-and-timestep-precision-correction). CTest18/18, Python13/13, physical fields/ports4+4 plus DC on both backends, native4+4, CPU base+six groups, two unchanged frozen dipoles and ordinary3/3. Prior DC failures retained; no arbitrary-source, device antenna/absorber/loss or scalability qualification. |
| 2026-11-15 | Extended partial B02/B03/B04/B05 and U04/U05/U06: device PML and public dipole | This PR, paired device PML increment (parent `ca6f755`) | [Core PML controls](tests/test_cartesian_subgrid_device_pml.cpp), [native PML](python/Tests/Cartesian_Subgrid_Device_PML.py), [native dipole](python/Tests/Cartesian_Subgrid_Device_Radiation.py), [exact build/scope](CARTESIAN_SUBGRID.md#scoped-device-pml-and-dipole-radiation). CTest19/19, Python14/14, long reference/Metal packet/corners, native18 PML/four dipole comparisons, two frozen CPU dipoles unchanged, ordinary3/3. CUDA harness/guards are prepared, not compiled or physically qualified on NVIDIA. No representative PCB, losses/conductor, structured or multi-device qualification. |
| 2026-11-15 | Further partial C03/C05 and U05/U06: exact device curl patterns | This PR, device-pattern increment (parent `bb2edbb`) | [Bit-preserving metadata/gather controls](tests/test_cartesian_subgrid_device_patterns.cpp), [qualification and measured costs](CARTESIAN_SUBGRID.md#exact-device-curl-patterns). CTest20/20, Python14/14, full prior CPU/reference/Metal matrix; 11,255 prior numerical arrays unchanged. Small local24 batch32.6% shorter and accounted device memory77.2% lower; process RSS3.6% lower. Expanded host upload, structured kernels/guards, PCB/lossy and physical CUDA/multi-device remain open. |
| 2026-11-15 | Further partial C01/C03/C05: structured edge addressing | This PR, address-index increment (parent `f5aedf2`) | [Address/oracle regression](tests/test_cartesian_subgrid_edge_index.cpp), [build/evidence/benchmarks](CARTESIAN_SUBGRID.md#structured-independent-edge-address-index). CTest 21/21, Python14/14, standalone ASan/UBSan and full prior CPU/reference/Metal matrix; 11,255 numerical arrays unchanged. Seven interleaved core cases show setup 43–51% shorter; update unchanged. Structured field kernels, larger guards, PCB/lossy, physical CUDA and multi-device remain open. |
| 2026-11-15 | Extended partial G01/G03/G04/G06: native planar/extruded PEC CAD | This PR, polygon increment (parent `a846ebf`) | [Native fields/contacts/alignment/shape controls](python/Tests/Cartesian_Subgrid_Polygons.py), [build, results and limits](CARTESIAN_SUBGRID.md#native-planar-and-extruded-pec-polygons). All17 qualification stages pass; CTest21/21, Python14/14, full previous CPU/reference/Metal and ordinary matrix. Eleven prior core/NF2FF prefixes and11,255 numerical arrays retained. No polygonal material, transform, lossy, representative PCB or physical CUDA/multi-device qualification for this increment. |
| 2026-11-15 | Further partial C02/U02/U05/U06: bounded larger interfaces | This PR, larger-interface increment (parent `3a502cc`) | [Resource-policy and CPU controls](tests/test_cartesian_subgrid_factor.cpp), [three-orientation/device source controls](tests/test_cartesian_subgrid_device_large.cpp), [runtime and exact limits](CARTESIAN_SUBGRID.md#bounded-larger-interface-factors). All 18 stages pass, CTest 22/22 and Python 14/14; twelve prior prefixes and 11,255 arrays unchanged. Actual 16,384-edge CPU/Metal controls pass, 16,386-edge construction rejects. No new arithmetic, CUDA, lossy, structured-field or representative antenna qualification. |
| 2026-11-15 | Extended partial U02/U04/U05/U06: single physical CUDA | Numerical source `02b94ef` + test-only `a846ebf` include | [Pinned source/build, errors and limits](CARTESIAN_SUBGRID_CUDA.md#qualified-single-rtx3080-lane). All eleven original gates pass across four immutable bounded attempts: CTest 21/21, Python 14/14, fields/ports/DC, native ports/PML, public dipole/NF2FF and ordinary controls. Earlier timeouts preserved; final identities pass. RTX3080 SM86 only, built-in HIP/nvcc, no hidden fallback. Later CAD/larger-interface CUDA, RTX4090 hardware, plugin mode, losses/PCB and multi-device remain open. |
| 2026-11-16 | Partial D01–D05: compact-domain ownership and host transport | This PR, domain increment (configure-time parent `1a8f48f`) | [Contract, tests, runtime and limitations](CARTESIAN_SUBGRID_DOMAINS.md). Local 19-stage suite, CTest 23/23 / Python 14/14; 13 older prefixes and 11,255 arrays unchanged. Reference, one physical Metal GPU and one physical RTX3080 support 1–4 logical domains and mixed CPU-reference/GPU workers with identical one-device fields. ASan/UBSan domain-source and actual 16,384-edge Metal controls pass. Physical multi-GPU, domain ports/PML/losses/DFT, native selectors, scalability and representative antenna qualification remain open. |
| 2026-11-16 | Extended partial D01–D05 and P05: coupled domain ports/work | This PR, port-domain increment (configure-time parent `15daf56`) | [Source, gates, hashes and limits](CARTESIAN_SUBGRID_DOMAINS.md#coupled-resistors-and-source-ledgers). All 20 local stages, CTest 24/24/Python 14/14, CUDA six gates/CTest 24/24; 14 older prefixes and 11,255 arrays unchanged. Three logical workers on reference/one Metal GPU/one RTX3080, including mixed CPU-reference/GPU and source-only workers, warm ledgers and strong DC. Identical one-device fields; maximum double error 1.941e-5. Physical multi-GPU, native domain selectors, domain PML/bulk loss/DFT and representative antenna/scaling qualification remain open. |
| 2026-11-16 | Further partial D02/D03/D05: physical mixed-device contexts | This PR, context increment (configure-time parent `458e612`) | [Reproduction, fix, physical gates and hashes](CARTESIAN_SUBGRID_DOMAINS.md#physical-mixed-device-contexts). Actual RTX4090/RTX3080 pair, both orders, nine field/domain and six coupled-port cases per lane plus retained single-device/fields/ports/DC and CUDA CTest 24/24. Fields identical to one-device; maximum double error 1.941e-5, balance 1.706e-6. All 20 local stages, CTest 24/24/Python 14/14; 15 older prefixes and 11,255 arrays unchanged. ASan/UBSan domain-source controls pass. Host-only activation timing is not FDTD throughput evidence. Native distributed selectors/PML/losses/DFT, collective pause/resume, production scaling and antenna acceptance remain open. |
| 2026-11-16 | Partial extension/parity infrastructure | This PR, paired-advance increment (configure-time parent `c640333`) | [Boundary, tests, source/build and limits](CARTESIAN_SUBGRID.md#paired-compact-advance-extensions). All21 local stages, CTest25/25/Python14/14;16 older prefixes and11,255 arrays unchanged. Physical RTX4090/RTX3080 retained domains/ports plus CTest25/25; canonical real-drive CPU/reference/Metal/RTX4090/RTX3080 supplement separates the repaired test's earlier initial-only coverage. No conductivity/native lossy admission, distributed histories or antenna qualification. |
| 2026-11-16 | Further partial P01/U05: native source tails | This PR, waveform increment (parent b48e819) | [Contract and evidence](CARTESIAN_SUBGRID.md#native-source-waveform-tails). Step/custom sample-zero repair;21 local stages, CTest 26/26, Python 14/14 and17 older prefixes/11,255 arrays unchanged. Physical RTX4090/RTX3080 waveform/Step plus retained domain/port controls and CUDA CTest 26/26 pass. Native periodic/SSD, delayed/H/TFSF and antenna qualification remain open. |
| 2026-11-16 | Partial F01 and native G01/U04/U05: constant conductivity | This PR, conductivity increment (configure-time parent b48e819 plus qualified waveform patch) | [Contract, hashes and evidence](CARTESIAN_SUBGRID_CONDUCTIVITY.md). Independent dense and native analytic/priority/unit controls; all 25 local stages, CTest 29/29/Python 14/14; physical RTX4090/RTX3080 single-device loss/native controls and retained pair fields/ports, CUDA 13-gates/CTest 29/29 pass. 21 older prefixes and 11,255 arrays unchanged. PML/loss, skin/surface/ADE, distributed histories, scalability and antenna qualification remain open. |
| 2026-11-16 | Further partial F01/B02/B03/U04/U05: disjoint conductivity/PML | This PR, prepared-support increment (parent `ec73f8d`) | [Joint contract and build hashes](CARTESIAN_SUBGRID_CONDUCTIVITY.md#disjoint-conductivity-and-pml). All26 local stages, CTest30/30/Python14/14; all15 physical CUDA gates/CTest30/30. Twenty-four ledger cases/eight field comparisons per RTX4090/RTX3080, native scripting, zero-loss/order/overlap guards and retained lossless pair controls. Maximum joint field error7.431e-7 and balance1.805e-7. Twenty-two older prefixes/11,255 arrays unchanged. No overlapping/Mur loss, surface/ADE, distributed histories, representative PCB or scaling qualification. |
| 2026-11-16 | Partial F02, retained F01/B02/U04/U05 | This PR, surface-state increment (parent `08e2274`) | [Equations, source/build identities and limits](CARTESIAN_SUBGRID_SURFACE.md). All30 local stages/all19 physical CUDA gates, CTest33/33/Python14/14;11,255 arrays unchanged. CPU/reference/Metal and RTX4090/RTX3080 scalar/path surface controls, disjoint PML, unchanged retained lossless pair controls in both orders. Native conductor geometry, RF line/array, distributed histories and scaling remain open. |
| 2026-11-16 | Further partial F02/G01/P02/P05/U04/U05 | Native volume/path increment `c006caa` | [Contract, source/build identities and limits](CARTESIAN_SUBGRID_CONDUCTOR.md). Actual r03 source matches the published numerical/test files; all 34 local stages, CTest 37/37/Python 14/14, all 27 physical CUDA gates and CTest 37/37 pass. Native Cu18/Cu35 opaque volumes, same-metal unions/priority cuts and ordinary per-cell R/N paths run on CPU/reference/Metal and each RTX4090/RTX3080. Retained lossless pair controls pass both orders; 11,255 older arrays remain bit-identical. Actual older-CSXCAD build passes 35 solver tests with two explicit absent-feature skips. Point-port RF mesh changes remain unresolved; arbitrary CAD, scalable interiors/interfaces, distributed histories/DFT and antenna acceptance remain open. |
| 2026-11-16 | Partial C01/C03/U01–U03/U05, focused spatial interiors | This PR, exact-pair spatial-curl increment (parent `d026d92`) | [Source contract, controls and limits](CARTESIAN_SUBGRID_BLOCKS.md). Host-only r02:13 stages/CTest38/38/Python14/14 and11,255 exact prior arrays. Later focused device geometry passes on reference/Metal and each physical RTX4090/RTX3080; composed local replay passes14 stages/CTest39/39/Python14/14 with11,255 exact prior arrays. Full CUDA replay passes33 gates and CTest39/39 on the actual RTX4090/RTX3080 pair. Interleaved Metal layout timing is1.28–1.40x slower; no speedup is claimed. Expanded upload/reference field organization, distributed histories, general scaling and antenna acceptance remain open. |
| 2026-11-16 | Partial C01/C03/C05/U01–U03/U05; no guard increase | This PR, exact order and streamed preparation | R06 CUDA33 gates/CTest39; composed r07 local14 stages/CTest39/Python14 and11,255 unchanged arrays. Bounded491,520-leaf peak RSS−23.60%, setup+3.34%; r07 CUDA33 gates/CTest39 pass with immutable runtime/source and terminal exit0/no OOM. See [block evidence](CARTESIAN_SUBGRID_BLOCKS.md#exact-ordering-and-streamed-preparation). Local r07 core `4af61e6916978ba5d017b802e1bc64355cd6b3ce133d3753568e4dde896ece1f`.|
| 2026-11-16 | Further partial F02/G01/G04/G06/U04/U05, native shaped conductors | This PR, r08 (parent `549bf55`) | [Mapping, evidence and limits](CARTESIAN_SUBGRID_CONDUCTOR.md#native-shaped-conductors). CPU/Metal15 stages/CTest40/Python14; physical CUDA35 gates/CTest40;11,255 prior arrays unchanged. Exact occupied-box decomposition, hollow-bore and positive loss/energy controls; six scoped older-dependency tests with three explicit skips. Local core `5e9b2a6d385ecb7d98c8ba306dd63e2fe139475c1f79685716dd2c2dbe7e2203`; CUDA `69b8fa61dead03bf134f3aaf815dffcd997084bbabdd8529d42faafade2386a5`. Conformal/curved RF, distributed histories and antenna acceptance remain open. |
| 2026-11-16 | Further partial F02/G01/U05: native closed-ring admission | This PR, r09 (parent `56be081`) | [Conductor contract](CARTESIAN_SUBGRID_CONDUCTOR.md#native-shaped-conductors). Eight positive/six negative shapes; CPU/Metal15 stages/CTest40/Python14; physical CUDA35 gates/CTest40;11,255 arrays unchanged. Local core `e7f9501ef239d13f657e869fd9a04e6dade05cb9b40a5f47a03b936a3d995955`; CUDA `5f273237e245d42dcde5340fc6be837b6b4cb1ad36fa4c43c1dd09358f4f60bd`. Existing CAD retained; RF/antenna/distributed qualification stays open. |
| 2026-11-16 | Diagnostic only; F02/A02/A03 remain open | Unchanged qualified r09 runtime | [Coarse curved-boundary diagnosis](CARTESIAN_SUBGRID_CONDUCTOR.md#curved-line-boundary-diagnosis). Native reconstruction checks126,390 leaves with zero mask mismatches. High conductivity approaches matching PEC mask within1.52e-5 complex S11/0.00335° S21; boundary classification reproduces12.627°. Consistent incremental Cu transmission differs0.02953°/0.00258dB. One coarse fixture only; no curved convergence, upstream-defect or antenna acceptance claim. |
| 2026-11-16 | Further partial D01/D02/F01/F02: cold C++ constitutive ownership | This PR, r13 (parent `9eb7291`) | [Prepared conductivity domains](CARTESIAN_SUBGRID_DOMAINS.md#prepared-conductivity-domains). Reference14/Metal15 constitutive cases; local19 stages/CTest41/Python14 and11,255 exact prior arrays. CUDA39 gates/CTest41, including58 constitutive cases on each physical RTX4090/RTX3080 and both pair orders. FP32 one-device fields are identical; maximum double field error5.830e-6 and balance2.638e-7. Local core `f355d86377b742520112925f8f38cd0dfeed66c3724510a477254b865f6bf880`; CUDA `21ab8476d5ba7845e6518bdc00743cbe3e947ffa57ea851872eb3cfc8a64c83f`. Native distributed PML/DFT/pause, RF/antenna and scaling remain open. |
| 2026-11-16 | Further partial G01/G04/F02: internal opaque refinement faces | This PR, r14 (parent `daca3b6`) | [Internal-face qualification](CARTESIAN_SUBGRID_CONDUCTOR.md#internal-opaque-refinement-faces). Twelve positive/all-normal face controls plus retained partial-face rejection; the new test fails against the actual R13 parent after two controls. Local19 stages/CTest42/Python14 and11,255 identical prior arrays pass. Contained native Cu18 setup passes Metal and eight actual CUDA steps; full CUDA39-gate/CTest42 compatibility replay passes with58 constitutive comparisons. Finer curved RF remains in progress, not accepted. No geometry/trace/update-equation change or production adoption. |
| 2026-11-16 | Further partial D01/D02/B02/U04: cold PML composition | This PR, r16 candidate (numerical parent `4663326`) | [Prepared PML domains](CARTESIAN_SUBGRID_DOMAINS.md#prepared-pml-domains). Six reference/seven Metal cases and local43/43 CTest pass; exact one-device FP32 fields, maximum double field error1.420e-5 and work/balance below4.387e-8. Full local20-stage replay passes with11,255 exact arrays; physical R16 replay reproduces the context failure corrected below. Native distributed scripting/DFT/pause, general PML parity, RF/antenna and scaling remain open. |
| 2026-11-16 | Partial D/B/U; compact PML owner context | R17 follow-up to `daa8596` | Reproduced R16 physical-pair failure; [own-device fix](CARTESIAN_SUBGRID_DOMAINS.md#compact-pml-owning-device-correction),20 local stages/Ctest43/Python14/11,255 exact arrays and both five-case physical PML pair orders. Physical43-gate/CTest43 replay passes, including24 PML and58 constitutive comparisons; no native distributed/antenna adoption. |
| 2026-11-16 | Further partial D/P/B/U: native fixed local domains | R18 adapter and R02 sentinel follow-up | [Native runtime contract](CARTESIAN_SUBGRID_RUNTIME.md): actual CLI/Python, reference/Metal, PML/loss, TD/FD/U/I and scoped dipole NF2FF; R01 full20-stage/CTest44/Python14/11,255 exact replay, R02 focused CTest44/native controls. R18-r01 physical48-gate/CTest44 and16-case native Python matrix passes; parser-only R02 physical follow-up passes7 gates/CTest44 and actual CLI/native Python controls. No representative-array/scaling or external collective-pause qualification. |

| 2026-11-16 | Further partial C01/C03/C05/D01/D05: global domain row stream | R19 over `389d372` | [Setup contract and measurements](CARTESIAN_SUBGRID_DOMAINS.md#streamed-global-domain-preparation). Reference/Metal12 exact metadata/field cases, native8/16 TD/FD/port cases and CTest45/45 pass. Actual frozen-parent491,520-leaf peak RSS−23.68%, whole setup+10.07%; local core `560fa47d714b125f0f9417f06692969422e64745da0eef732055033bef8d7fbb`. Physical CUDA pending; local uploads/topology/scaling and antenna acceptance remain open. |

| 2026-11-16 | Qualification audit and test-fixture correction | R19 fixed test over `ef5a17b` | The original local test binary lacks the late port edit; the original CUDA fixture rejects resistor-after-source setup before backend creation. Preserve both records. A corrected separately compiled fixture passes6 reference/6 Metal cases against unchanged core `560fa47d...`; this does not accept the failed CUDA run or qualify an antenna. |

| 2026-11-16 | Further partial C01/C03/C05/D01/D05: exact local row packages | R20-r03 over `21eff5b` | [Local package and measurements](CARTESIAN_SUBGRID_DOMAINS.md#streamed-local-domain-preparation). Fresh source-before-build reference/Metal and native8/16-case controls; CTest45/45 and13 rejection guards. Actual R19 comparison: peak RSS−30.37%, whole setup+7.93%; graded benefit2.88%.960k setup-only control passes. Local core `65309099920fd3f93e5e039ffde808ed77fb520910d9f1c17d672f1ce7de39c7`; physical CUDA pending, no antenna acceptance. |

| 2026-11-17 | Further partial C03/C05/D01/D05: shared owned-state domain addresses | R21, parent `648683d` | [Address/storage evidence](CARTESIAN_SUBGRID_DOMAINS.md#shared-owned-state-domain-addresses): focused sanitized map oracle, exact reference/Metal metadata/fields, native8/16 and CTest46 pass. Interleaved RSS/setup gains only; physical14-gate/CTest46 replay passes on4090/3080 in both orders, with native observations. No guard increase, antenna or production qualification. |

| 2026-11-17 | Further partial C03/C05/U04/D05: lazy CPU reference workspace | R22, parent `ef323a2` | [Workspace evidence](CARTESIAN_SUBGRID_DOMAINS.md#lazy-serial-reference-workspace): first-use/negative/capacity controls, reference/Metal exact metadata/fields, native8/16 and CTest47 pass. Interleaved cold-storage/RSS reductions only; physical14-gate/CTest47 replay passes in both pair orders after terminal review. No guard increase or production/antenna qualification. |
