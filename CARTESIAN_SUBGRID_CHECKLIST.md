# Cartesian subgridding — implementation and qualification checklist

Last reviewed: **17 November 2026**.

This is the canonical detailed checklist for extending the Cartesian subgridding prototype
to the functionality of the regular Cartesian solver. Check an item only after
its stated evidence exists. Implementation and numerical qualification are
separate tasks; a passing test on one engine does not qualify another engine.

## Current RF qualification correction —17 November2026

The previous R47 ordinary-Yee spatial-reference acceptance is **withdrawn
pending a safe-step replay**. R64's stricter preparation caught a forced2fs
step above the refined ordinary calculated1.1785fs bound. Audit then found
three preserved R47 ordinary logs warning that4/5fs exceeded their calculated
2.88481fs bound. Their finite raw results remain observations, not accepted
spatial-convergence evidence. The earlier review omitted this warning check;
this is a project qualification defect, not an established engine defect.

R62's native single/mixed-pair replay still establishes scoped runtime/domain
agreement on the unchanged native input. It does not validate those ordinary
references or asymptotic accuracy. R63 CPU77/physical82/host-memory evidence
is independent and remains accepted. R64-r02 retains all geometry/old knots/Z
with a common1fs candidate, strict per-case native CFL rejection, and new
setup/RF/time gates. No frozen data or numerical tolerances are replaced.
Evidence: `curved-xy-refinement-r64/prior-cfl-audit-r01.json`,
`curved-xy-refinement-r64/preparation-r02.json`.

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
16,384 base intervals per axis, default one-million base/total leaf counts,
and 16,384 coupled degrees of freedom per interface factor. An explicit checked
leaf-count budget can opt into larger counts; R50 characterizes15M homogeneous
host leaves without qualifying antenna capacity or changing defaults. Above1,024 edges,
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

A full-board input audit now confirms an additional native-port blocker:
rectangular sheet resistors with transverse extent are rejected by the current
line-path mapper. A minimal line-versus-sheet fixture reproduces the exact
rejection on the unchanged R32 runtime. Native multi-edge line support does
not qualify that distributed geometry. P04/P05 and representative-antenna
admission remain open until the ordinary geometry-normalized per-edge G,
cap electrodes, probe conventions, passive power and device/domain execution
are independently qualified. Do not collapse the antenna port to a line or
increase the64-resistor dense reference guard as a workaround.

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

Current C01/C03 preparation increment: [implicit leaf addresses](CARTESIAN_SUBGRID_BLOCKS.md#implicit-rootchild-leaf-addresses)
is implemented with Linux CPU48/48, exact operator/field and physical
4090/3080 controls, including a separately corrected native16-case temporary
directory followup. Local Metal qualification is pending after two
host-reserve build stops. Peak setup RSS falls1.60–3.41% in three bounded
interleaved cases, with longer setup. C01/C03/C05 remain incomplete; no
antenna or guard-increase acceptance follows. The direct-row R25
candidate subsequently passes CPU49/49, sixteen exact operator/field files,
twelve coefficient controls and native reference cases, plus thirteen physical
CUDA gates with24 exact domain cases and16 native cases in both pair orders.
Interleaved peak setup RSS falls18.64% on the graded case but grows
3.25–3.79% on regular cases, so it is not adopted as a general memory win.
R26 compresses only exact incoming face IDs, with
an allocation gate retaining raw CSR when smaller. Its Linux CPU50/50,
sixteen exact frozen-parent files and native8-case controls pass. Three
interleaved repetitions against R24 reduce peak RSS7.59%/22.63%/8.20% on
uniform491,520/graded61,440/uniform960,000 leaves, with whole setup
changes+10.92%/−11.92%/+2.39%. Thirteen physical CUDA gates pass on each GPU
and in both pair orders, including24 exact domain and16 native Python cases.
The exact R26 numerical/test bytes are now in this branch. Metal remains open;
this is not a general speedup or full-PCB adoption. Exact native child
coordinates and API/XML from R27-r02/R28 are now included with matching source
and test bytes; see [native coordinate contract](CARTESIAN_SUBGRID_RUNTIME.md#exact-native-child-coordinates).
Linux CPU53/53,16 exact frozen-parent files, eight native reference cases,
18 physical CUDA gates/48 new coordinate cases,24 exact domain cases and16
native Python cases pass on4090/3080 in both orders. Metal and actual
Matlab/Octave remain unqualified. An additional standalone R30 composition
fixture passes14 CPU/reference native-coordinate bulk/surface-loss/port/copper
cases plus58 physical CUDA cases across each GPU and both orders, with exact
single/domain fields and maximum double-field error2.033e-6. Its unchanged
standalone test is now supplied and registered; the combined registration passes in R31's56-test superset, while Metal for
the composition remains open. The original
degenerate field diagnostic and two R29 fixture failures remain retained.
R29-r02 implements implicit incoming curl support without global incoming CSR.
Its corrected24-case all-row oracle and full CPU54/54 pass, including16 exact
actual-R28 files, eight native reference and six exact domain cases. Three
interleaved repetitions reduce peak RSS2.28–3.82% against R28, with slower
whole setup by8.84–18.90%; this is a structural memory/overflow fix, not a
preparation speedup. Eighteen physical CUDA gates pass on4090/3080 in both
orders, including24 exact domain and16 native Python cases. The exact
numerical/test source is now promoted; Metal remains open. See the
[implicit incoming contract](CARTESIAN_SUBGRID_DOMAINS.md#implicit-incoming-curl-support).
The qualified R31 source now compresses only
worker ownership labels in256-state uniform/byte-mixed tiles, retaining the
dense oracle. Thirty storage/address controls, full CPU56/56,16 exact parent
files, native reference and the unified R30 test pass. Three interleaved
repetitions reduce peak RSS2.90–5.41% against R29, with measured whole setup
+2.36–4.73%. All22 physical CUDA gates pass on4090/3080 in both orders,
including24 exact domain cases,16 native Python cases and58 native-child loss
compositions. Exact source is promoted; Metal remains open. See the
[tiled-owner contract](CARTESIAN_SUBGRID_DOMAINS.md#tiled-domain-ownership-labels).
Qualified R32 metadata-only initial
state preparation removes global FP32 E/H copies while preserving global FTZ
qualification and constructor-local borrowing. Its original misplaced-flag
compile failure is retained; R32-r02 corrects only the header placement and
passes16 metadata/tail/domain-state controls, full CPU57/57,16 exact
actual-parent files and native observations. Three interleaved measurements
reduce peak RSS2.73–6.05% against R31, with whole setup+0.23–3.08%;
all22 physical CUDA gates pass on4090/3080 in both orders. Exact source is
promoted; Metal remains open. See the [initial-state contract](CARTESIAN_SUBGRID_DOMAINS.md#metadata-only-initial-state-preparation).
Qualified R33 owns finalized labels and precomputes canonical owned ranks using
uniform-tile bases or exact mixed/dense uint32 ranks, removing mutable global
address slots. Its frozen implementation builds successfully. Fifty rank size/pattern
cases, five independent map comparisons, ASan/UBSan, full CPU58/58,16 exact
actual-parent files and native observations pass. A separately named collector
corrects only a copied19-versus20 wrapper-command count on the same archive;
no numerical rerun. Three interleaved repetitions reduce RSS2.43–5.98% against
R32, with whole setup−1.04% to+1.67%. All22 physical gates pass on4090/3080
in both orders, including24 exact domain cases and16 native Python cases.
Exact source is promoted; Metal remains open. See the [owned-rank contract](CARTESIAN_SUBGRID_DOMAINS.md#immutable-owned-local-ranks).
A separately frozen R34 native
rectangular-port candidate composes geometry-normalized independent G through
the existing admittance extension, preserving cap electrodes and the64-port
dense guard. Its original new-test compile error is preserved; r02 changes
only that test declaration. Its next retained test fixture violates the required
region margin; separately corrected fixtures preserve the margin and compare
ordinary fields on matching grids, ordinary G on both layouts and actual
refined CPU/device fields. Fresh r03 also keeps the shared timestep below the
ordinary CFL. Full CPU59/59,48 new port cases,16 exact actual-parent files
and native observations pass. A separately named collector corrects only a
renamed retained rank-log marker on the same immutable archive; no numerical
rerun. The26 physical gates, Metal, native RF and full-antenna acceptance remain
open. Separate R35 retains already verified local-edge support to avoid a
whole-grid voltage-path scan per port edge. Its original new test assumed
identical support; instrumentation finds a tiny extra general-quadrature
coefficient at a rounded cut (3.63e-16 relative), with identical main weight.
Fresh r02 changes only that oracle, keeps the original2e-14 numerical bound
on the complete functional and adds wrong-main/wrong-support negatives. Its
first48 path cases pass; full CPU60 and physical qualification remain open.
No R34/R35 numerical source is promoted.
No full-PCB, guard increase or production adoption follows.

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
passes the same local replay and five physical 4090/3080 PML controls in each
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

The separate [compact two-resolution U diagnosis](CARTESIAN_SUBGRID_CONDUCTOR.md#compact-curved-line-two-resolution-diagnosis)
completes eight physical CUDA cases with unchanged predeclared gates. At
16/16/1, local/ordinary copper complex-S11 difference falls to0.006944;
the8/8/1 copper case still misses0.01. Both engines change phase by about
3.4deg across those XY levels, so spatial convergence remains unqualified
despite negligible recorded late-time tails. F02/A01–A07 remain open.

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

Native rectangular resistor preparation is separately qualified on CPU/CUDA;
see the R34 completion row and [contract](CARTESIAN_SUBGRID_RUNTIME.md#native-rectangular-resistive-ports).
Driven RF and representative PCB qualification remain open.

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

| 2026-11-17 | Further partial C01/C03/C05: implicit leaf geometry | R24-r03, parent `c003761` | [Address/snapshot contract and evidence](CARTESIAN_SUBGRID_BLOCKS.md#implicit-rootchild-leaf-addresses). CPU48/48, sixteen exact operator/field files,126 address cases/723,354 cells,128 snapshot controls, twelve physical CUDA gates plus a separate corrected native16-case followup. Original TMPDIR collision and rejected R01/R02 prototypes retained. Peak RSS−1.60–3.41% with longer setup. Metal pending after resource stops; no scaling/antenna acceptance or guard increase. |
| 2026-11-17 | Further partial C01/C03/C05/D01/D05: direct rows and integer incoming addresses | R26, following unpromoted R25 measurement | [Domain preparation and three-prefix measurements](CARTESIAN_SUBGRID_DOMAINS.md#direct-canonical-rows-and-incoming-face-addresses). CPU50/50,16 exact frozen-parent files,12 coefficient controls,13 physical CUDA gates/24 exact domain cases and16 native cases in both pair orders. Peak RSS−7.59–22.63% against R24; whole setup+10.92/−11.92/+2.39%. Exact numerical/test bytes promoted; Metal, full-array scaling and antenna RF remain open. |

| 2026-11-17 | Further partial M/G/P/U/D: exact native child coordinates | R27-r02/R28; exact qualified source promoted | [API and coordinate contract](CARTESIAN_SUBGRID_RUNTIME.md#exact-native-child-coordinates). CPU53/53,16 exact frozen-parent files,18 physical CUDA gates/48 child-coordinate cases,24 exact domain cases and16 native Python cases pass on4090/3080 in both orders. Metal/Octave, large-scale and antenna qualification remain open; no guard increase. Separate R30 native-coordinate loss composition passes14 CPU/reference and58 physical CUDA cases; combined new CTest registration and Metal remain pending. |

| 2026-11-17 | Further partial C01/C03/C05/D01/D05: implicit incoming curl support | R29-r02, exact qualified source promoted | [Topology query and measured tradeoff](CARTESIAN_SUBGRID_DOMAINS.md#implicit-incoming-curl-support). CPU54/54,16 exact actual-R28 files,24 all-row controls and18 physical CUDA gates/24 exact domain cases/16 native cases pass. Interleaved peak RSS−2.28–3.82%, whole setup+8.84–18.90% against R28. Metal, large-array and RF convergence remain open; guards unchanged. |

| 2026-11-17 | Further partial C03/C05/D01/D05: tiled domain ownership | R31, exact qualified source promoted | [Contract and measured storage](CARTESIAN_SUBGRID_DOMAINS.md#tiled-domain-ownership-labels). CPU56/56,30 storage/address controls,16 exact actual-parent files and22 physical CUDA gates pass on4090/3080 in both orders. Native child losses, TD/FD/ports and exact domain fields are retained. Interleaved peak RSS−2.90–5.41% against R29; whole setup+2.36–4.73%. Metal, full-array scaling and RF remain open; allocation guards unchanged. |

| 2026-11-17 | Further partial C03/C05/D01/D05: metadata-only initial state | R32-r02, exact qualified source promoted | [Contract and preparation measurements](CARTESIAN_SUBGRID_DOMAINS.md#metadata-only-initial-state-preparation). CPU57/57,16 exact actual-R31 files,16 initial metadata/tail controls and22 physical CUDA gates pass on4090/3080 in both orders. Global FTZ bounds and local owned/ghost initial values remain exact. Interleaved peak RSS−2.73–6.05% against R31, whole setup+0.23–3.08%. Earlier original compile failure retained; Metal, full-array scaling and RF remain open. |

| 2026-11-17 | Further partial C03/C05/D01/D05: immutable local owned ranks | R33, exact qualified source promoted | [Owned-index contract and storage](CARTESIAN_SUBGRID_DOMAINS.md#immutable-owned-local-ranks). CPU58/58,50 rank cases/five map oracles, ASan/UBSan,16 exact actual-R32 files and22 physical CUDA gates on4090/3080 in both orders. Interleaved peak RSS−2.43–5.98% against R32; whole setup−1.04% to+1.67%. Collector-only19/20-count failure preserved and corrected against the same archive without rerunning CPU. Metal, full-array scaling and RF remain open. |

| 2026-11-17 | Further partial P01/P02/P05/G03/U04/U05: native rectangular resistors | R34-r03, exact qualified source promoted | [Port contract and evidence](CARTESIAN_SUBGRID_RUNTIME.md#native-rectangular-resistive-ports). CPU59/59,16 exact frozen-R33 files and26 physical CUDA gates pass on4090/3080 in both orders. New48-case geometry/coefficient/passive-field controls per lane retain native sheet/volume bounds and optional PEC end caps; retained domains/loss/PML/native observations pass. Earlier test/collector failures preserved without numerical reruns. Metal, driven RF, scaling and antenna acceptance remain open; guards unchanged. |

| 2026-11-17 | Further partial P01/P02/P05/C05/U04/U05: resolved native port support | R35-r02, exact qualified source promoted | [Descriptor contract, diagnosis and measurements](CARTESIAN_SUBGRID_RUNTIME.md#reusing-resolved-native-port-support). CPU60/60,48 descriptor controls,16 exact actual-R34 files and26 physical CUDA gates on4090/3080 in both orders pass. Original support-cardinality test failure retained; unchanged2e-14 complete-functional bound and wrong-support negatives pass. Actual native294-trace setup improves12.83x/13.24x at16,192/118,720leaves with identical prepared fingerprints. This is scoped setup timing, not full-array/FDTD speedup. Metal, RF and scaling remain open; guards unchanged. |

| 2026-11-17 | Further partial P01/P02/P05: graded native line normalization | R36, exact qualified source promoted | [Physical-length contract](CARTESIAN_SUBGRID_RUNTIME.md#physical-length-native-line-normalization). CPU61/61,16 exact actual-R35 files and28 physical CUDA gates pass; the new12-case ordinary-coefficient oracle fails against the actual old library. Both4090/3080 orders,24 exact domains, losses/PML and16 native Python cases are retained. Metal, driven RF and representative antenna qualification remain open. |

| 2026-11-17 | Further partial C03/C05/D01/D05: sequential worker uploads | R40, exact qualified source promoted | [Lifetime contract and measurements](CARTESIAN_SUBGRID_DOMAINS.md#sequential-numeric-worker-upload). CPU62/62,16 exact actual-R36 files and32 physical CUDA gates pass;24 new exact upload/lifetime cases, both pair orders, native observations and losses/PML remain. GPU-free preparation peak RSS falls11.04–18.59% in three interleaved controls. Metal, end-to-end native device RAM, full-array scaling and antenna RF remain open; no guard increase. |

| 2026-11-17 | Further partial C03/C05/U05/D05: paired magnetic metrics | R42, exact qualified source promoted | [Contract and measurements](CARTESIAN_SUBGRID_DOMAINS.md#exact-paired-magnetic-metric-storage). CPU62/62,16 exact actual-R40 operator/field files and32 physical CUDA gates pass in both4090/3080 orders. All six field/history states remain distinct. Four-mode/four-repeat preparation RSS falls a further2.47–5.22% against consumed R40; Metal, steady device RAM/full-array scaling and RF remain open. No guard or production-runtime change. |

### Pending port follow-ups — 17 November

R35-r02 resolved local port descriptors pass CPU60/60,16 exact actual-parent
files and all26 physical CUDA gates after terminal/source/runtime review. An actual native-port
setup benchmark on read-only R34/R35 libraries times complete preparation,
with294 traces and identical prepared-coefficient fingerprints:16,192leaves
2.803563s→0.218491s and118,720leaves20.918891s→1.580225s (three interleaved
repetitions). These12.83x/13.24x small-case setup ratios do not establish
full-array or FDTD speedup. Exact qualified R35 source is now promoted.
Metal and representative RF/scaling remain open.

The [retained line-adapter diagnosis](CARTESIAN_SUBGRID_RUNTIME.md#retained-native-line-normalization-diagnosis)
reproduces graded longitudinal R/N normalization error against actual ordinary
coefficients. The physical-length correction R36 now passes CPU61/61 and28
physical CUDA gates, with16 exact actual-parent files and both4090/3080 orders.
Exact qualified source is promoted; Metal and representative RF remain open.
No frozen result, production runtime or allocation guard changed.

R37 prepares four native finite-width Cu18/PEC U-line models at8/8/1 and
16/16/1, preserving the previous copper, grids and2ns cap while deliberately
replacing point launches with full-strip-width sheets and broad fixed current
contours. All four native setups pass without FDTD; expected CLI setup-only
exit1 is checked with completed-operator diagnostics, not accepted by itself.
All eight local/ordinary RF controls complete on one authorized4090 using
actual accepted R35 and frozen ordinary PR242. Independent raw/timestamp DFT
review finds both PEC levels and fine copper pass the original band gates.
Coarse copper misses complexGamma (0.015748), as does local copper two-level
sensitivity (0.013325). Time-prefix tails are negligible. The full RF gate
remains failed; fixed Z and two XY levels do not qualify 3D or antenna convergence.

### Current scaling and RF follow-ups

R37 finite-width U-line controls improve the former point-launch sensitivity,
but remain a failed full RF qualification: coarse copper/local two-level
complexGamma exceeds0.01. Fine copper agrees with ordinary PR242 within
0.001784 complexGamma/0.16851°/0.00443dB over24–24.25GHz. All time-prefix gates
pass; no duration extension, failure erasure or threshold relaxation follows.

R38 profiles the unchanged R35 library's C++ allocations on three bounded
cases. A960k-leaf grid plus three-domain plan requests about528MB; requested
capacity includes uncommitted pages and is not measured antenna RAM. R39
tests cold capacity reduction with an interleaved resident-memory benchmark;
CPU61/61 and16 exact parent files pass, but three-repeat RSS changes against
R36 are−0.26%,−0.21%,+0.34%; no useful memory reduction is established.
R39 remains unpromoted. R40 separately implements
final topology followed by one-domain prepare/upload/release, with retained
C++ plans as the default independent path. CPU62/62,16 exact parent files
and32 physical CUDA gates pass; exact source is promoted. Three-repeat
GPU-free preparation peak RSS falls11.04–18.59% against actual R36. Metal,
end-to-end native device RAM, full-array scaling and RF remain open; there
is no allocation-guard increase or production adoption.

R41 prepares fixed-fine-XY Z sensitivity without increasing guards: eight
dielectric plus two physical Cu interval midplanes raise587,384 to906,304
local leaves. Both native PEC/Cu setups pass on accepted R35, with unchanged
geometry, ports, current contours, XY,5fs timestep and2ns cap. All four new
local/ordinary RF cases complete, exit0/noOOM, with independently rechecked
raw/input/runtime and own-timestamp DFT. Same-level Z2 engine comparisons
pass (copper0.001172 complexGamma/0.14678degree/0.000736dB), but Z1-to-Z2
complexGamma0.011015–0.011168 fails the0.01 sensitivity gate in both engines
and both conductors; late-time tails pass. RF convergence remains open. This is neither an antenna launch nor asymptotic3D
convergence. Frozen R37 ratio16 is the separately hash-bound Z1 comparison.

R42 was separately frozen over R40 and is now promoted after the CPU/physical
qualification recorded below: identical opposite-face metric coefficients are stored once per pair while all six magnetic
field/history states remain separate. Its build, independent face-view
oracle, CPU62/62 and16 exact actual-R40 operator/field files pass. A four-mode
resident-memory comparison and physical replay remain pending. It changes no equations,
precision, mesh or size guards.

R43 full-model type inventory identifies a separate required PCB-CAD gap:
rounded laminate/backing volumes and cylindrical air bores are outside the
native material adapter's strict box-only contract. Small unchanged-runtime
box-positive/shaped-negative controls reproduce that explicit feature limit.
Do not remove these shapes or silently use bounding boxes to admit an antenna.
Shaped material integration, priority and convergence must be implemented and
qualified separately; memory improvements alone do not make the model usable.

R43 independent archive/runtime/log review confirms the explicit shaped-material
limit: one tiny box setup completes, extruded material polygon and priority air
bore reject, with zero FDTD steps. The failed original XML-writing harness is
retained separately; the full antenna was only inventoried, not launched.

R44 is a separately frozen unpromoted candidate over R40, excluding R42.
Optional material samples2/4/8 select positive tensor-midpoint volume averages
of epsilon/mu and electric/magnetic conductivity, resolving native CAD priority
at each point. Strict aligned boxes remain the default. Original contours are
retained; conservative contour-segment bounds limit quadrature to boundary
leaves. This is a scalar approximation, not conformal integration. Independent
analytic fraction/energy and six native axis/unit/priority controls plus
C++/Python/XML/Matlab interfaces are prepared; build and qualifications remain
pending. Thin inclusions, quadrature/mesh RF convergence and actual shaped-
material device/domain controls remain open. No full-array admission follows.

R42's four-mode setup benchmark completes56 commands with four interleaved
repetitions, retaining exact counts/curls/timestep and24 fewer logical metric
bytes per leaf. Measured additional peakRSS reduction against consumed R40 is
2.47–5.22%; native GPU steady RAM/FDTD/full-array savings are not established.
Physical32 replay is now running after R41's actual terminal hardware release.

R44's C++ candidate builds; the original Cython wrapper fails on the imported
libcpp.bool shadowing the Python builtin and is preserved. Fresh R02 corrects
only that wrapper. Its new analytic/native test detects a test-oracle constant
mismatch even in air leaf0: hardcoded permeability differs from native MUE0 by
4.48883e-10, above the unchanged2e-13 comparison. A diagnostic against unchanged
R02 core reproduces zero epsilon error and exactly that mu discrepancy. R03
uses native EPS0/MUE0 in the independent analytic CAD oracle, retaining the
bound and numerical code, and is separately frozen/building. This does not
qualify shaped-material physics, RF, CUDA/Metal or representative PCB use.

R44 current scoped evidence: R03 builds successfully and its six native
polygon/priority-air fixtures pass across all axes, mm/micrometre units and
reversed contours/extrusions. A small octagon-minus-bore volume differs0.051052%
from analytic CAD;2000-step energy drift is1.55431e-15. R04 corrects only the
external XML-test fixture's missing excitation/CSX and both API tests pass on
unchanged R03 library/bindings. Full same-build CPU 63/63 and 16 exact R40
operator/field files now pass, independently reviewed after terminal exit0/noOOM.

R06 independently reviewed reference-device controls pass12 native material
cases (six geometries on one/three logical workers),256steps each. FP32 fields
are identical; maximum normalized difference from double is5.79308e-7 under
retained2e-4/2e-6 double/FP32 gates. Test-only missing-header/oracle/API failures
are preserved. R07 completes 36 material-field cases plus 32 retained same-build
CUDA gates after accepted CPU63 and actual R42 physical release. Independent
review checks actual backend models, both pair orders, raw metrics and 763
source hashes. FP32 fields agree exactly; the standalone R44 source remained unpromoted until
the separately qualified combined R46 integration below.
Metal, mesh/RF, resource/scaling and representative antenna acceptance remain open. The production runtime and frozen antenna inputs remain unchanged.

R42 physical qualification is complete:32 gates, native16 Python cases,
retained loss/PML/ports/fields and upload lifetime controls pass, terminal
exit0/noOOM on actual4090/3080 pair. Exact source is promoted after independent
759-source/385-result rehashing. Remaining Metal/scaling/RF status is unchanged.

R45 is a separately frozen unpromoted candidate over R42 (R44 excluded).
It adds an explicit count-only budget across C++/native/Python/XML/Matlab;
default construction retains one million leaves. Checked products/region sums
and a conservative32-bit edge-address ceiling reject before overflow/allocation.
The index ceiling is not a RAM guarantee or qualified problem size. The isolated
build completes exit 0/no OOM. Small boundary/overflow/unchanged-field/API
controls and full CPU 63/63 replay pass, with 16 exact R42 operator/field
files. Externally bounded 2M/4M characterization completes all six cases and
is independently reviewed. No full antenna is admitted and
interface-factor limits remain unchanged.


R44 external R10 sampled-loss composition passes 36 reference cases: three
extrusion axes, two drawing-unit/contour directions and electric/magnetic/both
constant conductivities, on one and three logical workers. An independent
analytic octagon-minus-bore midpoint sum, without CSXCAD or the production
quadrature helper, gives identical native/manual fields and loss ledgers.
Across 512 steps/case, maximum CPU energy imbalance is 1.28750e-14 and reference
FP32 field difference from double is 1.48051e-6 under unchanged 2e-4/2e-6 gates.
R09's external include-order/PI macro compile failure is retained; R10 changes
only test includes, not the built R03 library. Physical/Metal sampled-loss,
material-interface RF convergence and full-antenna qualification remain open.


R44 external R11 completes 108 sampled bulk-loss cases: reference single/three
workers, each physical 4090/3080 and both physical pair orders. The independent
analytic native/manual CPU fields and ledgers agree exactly; physical/reference
FP32 fields also agree exactly. Maximum normalized FP32 loss-ledger difference
is 2.13088e-14 and energy imbalance 1.79962e-7, under unchanged bounds. Actual
terminal exit0/noOOM, GPU UUID/model observations, source/core/binary identities
and raw reviews are retained in the R44 study. These are small scalar-volume
bulk-loss controls, not conductor-surface, RF or full-antenna acceptance.

R46 is a separately frozen, now CPU/CUDA-qualified integration: R42 metrics,
R44 native shaped-material/bulk-loss support and R45 explicit count budget.
The previous external R04/R07/R11 test corrections/controls are included as
normal source, with a new combined API/XML/native setup regression. Both options
retain their defaults and independent reset behavior; no production/frozen
runtime changes. Its isolated build and CPU65/65 complete exit0/noOOM with16 exact frozen-R42
operator/field files. All40 physical CUDA gates pass on4090/3080 in both orders,
including sampled bulk losses, child knots and native TD/FD/ports/PML. Source/
runtime and396 archived files were independently rehashed. Exact numerical/API/
test source is now promoted; prior canonical evidence sections are preserved.
Metal, native PCB resource scaling and representative RF remain open.


R45 staged host characterization completes six GPU-free homogeneous-grid cases
at960,000 /2,109,120 /3,932,160 leaves, two curl layouts and four consumed worker
packages. Independent raw/source/runtime review confirms default rejection above
one million, exact owned counts,64 sampled curl rows and12 sampled E/H points
per layout pair. Maximum first-step energy drift is4.66294e-15. Peak process RSS
at3,932,160 leaves is1,710.26MiB for patterns and1,630.13MiB for spatial blocks,
under the unchanged6GiB external cap. Single observations do not establish
statistical performance, long-time stability, device memory or native PCB/RF
scaling. Larger/full-array admission remains open; do not extrapolate those
numbers into accepted antenna resource requirements.


R47 retains the failed native5fs Z3 copper setup: its requested timestep exceeds
its coupled bound. A separate GPU-free setup-only diagnostic reports half-bound
2.26377fs, motivating explicit4fs copper controls in both engines and both Z
levels, with unchanged2ns duration. PEC/runtime5fs controls remain. The named
R04 derivative must pass actual native setup before RF; no stability bound or
RF acceptance threshold is weakened. This is not an RF acceptance result.

R48 adds selected normalized-state readback and lazy complete-field host mirrors
as a separately frozen, now CPU/CUDA-qualified increment over R46. Its original port test
used an independent FP32 oracle without resistor support; that same mismatch
is reproduced on unchanged R46. Corrected tests retain source-only2e-6 and
single/domain2e-6 gates, using the independent double reference for port solves
under the established2e-4 bound. Original failure and test-only correction are
retained. Fresh CPU66/66 and44 physical CUDA gates pass exit0/noOOM, including each GPU
and both pair orders, retained native/loss/PML and24 new sparse-state cases.
Source/runtime and401 archived files are independently rehashed before exact
numerical/API/test promotion. Metal remains open. Native processing still requests complete fields;
no native DFT/RAM saving or representative antenna acceptance is claimed.


R49's exact physical E/H/U/I bridge is now promoted after CPU67/67,16 exact
frozen-R48 operator/field files and48 physical CUDA gates on4090/3080 and both
pair orders. Direct reconstruction and selected/full single/domain values agree
exactly; maximum relative independent-double difference is3.23037e-7. All773
source and406 archived files were independently rehashed. Native Processing,
energy and DFT still read complete fields; Metal/native savings remain open.
Evidence: `physical-observations-r49/`.

R49's external owner inspected an old completed container instead of R47.
Small physical controls consequently overlapped our own RF series. Fresh VRAM/
host guards still passed; independent numerical buffers were retained, with no
foreign preemption. No performance/timing-isolation claim is accepted. Original
owner/admission/audit remain immutable. The R51 successor now verifies exact
container Name+ID before admission and again remotely; regression rejects old
completed IDs, reused names, live/failed/OOM/dead predecessors.

R50 completes all six GPU-free homogeneous controls at3.93/7.68/15M leaves on
immutable R48-r04, under4CPU/12GiB and the unchanged80GiB cap/12GiB reserve.
Actual terminal exit0/noOOM and raw source/runtime were reviewed. At15M leaves,
peak RSS is6,047.98MiB for blocks /6,355.39MiB for patterns; corresponding grid
setup is78.29/73.16s and four-worker planning146.24/75.27s. Exact ownership,
64 curl/12 field samples per layout pair and first-step energy drift<=6.22e-15
pass. These are single observations of homogeneous host preparation, not native
PCB/GPU/RF/long-time or full-array admission. Public default remains1M.
Evidence: `scaling-budget-r50/review-r01.json`.

R51 is a separately frozen, now CPU/CUDA-qualified optional device field-energy increment
over R49. It includes off-diagonal Galerkin electric mass, unique owned states,
scaled compensated reductions and complete-epoch guards; native processing is
unchanged. The original test configured a resistor after unused soft sources
and correctly rejected. A preserved test-only correction passes185 energy/
loss/PML checks on the unchanged core, maximum relative error1.15141e-8. The
fresh corrected build and CPU68/68 pass, with16 bit-identical actual-R49 arrays
and779 source hashes rechecked. All52 physical gates finish exit0/noOOM after exact R47/R49 release and fresh
admission checks;779 source/413 archived files were independently rehashed.
Exact numerical/API/test source is promoted. The original external collector
confused52 commands with48 inherited native-box rows; corrected row cardinality
only, no numerical restart or threshold change. Original failure/audit retained.
Metal/native use remains open.
Evidence: `device-energy-r51/`; R01 failure is not overwritten or solver-fixed.


R47 R04 now completes all eight matched-step/runtime cases exit0/noOOM and passes
its original scoped fixed-XY RF gates. Independently rehashed228 raw files and
17-point own-timestamp DFT across12 cases reproduce spectra within6.02e-16;
primary temporal/band analysis covers501 frequencies. At Z3, local/ordinary
complexGamma difference is0.006073 PEC /0.000939 copper, S21 phase0.01791/
0.13811deg and magnitude0.03493/0.00334dB. Local Z2-to-Z3 differences are0.004142/
0.004038 complexGamma, under0.01. Runtime replay, late prefixes and5-to4fs
sensitivity also pass. R37/R41 failures remain; fixed fine XY/two finer Z
levels do not prove asymptotic3D, distributed RF or antenna acceptance. Original
5fs Z3 copper setup rejection remains. No isolated timing claim after R49
control overlap. Evidence: `finite-width-z3-r47/review-r05.json`.

R52 separately freezes an exact bounded host material palette over promoted R49,
excluding unpromoted R51. Repeated epsilon/mu retain exact bits; homogeneous
storage needs no leaf IDs, up to4096 pairs use32-bit IDs, then dense fallback
bounds dictionary growth. The dictionary is discarded before field construction;
fallback peak includes temporary IDs. Focused lifecycle, nextafter, copy/move,
dense and1472-leaf sample controls pass; full CPU68/68 passes with16 exact frozen
parent arrays. Matched3.93M host characterization also passes:62,914,424 fewer
accounted bytes, peak RSS1711.08->1650.57MiB patterns /1632.61->1572.60MiB blocks
(3.54/3.68%). Independent source/raw/terminal review is retained. These single
homogeneous measurements are not variable-medium/native PCB/GPU capacity. All48 physical CUDA gates now pass exit0/noOOM on4090/3080 and both orders,
with775 source and409 archive files independently rehashed. Exact R51 physical
terminal and fresh resource/identity guards were checked before admission.
This separately qualifies the R49-based candidate. R54 below completes its
composition with current native observations before exact promotion. Metal
and native PCB/antenna scaling remain open.
No production/default guard change. Evidence: `material-palette-r52/`.


R53 separately stages native startup U/I registration, selected batch readback
and cached scalar device energy over exact R51-r02 (R52 excluded). The native
application seals after processing setup; manual C++ engines remain eager until
explicit seal. Generic E/H dumps and explicit ObservationGrid refresh full fields
once per epoch, guarded for concurrent readers; raw setup Grid references are
not live mirrors after opting into sealing. Initial double host fields still
exist; prepared E/H/DFT gathering is required follow-up. A new12-case control
covers single/two workers, all cuts, PML/bulk loss, U/I/SI functionals, scalar
energy/cache accounting, exact full fields and concurrent fallback. Frozen build and12 focused reference cases pass: functional discrepancy0,
maximum energy error7.60313e-9, no full reads across four processing batches,
and exactly one refreshed full fallback per observed epoch under concurrent
readers. Full CPU69/69 and16 actual-parent operator/field files now pass;780 source files were rehashed. Physical56 now passes exit0/noOOM;780 source and418 archived files were independently rehashed. Exact numerical/API/test source is promoted. Metal/native RF remain pending. The
host-observation byte counter excludes registry/reader/caller memory.
Evidence: `native-sparse-observations-r53/`. No production runtime is changed.


R54 freezes the exact R53 native sparse observation/energy changes with the
separately CPU/CUDA/RAM-qualified R52 material palette. CPU70/70,16 exact
actual-R53 files and56 physical CUDA gates pass exit0/noOOM. Independently
rehashed782 source and419 archive files before exact numerical/API/test
promotion. R52 memory numbers remain separately scoped homogeneous evidence;
no combined native PCB/antenna capacity or Metal qualification is inferred.
No production runtime changes. Evidence: `native-material-composition-r54/`.

R55 exact implicit initial +0 E/H storage passes CPU71/71,16 bit-identical
actual-R54 operator/field files and60 physical CUDA gates, exit0/noOOM. All784
source and424 archive files were independently rehashed before exact numerical/
API/test promotion. Read-only metadata, domain setup, selected observations and
energy retain zero value-array allocation; CPU writes/full imports materialize
fields. Four matched3.93M homogeneous host controls save281,743,856 accounted
initial bytes (268.69MiB). Plan peak RSS falls19.08% patterns/20.25% blocks;
post-write storage is comparable. These single host measurements are not native
PCB/GPU/antenna admission or FDTD speedup. Metal and representative RF remain
open. Evidence: `implicit-zero-state-r55/independent-qualification-r01.json`.


R56 separately stages native E/H dump gathers over frozen R55. It reuses the
existing raw-region FDTD interpolator, with SI-to-raw float conversion before
interpolation, bounded table preparation and once-per-epoch caches. Application
sealing is after processing initialization; explicit manual sealing/eager full
fallback remain. The new12-case control compares all modes, both fields,
graded/unequal child coordinates, source/PML/loss, drawing transforms, four
concurrent output partitions and complex DFT sums against direct full fields
with a predeclared2e-6 gate. SI stencils can reassociate double operations; no
bit-identity claim. CPU72 and physical64 are queued behind actual predecessor
acceptance/release. Device DFT, generic E/H probes, Metal, collective pause,
RF/NF2FF and antenna qualification remain open. No production runtime changes.
Evidence: `native-field-gather-r56/preparation-r01.json`.


R56 R01's isolated build failed before numerical execution because a borrowed
const interface called an existing nonconst interpolation getter. Preserve its
exact failed terminal/source/diagnostics. R02 adds one const getter overload,
retaining the existing nonconst signature and every adapter/test/physics/gate
unchanged. Its isolated build/CPU72 pipeline is admitted after accepted R55;
physical64 waits for CPU acceptance and the exact accepted R55 pair release.
This remains an unqualified source, not a promoted solver fix or antenna build.
Evidence: `native-field-gather-r56/build-failure-review-r01.json`,
`native-field-gather-r56-r02/preparation-r01.json`.


R56 R02's twelve focused reference cases pass with zero measured field/complex
DFT discrepancy, but a separate minimal corner reproducer finds SIGSEGV when
the raw region has no support. This is a new adapter defect, not an upstream
FDTD defect. Only its waiting physical-admission watcher is held; CPU evidence
continues unchanged. R03 keeps a budgeted one-float region sentinel and adds
empty upper-E-cell/lower-H-node cases in each of twelve source/PML/loss/single/
domain controls. Its isolated build passes; full CPU72 is running and physical64
remains gated on acceptance and exact qualified R55 release. No R56 source is
promoted. Evidence: `native-field-gather-r56-r02/zero-region-review-r01.json`,
`native-field-gather-r56-r03/preparation-r01.json`.

The exact independently compiled R55 memory driver and reproduction steps are
also available in this fork as `benchmark/subgrid_initial_memory.{cpp,md}`;
this is a host-memory characterization, not an FDTD speed benchmark.


R56 R03's twelve focused CPU/reference cases now pass, including24 empty E/H
boundary controls. Maximum measured field and four-batch complex-DFT discrepancy
is0; no full read occurs before the explicit fallback. All787 source files
and the executing test/library/log identities were rechecked. This focused
receipt does not replace full CPU72, physical64, Metal or RF/antenna acceptance.
Evidence: `native-field-gather-r56-r03/focused-reference-r03.json`.


R56 R03 full CPU72/72 and16 exact actual-R55 parent arrays now pass, exit0/noOOM.
Independent review rehashes787 source files and the compact raw archive;12
focused cases include24 empty regions and maximum measured field/DFT error0.
Physical64 is running after the exact accepted R55 release; no R56 source is
promoted yet. Evidence: `native-field-gather-r56-r03/independent-cpu-review-r01.json`.

R57 separately implements owned normalized-state DFT through the shared point
operations and device/coordinator/worker protocols. Inline weights avoid blocking
uploads; selected states have unique owners, bounded aggregate memory and
complete-epoch/duplicate-sample checks. The original99,999-sample cancellation
control exposes correction-accumulator drift (relative2.70533e-4); preserve its
exact negative/source/build. R02 separately renormalizes the pair without gate
relaxation. Its isolated build and focused reference control pass (six domain/
port cases, max complex error1.35356e-14, exact unchanged single-device fields).
Full CPU73 and physical68 remain pending. Native FieldDFT/DFT-only registry
integration, Metal, collective pause, matched RF and antenna qualification are
still required. No R57 source is promoted or production runtime changed.
Evidence: `device-dft-r57/negative-review-r01.json`,
`device-dft-r57-r02/focused-reference-r01.json`.


R56 R03 now passes physical64/64, exit0/noOOM, after the exact R55 release.
Independent review rehashes787 source and429 archive files. All48 physical
native gather cases include96 empty-boundary controls; maximum measured raw
field/complex-DFT error is0. Combined with CPU72/72 and16 exact parent arrays,
the exact nine numerical/API/test/build files are promoted. All583 source and
build files from the frozen snapshot match the canonical branch. The older
R01 compile failure/R02 corner crash remain preserved. A collector summary
count was corrected56->64 to match its unchanged64-command assertion; no
numerical rerun or gate change. Metal, device DFT, collective pause, native RF
and representative antenna acceptance remain open. Production is unchanged.
Evidence: `native-field-gather-r56-r03/independent-qualification-r01.json`,
`native-field-gather-r56-r03/promotion-r01.json`.


R58 separately prepares native device FieldDFT over frozen R57 R02, retaining
the existing ProcessFieldsFD API/weights/scheduling. Owning recipes reconstruct
complex raw regions from owned-state sums; dedicated device-only recipe IDs
are filtered from batch host selection, with genuine U/I/TD support retained.
Explicit later TD access uses live full fallback. Whole-budget rejection is
pre-upload; identity/shape/frequency/epoch checks and failed-state propagation
are explicit. Twelve new native controls cover all modes, E/H, child transforms,
source/PML/loss, single/domains, tiny timestep-scaled weights, empty support and
repeated snapshots. Source794 files is frozen; build/CPU74 waits for actual
R57 CPU73 acceptance, physical72 for exact accepted R57 physical68 release.
Implementation is staged, not qualified or promoted. Recent Metal replay is
also queued with the unchanged1GiB compiler/1.5GiB test,12.5GiB reserve and
15.5GiB stable60s admission guards. No production/antenna runtime changes.
Evidence: `native-device-dft-r58/preparation-r01.json`,
`native-device-dft-r58/DESIGN.md`.


R57 R02 now passes full CPU73/73 and16 exact actual-R56 operator/field arrays,
exit0/noOOM. Independent review rehashes792 source files and the compact raw
archive. Six focused DFT domain/port cases have maximum complex error
1.35356e-14, with exact unchanged single-device fields; long cancellation and
empty/duplicate/order/budget/epoch controls pass. Physical68 is running after
the exact accepted R56 release. No R57 source is promoted yet. R58 native
integration build is now running behind that accepted CPU parent, not yet
qualified; recent Metal still waits for unchanged headroom admission.
Evidence: `device-dft-r57-r02/independent-cpu-review-r01.json`.


R57 R02 physical68/68 now passes exit0/noOOM after exact accepted R56 release.
Independent review rehashes792 source and434 archive files. Twenty-four physical
DFT domain/port cases have maximum complex error1.35356e-14; fields remain
unchanged and accumulation performs no readback. Combined with CPU73/73 and
16 exact parent arrays, thirteen exact code/API/test/build files are promoted;
all587 frozen source/build files match. Native FieldDFT remains separate R58: its
build/twelve focused reference cases pass (max whole-complex error2.85387e-7),
with eleven cases having no host E/H support and no full reads before explicit
fallback. Full CPU74 is running; physical72 remains gated. A matched public
dipole RF/NF2FF host/device/single/pair control is separately staged as R59 R02,
waiting for actual R58 acceptance. Recent Metal still waits for unchanged
headroom; no production/antenna runtime changes. Evidence: `device-dft-r57-r02/
independent-qualification-r01.json`, `device-dft-r57-r02/promotion-r01.json`,
`native-device-dft-r58/focused-reference-r01.json`.


R58 now passes full GPU-free CPU74/74 and16 exact actual-R57 operator/field
arrays, exit0/noOOM. Independent review rehashes794 source files. Twelve native
complex-array cases have maximum relative error2.85387e-7; eleven retain no
host E/H selection. Actual physical72 is separately admitted and observed by
container b1579d3a2eba27269fbcbef56957e3bf5807e51fd369b3ef6edd478b6e0e9314;
no R58 source is promoted before physical acceptance. The matched public-dipole
host/device/single/pair NF2FF harness is R59 R03. R59 R02 was still unlaunched
when static review corrected its final percent-comment port-file reader;
its driver/public fixture and gates are unchanged, with the review preserved.
Evidence: `native-device-dft-r58/independent-cpu-review-r01.json`,
`native-device-dft-r58/cpu-qualification-r01.json`,
`native-dft-radiation-r59-r02/static-harness-review-r01.json`,
`native-dft-radiation-r59-r03/preparation-r01.json`.

R60 stages an explicit completed-epoch owner-thread fence over exact R58.
A worker selects/waits for its backend before acknowledging that exact epoch;
the coordinator returns only after every logical worker succeeds. Wrong or
missing acknowledgements/backend failures prevent continuation. No field or
DFT download, migration, ordinary-step barrier or competing control-file
protocol is introduced. Two focused targets include actual backend-wait call
verification before return, half-step/epoch/failure negatives, unchanged resumed
fields/DFT and native PML/bulk-loss cases.799 source files are frozen; the
isolated GPU-free build is admitted as
edd09d25cb9de520434602a0376ee91a03e934c89784c7df75b90b4f246a431f.
Expected fullCPU76, physical80 and external pause/snapshot integration remain
pending. Physical admission must also follow the exact accepted R59 job release.
This is staged implementation, not numerical or antenna qualification.
Evidence: `collective-fence-r60/preparation-r01.json`,
`collective-fence-r60/DESIGN.md`.


R58 physical72/72 now passes exit0/noOOM after the exact accepted R57 release.
Together with CPU74/74 and16 exact parent arrays, independent review rehashes
794 source/439 archived files. Forty-eight physical native complex-array cases
have maximum relative error2.85387e-7, with eleven of twelve cases per lane
having no host E/H selection. Eight exact code/test files are promoted after
checking all source/build files against the frozen accepted tree. Native RF
remains a separate R59 harness: R03 failed linking its external driver before
FDTD; R04 links existing CSXCAD explicitly with unchanged driver/fixture/gates.
R60 isolated build and six direct/twelve native reference fence cases pass;
full CPU76 and physical80 remain pending. Recent Metal admission expired
without starting a compiler because unchanged headroom guards were not met;
no Metal acceptance is inferred. No production/antenna runtime changes.
Evidence: `native-device-dft-r58/independent-qualification-r01.json`,
`native-device-dft-r58/promotion-r01.json`,
`native-dft-radiation-r59-r03/failure-review-r01.json`,
`collective-fence-r60/focused-reference-r01.json`.


R59 R06 matched native RF/NF2FF now passes all four actual4000-step public-dipole
runs and independent raw review, wrapperexit0/noOOM. Same byte-exact native XML,
qualified R58 library, one4090 and mixed4090/3080 pair, each host/device DFT.
Port U/I/Z/Gamma are exact; maximum relative far field1.80653e-7, near field
2.15114e-7 and radiation quantities8.067e-8. Single/pair exports are exact.
Independent review rehashes172 raw artifacts and retains both enclosing surfaces
and public power/pattern/contact checks. Earlier harness failures are preserved:
external link, constructor-before-parser order, and an overstrict random-display-
colour identity gate. R05 independently passes unchanged numerical gates; R06
also restores the original strict byte gate with deterministic display XML.
None was an established solver defect. Physical R60 remains gated by fullCPU76;
its waiting owner now requires actual accepted R58 and R59 R06 identities.
Evidence: `native-dft-radiation-r59-r06/independent-qualification-r01.json`,
`native-dft-radiation-r59-r06/README.md`,
`native-dft-radiation-r59-r05/independent-qualification-r02.json`.

R61 is separately staged, unbuilt and unfrozen over R60: it reuses the existing
runtime/startup/snapshot file protocol and Python/Matlab request helpers. Native
workers are fenced after processing drains and before paused/snapshot publication;
acknowledged logical-worker count/epoch are recorded. Nine pure Python request
helper tests pass, but establish no native runtime qualification. Matched actual
pause/resume/snapshot controls, full build/replay and physical/Metal gates remain
required before any promotion. Production/frozen runtimes are unchanged.
Evidence: `external-pause-r61/reuse-r01.json`, `external-pause-r61/DESIGN.md`,
`external-pause-r61/request-helper-tests-r01.json`.


R60 full GPU-free CPU76/76 now passes exit0/noOOM,42 commands and16 exact actual-
R58 operator/field arrays. Independent review rehashes799 source files; focused
six direct/twelve native fence controls include actual backend-completion calls,
resumed fields/DFT and malformed/missing/failed acknowledgement guards.
Initial physical admission held before docker run because GPU5 had57MiB free;
no numerical R60 physical job ran there. A named resource-only R02 harness keeps
all binaries/gates and uses measured GPU2RTX4090 plus free GPU6RTX3080, with exact
accepted R58/R59 release and fresh guards. Physical80 is now running as
818d3be7cad5d00fd9b712bbb8675dc4b9edbc40d63cc8db254070362b985470.
No physical or Metal acceptance is inferred from CPU76.
Evidence: `collective-fence-r60/independent-cpu-review-r01.json`,
`collective-fence-r60/alternate-pair-r02.json`,
`collective-fence-r60/alternate-pair-collector-review-r01.json`.

R61 source811 files and its byte-exact public dipole/application harness are now
frozen separately. The admitted GPU-free build waits for and follows actual R60
CPU76 acceptance; build ID
d72ab3674604205cccb7418e139712b0ea70a0ea5fbb31e9b857950f84628248.
The new application controls cover native basic/reference/split-reference and
ordinary basic/reference, startup/manual pause, repeated same-epoch snapshots,
resumed final-array identity and forged controller IDs. Implementation and nine
pure request-helper tests do not qualify the pending native application gates.
Source, device controls, recent Metal and representative antenna remain separate.
Evidence: `external-pause-r61/preparation-r01.json`,
`external-pause-r61/application_control.py`.


R60 physical80/80 now passes exit0/noOOM on measured GPU2RTX4090/GPU6RTX3080,
including both pair orders. Combined with CPU76/76 and16 exact actual-parent
arrays, independent review rehashes799 source/449 archived files.24 physical
direct fence cases verify completion waits before returning;48 native cases
preserve resumed fields/DFT exactly with no fence downloads. Exact qualified
core/API/test files are promoted after complete source/build comparison.
Recent Metal, external application control R61 and antenna qualification remain
open. The resource-only pair derivative does not change any source or gate.
Evidence: `collective-fence-r60/independent-qualification-r01.json`,
`collective-fence-r60/promotion-r01.json`.


R61 isolated build passes. The first full CPU application attempt preserves a
harness-only same-epoch assertion failure: queued manual and startup pauses
correctly acknowledge epoch132. A named corrected harness permits equality and
also checks post-pulse -10/-20 checkpoints. A separate small protocol fixture
correctly rejected missing PML one-cell separation; R04 adds the required guard,
without changing source/ports/metals/timestep/4000step count. Five3180-cell CPU/
reference/ordinary protocol modes now pass with16 exact final arrays each.
This small control does not replace the full RF fixture.

Independent physical qualification uses the unchanged full152204-cell public
dipole on one4090 and actual4090/3080 pair after accepted R60 physical80 and
focused CPU protocols. Both lanes pass exit0/noOOM:160 final raw arrays each
match uninterrupted execution AND accepted R59 arrays exactly. Repeated snapshots
retain sample metadata; forged IDs reject.618 raw artifacts are rehashed.
Full R61 CPU76/16 actual-parent replay remains running; code is not promoted.
Timeout/stop/write-failure application checks, recent Metal, native PCB scaling
and representative antenna acceptance remain open. No production changes.
Evidence: `external-pause-r61/focused-protocols-r03.json`,
`external-pause-r61/physical-application-independent-r01.json`,
`external-pause-r61/application-failure-review-r01.json`,
`external-pause-r61/protocol-fixture-r04.json`.


R61 full GPU-free replay now passes exit0/noOOM:76/76 CTest,50 scripted
commands,16 exact actual-R60 arrays,811 frozen source files and531 archived
files independently rehashed. Five protocol modes preserve startup/manual/
post-pulse checkpoints and repeated DFT metadata. Full152204-cell public dipole
controls on one4090 and an actual4090/3080 pair preserve160 final raw arrays
per lane, including accepted R59 outputs. No new RF extraction is inferred.

Additional failure controls pass10 child cases each on native CPU, two logical
reference workers and the actual physical pair. Runtime timeout resumes exactly;
graceful stop acknowledges epoch132 without convergence; startup reject/timeout
fail; stale control files reject before operator creation. Reused snapshot IDs
and blocked staging do not overwrite published data, do not publish partial
results, and allow a fresh same-epoch request followed by exact resumed arrays.
198 CPU/reference and342 physical raw artifacts are independently rehashed.
The initial physical fault launcher hit an existing name before any test
started; the named R02 admission changes only its fresh directory/container name.
The background CPU owner fetched first; a second observer's exclusive-output
rejection did not rerun or replace numerical evidence.

Exact R61 core/API/helper/test bytes are promoted after whole-code comparison.
Public reproduction drivers/fixtures and semantics are in
[CARTESIAN_SUBGRID_RUN_CONTROL.md](CARTESIAN_SUBGRID_RUN_CONTROL.md). This extends
partial D04/U04 only. Recent Metal, application-level backend-fault injection,
opt-in overhead, native PCB scaling and representative antenna acceptance
remain open; no production/frozen runtime changes. R60 core failure guards do
not substitute for application fault injection.
Evidence: `external-pause-r61/cpu-independent-r03.json`,
`external-pause-r61/physical-application-independent-r01.json`,
`external-pause-r61/fault-controls-independent-r01.json`,
`external-pause-r61/physical-fault-independent-r02.json`,
`external-pause-r61/promotion-r01.json`.


R62 now executes the unchanged accepted R47 Z3/t4 native smooth-copper U-line
(1,544,144 leaves,4fs,500k steps) on the exact R61 installed runtime: one
RTX4090, then actual4090/3080 domains in both worker orders. Fresh4CPU/6GiB,
80GiB project/12GiB reserve and physical VRAM checks pass; actual single-device
FDTD progress is verified. The original observer quoting issue affected only
status reads; its named replacement retains the same running container ID.
Predeclared runtime/domain RF gates are complexGamma1e-5,S21 phase0.01deg and
magnitude0.001dB over24.000–24.250GHz. Running is not numerical acceptance,
and this small-network replay does not qualify a representative antenna,
large-PCB resources or general three-dimensional curved convergence.
Evidence: `pcb-pair-rf-r62/preparation-r01.json`, `pcb-pair-rf-r62/DESIGN.md`,
`pcb-pair-rf-r62/admission-r01.json`, `pcb-pair-rf-r62/owner-r02.json`.


R63 is separately frozen/unpromoted over exact R61: an owning64-bit integer
edge-address buffer replaces five host integers per edge, preserving canonical
IDs and exact lattice geometry. The current16384-interval/refinement16 guards
fit19 bits per coordinate; malformed values reject. Source/build and the
3241-value/11-negative focused test pass. Full CPU77/16 parent-array replay
and three interleaved3.93M-leaf host comparisons per layout/prefix are running.
The expected address-buffer reduction is not a measured total-memory win or
antenna admission. No numerical/production source is promoted; device and
recent Metal acceptance remain separate. A fresh local headroom sample is
below the unchanged15.5GiB admission threshold; no compiler/test was launched.
R62 continues on its unchanged R61 runtime.
Evidence: `packed-edge-addresses-r63/DESIGN.md`,
`packed-edge-addresses-r63/preparation-r01.json`,
`packed-edge-addresses-r63/cpu-admission-r01.json`,
`packed-edge-addresses-r63/memory-admission-r01.json`,
`packed-edge-addresses-r63/metal-headroom-r01.json`.


R63's independent host-memory comparison now accepts all12 interleaved
3,932,160-leaf cases:192MiB less accounted storage before/after a CPU write,
exact geometry/curl/fields/timestep/ownership and unchanged default guards.
Median grid/domain preparation peak RSS falls15.85/11.26% (patterns) and
17.67/12.57% (blocks); including the CPU write gives5.66/8.42% lower peak.
Whole preparation/control medians change-0.08/+0.41% on a shared host; no
speedup is claimed. Full CPU77 remains running and physical82 is staged only
after independently reviewed CPU/memory/R62 release. Code stays unpromoted;
this homogeneous host result does not qualify native PCB/GPU/antenna capacity.
Evidence: `packed-edge-addresses-r63/memory-independent-r01.json`,
`packed-edge-addresses-r63/physical-preparation-r02.json`,
`packed-edge-addresses-r63/physical-enumeration-r01.json`.


R63 full GPU-free qualification now passes exit0/noOOM:77/77 CTest,51 scripted
controls and16 exact actual-R61 operator/field arrays. Independent review
rehashes819 source and532 archive files. With the separately accepted memory
comparison this qualifies the scoped CPU/reference storage change; physical82
still waits for accepted R62 RF release, and recent Metal remains unadmitted.
The candidate code is not promoted and no antenna/scaling qualification follows.
Evidence: `packed-edge-addresses-r63/cpu-independent-r01.json`.


R62's unchanged native1,544,144-leaf smooth-copper U-line replay is now accepted
on one4090 and the4090/3080 pair in both orders. Independent raw own-timestamp
DFT over501 band points matches accepted R47 complexGamma within3.822e-16;
fresh single/pair extractions are identical. Original1e-5Gamma/0.01deg/0.001dB
runtime gates and1.5ns/2ns late-prefix gates pass;81 artifacts are rehashed.
Single/pair/reverse wall times543.08/892.80/964.77s show no small-case speedup.
This is fixed-XY/Z3 network evidence, not representative antenna convergence.

After actual R62 terminal release and accepted R63 CPU77/memory receipts,
R63 physical82 was admitted on GPU2RTX4090/GPU6RTX3080 as
`e832c2b381fccf88ae8645dff1dd9db106c03f134cf8eee2d64156c9519fb508`.
Source/storage code remains unpromoted until this gate and independent review
pass. Recent Metal and representative antenna qualification stay open.
Evidence: `pcb-pair-rf-r62/independent-review-r01.json`,
`packed-edge-addresses-r63/cuda-owner-r01.json`.


R63 now passes all82 physical controls on actual GPU2RTX4090/GPU6RTX3080,
including both orders, native ports/materials/PML/DFT and full public application
pause/snapshot/resume.160 raw final arrays per application lane match both
uninterrupted execution and actual R61 outputs exactly. Independent review
rehashes819 source/1062 archive files. The frozen collector's missing CLI-map
key and inherited CPU archive/test-count assumptions caused an observer-only
failure after successful numerical completion; the named R04 adapter checks the
actual application start/end CLI hashes against CPU and all77 tests, retaining
all other acceptance gates. No numerical rerun or tolerance relaxation occurred.

Exact packed-address core/header/test bytes are promoted only after CPU77,
16 exact actual-parent arrays, physical82 and the accepted12-case interleaved
host-memory comparison.192MiB accounted savings at3.93M leaves and measured
preparation RSS reductions are scoped, not antenna capacity or speed claims.
Matching C++ callers must be rebuilt. Recent Metal, native large-PCB resources,
further curved spatial RF convergence and representative antenna remain open.
Evidence: `packed-edge-addresses-r63/independent-qualification-r04.json`,
`packed-edge-addresses-r63/collector-review-r04.json`,
`packed-edge-addresses-r63/promotion-r01.json`.


R64 prepares a further XY sensitivity level of the exact accepted curved PCB
U-line, retaining all old drawing knots, properties/ports, physical copper,
waveform, boundaries and Z3. Native leaves rise1,544,144 to5,933,200; matched
ordinary Yee has11,565,280 cells. Six PEC/copper/common2fs controls are frozen.
Only GPU-free setup is admitted under4CPU/8GiB/reserve12GiB; actual native CFL,
resources and completed preparation must pass before any RF launch. This is
not a contour/antenna change or accepted spatial refinement. Runtime is the
separately qualified R63 candidate; no production binary is modified.
Evidence: `curved-xy-refinement-r64/DESIGN.md`,
`curved-xy-refinement-r64/preparation-r01.json`,
`curved-xy-refinement-r64/setup-admission-r01.json`.


R64 common1fs native preparation passes four PEC/copper old/refined cases.
The two initial ordinary preparations used R63's regular path, which lacks
LossyMetal; their mere completed setup did not qualify a lossy reference.
The R02 RF diagnostic was explicitly stopped during its first native case
(exit130), preserving partial data and claiming no RF acceptance. R03 retains
all XML/grid/material/port bytes and restores the actual pinned PR242 ordinary
runtime. Both fresh ordinary setups now pass, with119 actual conductor fits
for copper; linked runtime and conductor-provider witnesses are mandatory.
Fresh RF runs on GPU2RTX4090 as
`d66abb78a0e5bcc7485c76ec00899fc05d20b6cc4833708dfa2daca22c2f2a8f`.
No RF or antenna acceptance follows from these setup/process results.
Evidence: `curved-xy-refinement-r64/reference-audit-stop-r01.json`,
`curved-xy-refinement-r64/reference-correction-r03.json`,
`curved-xy-refinement-r64/reference-setup-review-r01.json`,
`curved-xy-refinement-r64/rf-admission-r03.json`.

R65 independently stages an orchestration-only capability guard over exact
qualified R63: reject optional LossyMetal on this branch's unsupported regular
path before operator/backend creation; supported native and legacy PEC remain.
No update equation, kernel or production runtime changes. Focused eight negative
and four supported operator cases, actual R63 failure reproduction, CPU78 and
16 unchanged parent arrays are required before promotion. The isolated
4CPU/6GiB GPU-free build is running; no implementation/qualification completion
is inferred from source inspection. Recent Metal/physical acceptance remains
separate. Evidence: `ordinary-metal-guard-r65/DESIGN.md`,
`ordinary-metal-guard-r65/preparation-r01.json`,
`ordinary-metal-guard-r65/cuda-snapshot-r01.json`.
