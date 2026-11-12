# Cartesian subgridding — implementation and qualification checklist

Last reviewed: **12 November 2026**.

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
- Pinned master baseline: `1d710c210920908a143b573f2babca05455f5479`.
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
- Current capability is a **serial CPU numerical reference**, with homogeneous
  lossless background, cell-aligned isotropic lossless material boxes and PEC
  outer walls. It is not an antenna-ready solver.
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
one million base cells, one million total leaves, and 1024 coupled degrees of
freedom per dense interface factor. These are implementation guards, not
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
- [ ] **G03** Implement internal PEC objects, including finite physical metal
  height, sheets where supported, vias, gaps and finite-width contacts.
- [ ] **G04** Qualify dielectric and conductor crossings of coarse/fine
  boundaries. Any interim placement restriction must fail explicitly.
- [ ] **G05** Support the anisotropic material coefficients available in the
  pinned regular Cartesian solver, or document a remaining parity exclusion.
- [ ] **G06** Add layered-medium, dielectric-interface, cavity and PEC-object
  controls against analytic references and the regular fine-grid solver.

Acceptance: independently inspect material occupancy and contacts, verify
geometry priorities, and demonstrate convergence and bounded evolution.
The current outer PEC-wall tests do not complete internal PEC mapping.
G01/G04/G06 are **partial**: native volume-box priorities, units and leaf
occupancy, dielectric crossings of refinement faces, analytic layered-cavity
convergence and a graded dielectric/native-Yee control pass. General primitives,
internal PEC contacts and broader scattering qualification remain open.

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
- [ ] **P02** Implement lumped-port source/resistor coupling on local fields,
  plus voltage-path and current-contour integration across local cells.
- [ ] **P03** Support the remaining circuit/lumped RLC elements and qualify
  energy injection, stored energy and dissipation.
- [ ] **P04** Implement microstrip/coaxial/waveguide port workflows as
  applicable, including modal processing and reference-plane conventions.
- [ ] **P05** Validate single- and multiport power balance, complex impedance,
  reflection/transmission and passive terminations against regular fine grids.

Acceptance: qualify a short 50-ohm transmission line before an antenna.
Sampling reconstructed E/H on base observation coordinates is not, by itself,
a qualified physical port.

## 4. Boundaries and radiation postprocessing

- [ ] **B01** Implement and test PMC and supported symmetry boundaries.
- [ ] **B02** Integrate UPML and Mur boundaries. An initial requirement to keep
  refinement outside absorbing layers must be explicit and checked.
- [ ] **B03** Measure absorbing-boundary reflection over relevant frequencies,
  angles and polarizations, including nonuniform meshes.
- [ ] **B04** Qualify DFT accumulation, field timing, spatial interpolation and
  integration weights on NF2FF surfaces sampling the local grid.
- [ ] **B05** Compare small-antenna input response, near fields, radiation cuts,
  power and efficiency with regular fine-grid and analytic controls.

Acceptance: a radiating PEC antenna must work before treating the prototype as
an antenna solver. Correct E/H file creation alone does not qualify NF2FF.

## 5. Scalable representation and CPU execution

- [ ] **C01** Use structured storage/update kernels for regular block interiors;
  reserve compact constraint metadata and solves for the interfaces.
- [ ] **C02** Replace or scale the current dense interface-factor approach;
  measure setup cost, factor memory and per-step cost before raising guards.
- [ ] **C03** Reduce reference-only storage overhead, with regression coverage
  for any change to magnetic-face histories or numerical quadrature.
- [ ] **C04** Integrate the regular CPU execution options, threading and SIMD
  where applicable; document remaining engine differences explicitly.
- [ ] **C05** Benchmark peak setup/steady-state memory, setup time and time to
  the same simulated physical duration and accuracy, against regular grids.

Acceptance: do not infer production memory or speed gains from leaf counts.
Record minimum timestep, interface overhead, field precision and build identity.

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

## Milestone decisions

| Milestone | Completion condition | Status |
|---|---|---|
| Reference prototype | R01–R08 demonstrated on small CPU controls | Achieved, limited scope |
| Useful PEC antenna on CPU | Required geometry/mesh/port/PML/NF2FF items, with small-antenna evidence | Open |
| Useful PEC antenna on one GPU | Previous milestone plus scalable representation and a qualified device backend | Open |
| Representative lossy antenna | Required conductor support and A01–A07 | Open |
| Pinned-master Cartesian parity | Full applicable feature/engine matrix and regression suites | Open |
| Production multi-GPU subgridding | Relevant single-device physics qualified, then D01–D05 | Open |

Suggested first implementation sequence: isotropic materials and internal PEC,
50-ohm line/port, absorbing boundaries and small-antenna NF2FF, scalable
interiors and one GPU, then lossy antenna validation and remaining parity work.
Nonuniform base geometry is implemented and qualified on small homogeneous
controls; application compatibility still requires the material, port and
boundary gates above. These are dependency gates, not a
commitment that every master feature must precede the first useful PEC model.

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
