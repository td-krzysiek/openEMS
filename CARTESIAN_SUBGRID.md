# Experimental compact Cartesian XYZ subgridding

`CartesianSubgridBox` replaces selected base cells with `rx * ry * rz`
children. `CartesianSubgridRatio` independently selects the refinement on
X, Y and Z: each factor is 1, 2, 4, 8 or 16. The default `[2,2,2]` retains
the original eight-child behaviour; `[1,2,8]` makes 16 children while leaving
X unchanged, and `[2,4,8]` makes 64 children. Spacing on each axis is divided
by its selected factor inside each selected box. Several separate boxes
can have different axis ratios; outside them the coarse cells remain. This is a separate Cartesian numerical core, not
angular interpolation or a slab copied from cylindrical multi-grid.

It is a **serial CPU reference implementation**, suitable for numerical
review and small validation models. It is not yet a replacement for the
production antenna engine.

[Implementation and qualification checklist](CARTESIAN_SUBGRID_CHECKLIST.md)
tracks the completed controls and remaining feature/engine milestones in PR #29.

## Geometry, storage and execution

- One refinement level and 1 through 64 independent rectangular boxes, each
  with its own configurable axis ratios, strictly inside the domain. Every
  face coincides with a base-grid line. At least one base cell remains between
  each box and the outer boundary, and a full base-cell gap separates every
  pair of boxes on at least one axis. Face/edge/corner touching, overlap and
  nesting are rejected; touching blocks cannot yet be used as a transition.
- A Cartesian base mesh with finite, strictly increasing lines on each axis.
  Nonuniform intervals and a nonzero origin are supported. Each selected base
  interval is divided by its region's ratio; physical widths, faces, volumes,
  lookup and field reconstruction use those actual coordinates. Nested blocks
  and moving refinement remain unsupported.
- Region bounds must match supplied base lines. Only floating-point coordinate
  roundoff is tolerated (32 machine epsilons of coordinate scale, additionally
  limited to 1e-8 of adjacent spacing); geometric snapping is not performed.
  Unrepresentable child boundaries or sample centres fail at setup.
- Coarse cells covered by the box are removed. Only leaf cells and independent
  electric edges are stored. There is **no global fine parent E/H lattice**.
- One common timestep bounded by the coupled operator, not local subcycling.
  The smallest scale still limits the timestep.
- The reference uses double fields, sparse curl rows and interface mass
  factors. This favours inspectability, not production speed. Regular
  interiors need native structured kernels before a fair performance claim.

For base count N, M_i selected base cells and `R_i=rx_i*ry_i*rz_i` in region i,
the leaf count is `N + sum_i (R_i-1) M_i`. A fully fine control uses the maximum
requested factor on each axis: `N * max(rx_i) * max(ry_i) * max(rz_i)`.
The maximum factors define a common integer coordinate lattice; they do not
refine unselected cells or override any region's ratio. No full lattice of
field states is allocated. `LatticeRatio()` (and compatibility alias `Ratio()`)
reports this coordinate scale, while `Regions()` reports each local ratio. In the test with 8x8x8 base cells and a
4x4x4 selected block, the counts are 960 versus 4096. Accounted persistent
storage is about 1.10 versus 4.40 MB **for this same reference core**. These
numbers are not process RSS, peak setup memory or a comparison with the
optimized native Yee/GPU storage. Sparse topology, double fields, duplicated
face quadrature histories and map nodes have significant overhead.

## Numerical contract

Each hexahedral leaf has 12 electric edge values and six oriented magnetic
face flux histories. Electric unknowns are field amplitudes; line integrals
include the actual edge lengths. On the coarse/fine interface, fine
tangential values follow the neighbouring coarse hexahedron's edge basis:
constant along the edge, linear in the other face coordinate. At box edges
and corners this reduces to the shared coarse edge value.

Let P prolong independent edges to local leaf edges, C be the local oriented
face-edge incidence including edge lengths, D the diagonal local electric
quadrature and W the diagonal magnetic reluctance quadrature. The core uses

```text
M = transpose(P) D P       (electric mass, including off-diagonal terms)
G = C P
M dE/dt = transpose(G) W B
dB/dt   = -G E
```

The **same** constrained trace appears in G and its transpose. Electric mass
is not replaced by row-sum lumping after applying constraints. Away from an
interface it is diagonal. Connected interface blocks are SPD and are Cholesky
factored once at setup; timesteps perform triangular solves without allocation.
There is no iterative solve or user-selected residual tolerance in a half-step.

The face flux histories are local quadrature variables. Coincident faces have
separate storage but share the same constrained circulation; this is a deliberate
reference implementation overhead. `div(curl)` cancels per leaf, including
interface leaves. Arbitrary initial magnetic fields need not satisfy discrete
Gauss constraints; production sources/initialization require their own tests.

After diagonal normalization, write the equations as `S du/dt = K^T v` and
`dv/dt = -K u`. Symplectic leapfrog updates E then H. For a source-free run its
discrete invariant is

```text
0.5 (u^T S u + v^T v + dt v^T K u).
```

The code enforces a conservative bound below `dt * ||K S^(-1/2)||_2 = 2`,
using `sqrt(||K||_1 ||K||_infinity ||S^-1||_infinity)` and a 1.98 safety limit.
The inverse-mass bound is calculated at setup from the interface factors.
An explicitly unsafe timestep is rejected. Optional `TimeStepFactor` reduces
this coupled bound; native timestep-method selectors do not replace it.

This construction follows the hanging-variable/Galerkin line of work, rather
than claiming to reproduce a particular paper's complete algorithm:
[Venkatarayalu et al., 2007](https://doi.org/10.1109/TAP.2007.891809),
[Van Londersele et al., 2018](https://doi.org/10.1109/LMWC.2018.2840447).
The energy-based validation approach is also informed by
[Bekmambetova et al.](https://arxiv.org/abs/1705.02274).

## Native openEMS interface

Python, with coordinates in the base mesh's drawing units:

```python
fdtd = openEMS(CartesianSubgridBox=[-4, -4, -4, 4, 4, 4],
               CartesianSubgridRatio=[1, 2, 8], TimeStepFactor=0.8)
fdtd.SetBoundaryCond(['PEC'] * 6)
# Supply finite, strictly increasing Cartesian lines; spacing may vary.
# Add a zero-delay soft E box source and E/H probes or dumps.
fdtd.Run('simulation', engine='basic', numThreads=1)
```

Matlab/Octave uses its generic attribute serialization:

```matlab
FDTD = InitFDTD('CartesianSubgridBox', '-4,-4,-4,4,4,4', ...
                'CartesianSubgridRatio', '1,2,8', 'TimeStepFactor', 0.8);
FDTD = SetBoundaryCond(FDTD, repmat({'PEC'}, 1, 6));
```

Several separate blocks in Python (the example fits the -8:2:8 base mesh):

```python
fdtd.SetCartesianSubgrids([
    {'box': [-6,-4,-4,-2,4,4], 'ratio': [2,2,4]},
    {'box': [ 2,-4,-4, 6,4,4], 'ratio': [4,1,2]},
])
# Or append one region to an existing single- or multi-region setup:
# fdtd.AddCartesianSubgrid(box, ratio=[2,4,8])
```

`CartesianSubgrids=[...]` is also a constructor keyword. The list setter
replaces all regions; validation failure leaves the existing list intact.
Each entry requires `box` and defaults its own `ratio` to `[2,2,2]`.
The old singular setter explicitly replaces the list with one region.
`SetCartesianSubgridRatio` rejects multi-region mode as ambiguous. The
constructor rejects mixing singular and plural keywords.

Matlab/Octave provides a numeric matrix helper:

```matlab
FDTD = SetCartesianSubgrids(FDTD, ...
    [-6,-4,-4,-2,4,4; 2,-4,-4,6,4,4], [2,2,4; 4,1,2]);
```

It removes legacy single-box attributes before serializing the region list.
It checks box ordering, separation and ratios; the solver checks exact base
alignment, cell margins and resource guards. The helper has a Matlab test,
but Matlab/Octave execution is not part of the local macOS qualification.

For C++, `SetupCartesianSubgrids(boxes, ratios)` replaces all regions and
`AddCartesianSubgrid(box, ratio)` appends one. The compact core accepts
`CartesianSubgrid(lines, epsilon, mu, vector<Region>, uniformFine)` accepts
three vectors of physical coordinates, including nonzero origins. Existing
constant-spacing and single-box constructors remain source compatible and
construct their own base lines. `BaseLines()` preserves supplied coordinates;
`LeafStart()` / `LeafLengths()` expose actual physical leaf geometry. Integer
Leaf coordinates still index subdivisions of base intervals; dividing them
by LatticeRatio alone does not give physical coordinates on a nonuniform mesh.

Multiple regions use an optional FDTD child (one-box serialization remains
compatible with the legacy attributes):

```xml
<CartesianSubgrids>
  <Region Box="-6,-4,-4,-2,4,4" Ratio="2,2,4"/>
  <Region Box="2,-4,-4,6,4,4" Ratio="4,1,2"/>
</CartesianSubgrids>
```

Mixed legacy/list definitions, duplicate containers, unknown Region attributes,
empty lists, malformed coordinates, invalid ratios, touching, overlapping and
nested boxes are rejected. `ReadFromXML` resets previous region state. A
single Region without Ratio defaults to 2/2/2.

Python also supports `SetCartesianSubgrid(box, ratio=[2,4,8])` and the
separate `SetCartesianSubgridRatio([2,4,8])`. Changing the box alone preserves
the selected ratios. Ratio and box inputs are checked before a combined setter
changes either setting. A ratio of 1 is a valid no-refinement direction;
`[1,1,1]` is a regular coarse-grid control.

The C++ setters are `SetupCartesianSubgrid(vector<double>)` and
`SetupCartesianSubgridRatio(vector<unsigned int>)`. The compact core accepts
an optional final `Index ratio` argument after `uniformFine` (default 2x2x2).
Leaf sizes are now per-axis integer extents on that selected finest lattice;
electric edge length is still a scalar along its own direction.

XML uses optional FDTD attributes `CartesianSubgridBox="x0,y0,z0,x1,y1,z1"`
and `CartesianSubgridRatio="rx,ry,rz"`. A legacy file without a ratio retains
2x2x2. A ratio attribute requires a box. Malformed tokens,
trailing text and extra coordinates are rejected. Without this attribute the
normal operator and backend selection are unchanged. `Cartesian_Subgrid.py`
in `python/Tests` is a complete small model.

E/H probes and TD/FD E/H dumps sample reconstructed leaf fields through the
native processing interface. Their output sampling coordinates are still
the **base observation grid**, not an automatically expanded fine output grid.
The energy monitor integrates compact states; it does not estimate energy by
counting coarse observation samples. Soft E drive is added between the E and
H half-steps. A source that touches no unconstrained master edge is rejected.

Currently supported: positive isotropic, lossless background and cell-aligned
material boxes; outer PEC/PMC combinations; zero-delay nonperiodic soft E box sources; E/H observations.
Internal PEC and restricted one-edge resistive lumped ports follow the contracts below.
Unsupported requests fail at setup: distributed/LC/modal ports, conducting
sheets, dispersion/loss, PML/Mur, periodic sources, transformed
primitives, other field types and GPU/SSE/new CPU engine selection. Only
`engine=basic` selects this reference adapter; its execution is serial.

### Lossless material boxes

Use ordinary CSXCAD material properties, with no additional XML attributes:

```python
slab = csx.AddMaterial('slab', epsilon=3.556, mue=1.0)
slab.AddBox(start=[0, -8, -8], stop=[8, 8, 8], priority=2)
```

The native adapter resolves CSXCAD priorities at physical leaf centres after
checking that each intersecting box contains that whole leaf. Box faces may
follow fine local boundaries without lying on base lines, and dielectric
regions may cross coarse/fine interfaces. Boxes may overlap by priority.
Even a partially hidden box must be leaf-aligned: this conservative prototype
does not attempt Boolean subtraction to waive that restriction. Curved,
transformed, zero-volume and subcell objects are rejected, as are anisotropy,
spatial/non-unit weights, nonpositive/nonfinite coefficients and loss.

The core stores absolute SI epsilon/mu per leaf. Electric quadrature uses
that epsilon in the complete `P^T D P` mass, including interface terms;
magnetic reluctance uses the corresponding permeability. The coupled stability
bound is rebuilt from these coefficients. No material callback, geometry query
or extra allocation is added to the timestep update. Per-leaf material storage
adds 16 bytes per leaf on the tested platform; reference storage accounting
includes it, without claiming production RSS/VRAM savings.

Material-normal fields may jump. Reconstruction uses the positive-coordinate
leaf at internal cell faces (the last leaf at the domain's upper wall).
Roundoff-sized coordinate differences are treated consistently during lookup,
using a tolerance capped at `1e-8` of the local interval. The actual geometry
and observation position are unchanged; this does not relax material or
refinement alignment checks. Intentional offsets on either side remain distinct.
Probe integration and physical ports across discontinuities remain unqualified.

### Internal PEC geometry and contacts

Ordinary CSXCAD `Metal` properties now constrain local electric traces. No
new Python/Matlab geometry API or XML attribute is needed. Supported primitives
are axis-aligned volume/sheet boxes, solid cylinders and cylindrical shells.
Transforms, tilted cylinders, other curved primitives and wire/point boxes
remain unsupported. Volume/sheet boxes must follow actual leaf boundaries;
cylindrical end planes must align for all leaves overlapping their transverse
bounding box. Hidden boxes still require alignment. A setup with metal but no
resolved PEC edge fails explicitly.

Cylinder/shell radial contours use native-style edge-centre CSXCAD sampling,
including material/metal priorities. This is Cartesian EM sampling of the
continuous primitive, not a new fabrication contour or a conformal curved-cell
method. Thin walls, bores and gaps require actual occupancy/contact inspection
and mesh convergence. At least one PEC edge somewhere is not proof that every
primitive or contact is resolved.

The core eliminates PEC master E unknowns before forming curl and electric
mass; the same remaining trace is used in both. A dependent fine edge is
accepted only if its sampled PEC state agrees with the eliminated coarse trace.
An unresolved fine strip or an opening that the coarse side would close fails
at setup. Moving an interface, changing the base mesh or refining the whole
contact neighborhood may be required. This is a deliberate prototype restriction,
not support for arbitrary PEC/interface crossings.

`LeafEdgeIsPEC` exposes the actual local classification for independent audits;
`PECTraceCount` counts unique classified local edge keys, not independent
degrees of freedom. Magnetic flux rows enclosed by zero electric traces and
touching PEC retain zero initial normal B. `SetMagnetic` projects these rows
to zero; trapped initial flux is not part of this reference contract. The
time-step update adds no PEC lookup, allocation, geometry callback or new branch.
Reference set/face-mask storage remains accounted, not a measured speedup.

Qualification includes finite 18-um top/ground boxes connected by a hollow
barrel: all 31,488 local-edge samples agree with an independent geometric
classification. The connected wall is constrained while its bore remains open;
inserting an 18-um axial gap restores a free electric edge. This is a geometric
contact control, not a physical-port resistance or S-parameter result.

PEC boxes inside and across a refinement boundary, and a full sheet crossing
the boundary, preserve the source-free discrete invariant in 100,000-step tests
with a `1e-8` gate. A sheet splitting a cavity prevents field transmission into
the unseeded half. The remaining half-cavity frequency error falls from 3.11%
to 0.784% at base resolutions 6/12 (gate 1%, with improvement required).
Uniformly refined volume, sheet, cylinder and shell controls agree with regular
native Yee after 200 steps within `8.1e-7` of E/H peaks (gate `2e-4`).

Native integration controls cover internal field exclusion, a closed wall,
transmission through a higher-priority dielectric opening, unit/origin and
property-insertion invariance, sources covered by PEC and explicit rejection
of unsupported/unresolved geometry. E/H observations remain on the base grid;
dual-H snapping can put a nominal near-surface observation inside the metal.
Place observation stencils clear of the surface and inspect snapped coordinates.

### Local resistors and U/I integration

The serial adapter accepts `AddLumpedPort` with positive finite resistance when
its line spans **exactly one unconstrained independent local E edge**. Its two
transverse dimensions must be zero. This includes a local edge inside a refined
region that does not exist on the base observation grid. The existing zero-delay
soft E waveform/amplitude convention is retained; it is not reinterpreted as a
Thevenin voltage or Norton current. A zero-area `caps=True` line creates no
additional metal caps. Supply the actual electrode geometry separately.

Distributed boxes, unresolved/multiple longitudinal edges, hanging-only edges,
overlapping native resistors, series topology and additional L/C values are
rejected. This deliberately limits the native mapping rather than assigning an
unverified effective resistance to a multi-edge box. It is not yet the normal
finite-width microstrip-port workflow used by the antenna.

Voltage paths and rectangular current contours are prepared on actual leaves
once at probe definition and evaluated from cached sparse weights. U integrates
the constrained electric edge basis. I integrates the piecewise-constant **dual H
chain**, with segments split at local primal and dual boundaries. Integrating the
linearly smoothed H point reconstruction instead attenuates a Yee current
contour; point-field interpolation and circulation are distinct operations.
Zero-width current-probe spans expand to local dual bounds, tolerating only
coordinate roundoff. Output includes the actual physical integration bounds
before the normal column-name header; legacy base-index headers are informational.
E/H field probes and dumps retain their previous base-grid sampling contract.

The C++ reference also supports axis-aligned voltage-path resistors with optional
Norton current. If `V=b^T u`, `G=1/R`, and `C_p=sqrt(G_p)b_p`, midpoint coupling
adds `C C^T` to the normalized electric update. A setup-factored small Woodbury
system uses the existing interface-mass inverse; it never inverts a dense field
matrix. The fixed timestep and a maximum of 64 reference resistors are checked.
For the leapfrog invariant W, without the separately accounted soft E kick,

```text
W_new - W_old = dt * sum_p V_mid,p * (I_source,p - G_p * V_mid,p).
```

`SourceWork()` and `ResistorLoss()` accumulate source work and passive loss since
`SetResistors`; include the initial W when checking a nonzero initial state.
Soft-E work is accounted separately in the same ledger when resistors are
configured. Configure resistors before `SetSources`. Integration weights,
responses and factors are setup allocations; there is no per-step allocation.
A local-field scratch vector accounts for soft-source work; this reference is
not a production memory/performance claim.

Controls are in `tests/test_cartesian_subgrid_ports.cpp` and
`python/Tests/Cartesian_Subgrid_Ports.py`. The native fixture is a short PEC
microstrip with lossless epsilon-r 2.25 and two **50-ohm lumped terminations**;
it is not a calibrated 24-GHz laminate or a claim of ideal characteristic Z.
The substrate/port height is 0.2 mm. Frequencies are 20–28 GHz in 0.2-GHz steps,
with identical 50-fs timesteps and 24,000 steps. The controls preserve actual
point-port cell sizes: shrinking the point-port cell changes its contact and
local capacitance, even if its zero-width CAD line is unchanged.

| Control | max complex Z error | max complex S11 / S21 error | max S21 phase error |
|---|---:|---:|---:|
| Identical base grid, normal Yee / compact identity | 2.62e-6 ohm | 2.95e-8 / 4.54e-8 | not separately gated |
| Ports outside refinement, ratio 2/1/1 | 4.45e-4 ohm | 5.02e-6 / 3.72e-6 | 0.000212 deg |
| Ports inside refinement, ratio 2/2/1 | 0.0583 ohm | 4.06e-4 / 9.12e-5 | 0.00403 deg |

The outside control uses a conventional tensor Yee grid with the same fine X
interval and unchanged port cells. The inside control uses uniformly finer XY
Yee lines with the same local port cells; Z remains unchanged. Interface gates
are 0.5 ohm, 0.005 complex S and 0.1 degree. The maximum error in
`abs(S11)^2+abs(S21)^2=1` is 1.82e-5 outside and 1.78e-6 inside (gate 1e-3).
Polarity reversal leaves Z/S unchanged. Core controls include voltage additivity,
H orientation, coordinate-roundoff contours, parallel loads and rejected inputs.
A three-resistor XYZ-refined control, including a voltage path crossing the
interface, preserves injected-minus-dissipated-minus-stored energy within
1.46e-12 over 100,000 steps (gate 1e-10), with no passive energy increase.
These are scoped CPU controls, not arbitrary-interface, antenna, GPU or PML
qualification.

### Outer PEC/PMC boundaries

Use the existing Python or Matlab `SetBoundaryCond` workflow;
no new CSXCAD geometry or XML attributes are introduced. The native six-face
boundary codes retain PEC=0 and PMC=1. For example, in Python:

```python
fdtd.SetBoundaryCond(['PMC', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC'])
```

A PMC face retains tangential E unknowns and their one-sided physical mass;
the shared curl transpose enforces zero tangential H weakly. PEC eliminates
tangential E as before. This changes setup topology and coefficients only;
there is no new boundary branch, ghost lattice or allocation in the update loop.
The region-to-boundary base-cell margin is unchanged.

**Boundary planes are the physical outer base faces.** Ordinary Yee's PMC
planes are staggered on its magnetic lattice. Reusing identical outer line
coordinates therefore need not describe the same physical cavity. The native
identity control compares the half PMC domain against an independently
reflected, full ordinary-Yee PEC domain with matched sources and observation
positions. It does not silently count that staggered-plane offset as solver
error. Refined half domains compare against matching two-region full domains.

`Electric`/`Magnetic` remain local leaf reconstructions. At the exact PMC face,
tangential E is a well-defined constrained-space trace; normal E and tangential
H return unilateral cell samples, not symmetric boundary averages or imposed
point traces. The native controls compare tangential E on the plane and all
E/H components at three interior points. Do not place an NF2FF surface on that
boundary or infer exact pointwise boundary fields from these reconstructions.
Boundary-field reconstruction and broader radiation/absorber compatibility
remain separate qualification tasks.

All six face selections, mixed PEC/PMC and all-PMC analytical cavities, and a
100k-step two-region anisotropic-ratio invariant are covered by
`tests/test_cartesian_subgrid_boundaries.cpp`. Native tests in
`python/Tests/Cartesian_Subgrid_Boundaries.py` cover both half-domain sides
along X/Y/Z, with and without refinement. PML and Mur continue to fail at setup.
This is scoped serial CPU reference functionality, not absorbing-boundary,
antenna or GPU qualification. The six identity mirror controls agree with ordinary
Yee within 1.31e-6 of E peaks and 7.33e-7 of H peaks (gate 2e-4). Local and
graded mirror controls agree within 1.67e-15 of field peaks. The mixed-cavity
frequency error decreases from 0.2553% to 0.0640% on doubling base resolution;
all-PMC error decreases from 0.8569% to 0.2164%. The 100k-step two-region
invariant drift is 5.78e-15 (gate 1e-9). These values do not establish general
interface scattering or production performance.

## Tests and observed accuracy

```sh
cmake -S . -B build -DOPENEMS_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_subgrid --long
build/tests/test_cartesian_subgrid_mesh --long
build/tests/test_cartesian_subgrid_material --long
build/tests/test_cartesian_subgrid_pec --long
build/tests/test_cartesian_subgrid_ports --long
build/tests/test_cartesian_subgrid_boundaries --long
build/tests/test_cartesian_subgrid --mesh > mesh.json
python3 python/Tests/run_testsuite.py --unittests Cartesian
```

The C++ test verifies six-face tangential continuity, cancellation of divergence
of curl, compact allocation, and uniform-grid agreement with the native basic
Yee engine. At 200 steps the current three-component seed gives peak-normalized E/H
differences below 4e-6 for both cubic and anisotropic uniform grids (test
tolerance 2e-4, allowing float native fields versus this double reference). A random three-component field state
preserved the discrete invariant through one million steps within 1e-8.

For a fixed physical cavity, increasing the base resolution from 6 to 12 cells
per axis reduced analytic frequency error from about 1.10% to 0.28%; cyclic
coordinate permutations also agree. This is a mode/convergence check, not
antenna bandwidth or radiation qualification.

The PEC-compatible Gaussian-packet check compares time records with a fully
fine reference at an identical timestep. For base resolutions 6, 12, 24, the
peak-normalized errors are about 25.8%, 16.3%, 6.93% in each axis; at resolution
36 the X case is about 3.01%. Tests require monotone improvement, <10% at 24
and <5% for the resolved 36 case. **Coarse packet results are not accurate
enough for design.** No arbitrary-angle scattering or microstrip accuracy
claim follows from these tests. Earlier under-resolved packet experiments
are retained separately as development evidence, not relabelled as passes.

## Resource guards and remaining work

The prototype accepts powers of two from 1 through 16 independently per axis;
this ratio cap is an implementation/validation guard, not a physical maximum.
A large direct resolution jump can increase interface error and shrink the
common timestep, even if it conserves the discrete energy. Conservation and
allocation success are not accuracy qualification. Staged/nested transitions
and local timesteps remain future work.

The prototype has explicit allocation guards: at most 1024 base cells per
axis, one million base cells, one million leaves summed across all regions,
64 regions and 1024 coupled edges per
dense interface factor. These are implementation safeguards, not physical
limits of subgridding. Increasing them requires memory accounting and a
sparser/scalable interface solve, rather than merely removing assertions.

Additional ratio checks cover 1x1x1, 1x2x8 and axis permutations, 2x4x8,
4x4x4 and 1x1x16: leaf counts, six-face traces, div(curl) and energy invariants.
Uniform 1x2x4 (and permutations) recovers native Yee E/H propagation. A
1x2x4 packet family is compared against its own anisotropic uniform fine
reference at the same timestep, not against the cubic 2x2x2 dataset. In cyclic
axis permutations its error falls from 2.46% to 0.81% as the base resolution
goes from 6 to 12; the 2x4x8 cavity error falls from 1.18% to 0.30%. A
100,000-step 2x4x8 invariant check drifts by less than 1e-8 (observed ~5e-15).
Python tests cover ratio validation, XML round-trip/legacy reset and an actual
unequal-ratio E/H run. These are small-model controls, not a blanket accuracy
claim for all factor combinations or large interfaces.

The multi-region tests use two different anisotropic blocks: counts, local
sizes, all twelve box faces, div(curl), a 100k-step invariant and invariance
under region-order reversal. A 64-block control checks valid admission;
invalid overlap/nesting/touching and the combined leaf budget are rejected.
A packet traversing both blocks is compared with a fully fine control at the
same timestep. For base resolutions 6, 12, 24, 36, errors are about 61.0%,
37.2%, 9.42%, 4.37%; cyclic permutations agree at 6 and 12. Coarse results
remain unsuitable for design. Native multi-region soft-source/E/H probe and
FD dump checks plus API/XML round-trip/reset controls accompany the core tests.

Nonuniform controls (`test_cartesian_subgrid_mesh`) cover actual local geometry,
volume, twelve box faces, div(curl), coordinate-origin translation and long-run
invariant stability. An anisotropic, fully refined graded grid agrees with
native basic Yee after 200 steps within 2e-4 peak-relative gates. A valid thin
interval reduces the common step; malformed and unrepresentable grids fail.

The graded two-region packet errors at 6/12/24/36 base cells per axis are
58.1%, 36.6%, 10.9%, 5.21%. The 36 case initially missed the 5% gate and is
not relabelled as sufficiently accurate. Increasing propagation resolution to
42x36x36 reduces error to 3.77%, retaining the 5% resolved-case gate and the
coarser convergence evidence. This is a numerical control, not an antenna
accuracy claim. Native E/H source/probe and dump checks also compare a graded
fixture with its translated equivalent; broader physical-port qualification remains pending.

Material controls include an analytic two-dielectric TE cavity: the equation
`k1*cot(k1*a) + k2*cot(k2*b) = 0` imposes tangential E/H continuity for equal
mu. With epsilon-r 1 and 2.25, base resolutions 6 and 12 reduce the fundamental
frequency error from 0.939% to 0.238% (fine gate 1%, with improvement required).
The material interface intersects refinement faces. A 100,000-step random
heterogeneous epsilon/mu invariant control stays within the `1e-8` gate.
These are stability/dispersion controls, not broadband reflection qualification.

A graded fully refined dielectric slab agrees with native basic Yee after
200 steps to about `3.94e-7` / `1.02e-6` of E/H peaks (gate `2e-4`).
Independent native occupancy checks cover drawing units, nonzero origin and
overlapping priorities. Python controls reverse insertion order, translate the
same geometry, change drawing units, exercise local-only alignment and reject
unsupported media/geometries. They also guard field reconstruction at a
discontinuous-mu interface; the initial priority-order and boundary-roundoff
failures were in this experimental adapter/core, not established upstream bugs.

Before antenna use: broader material/PEC geometry and interface qualification, absorbing boundaries,
distributed/modal port and source qualification and broader scattering tests.
Before practical acceleration: structured interior kernels, compact interface
metadata, memory/RSS and interleaved native benchmarks. Then device kernels,
domain partitioning and communication can reuse this numerical contract.
Nested refinement levels and local subcycling are separate changes. No
running or frozen production model should adopt this prototype silently.
