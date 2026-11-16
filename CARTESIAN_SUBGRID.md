# Experimental compact Cartesian XYZ subgridding

`CartesianSubgridBox` replaces selected base cells with `rx * ry * rz`
children. `CartesianSubgridRatio` independently selects the refinement on
X, Y and Z: each factor is 1, 2, 4, 8 or 16. The default `[2,2,2]` retains
the original eight-child behaviour; `[1,2,8]` makes 16 children while leaving
X unchanged, and `[2,4,8]` makes 64 children. Spacing on each axis is divided
by its selected factor inside each selected box. Several separate boxes
can have different axis ratios; outside them the coarse cells remain. This is a separate Cartesian numerical core, not
angular interpolation or a slab copied from cylindrical multi-grid.

It contains a **serial CPU numerical reference** and scoped shared device
operations, qualified locally on the reference backend and physical Metal.
A separately pinned [single RTX3080 lane](CARTESIAN_SUBGRID_CUDA.md#qualified-single-rtx3080-lane)
also has scoped evidence; it predates the later CAD/larger-interface changes.
Both are for small validation models. This is not yet a replacement for the
production antenna engine. One-edge device RF ports, scoped device PML and
public vacuum PEC dipole controls are available; general accuracy/scalability
and representative PCB qualification remain open.

[Implementation and qualification checklist](CARTESIAN_SUBGRID_CHECKLIST.md)
tracks the completed controls and remaining feature/engine milestones in PR #29.

The [spatial-curl block increment](CARTESIAN_SUBGRID_BLOCKS.md) adds actual
structured interior addressing on the existing device backends with exact shared
coefficient pairs and sparse interface exceptions. Its full replay/scaling and
representative antenna qualification remain separate from focused controls.

The [native smooth-conductor/path-port candidate](CARTESIAN_SUBGRID_CONDUCTOR.md)
adds scoped volume-box mapping. Its implementation and completed scoped r03 CPU/Metal/CUDA evidence
are distinguished there; general PCB RF acceptance and
structured scaling remain open.

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
storage in the current banded reference is 1,105,364 versus 4,502,212 bytes
**for this same reference core**. These
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
The factors use their exact lower band in the existing breadth-first edge
order, with one reusable block-sized local RHS; see the measurements below.
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
material boxes; outer PEC/PMC combinations, scoped parallel Mur faces and
classical convolution PML with a checked homogeneous collar;
zero-delay nonperiodic soft E box sources; E/H observations.
Internal PEC and restricted one-edge resistive lumped ports follow the contracts below.
Unsupported requests fail at setup: distributed/LC/modal ports, conducting
sheets, dispersion/loss, mixed Mur/PML and adjacent Mur corners, periodic sources, transformed
primitives, other field types and SSE/new CPU engine selection.
`engine=basic` selects the full scoped serial reference. `engine=gpu-reference`
and `engine=gpu` have the narrower [device contract](#scoped-native-device-execution).

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

### Native planar and extruded PEC polygons

The native adapter now also accepts CSXCAD `Polygon` and `LinPoly` metal
primitives. It preserves their supplied vertices and uses native membership;
it does not replace a bend by rectangle unions or smooth a rasterized contour.
The Cartesian solver still samples that contour on its selected cells. This
is not a conformal-boundary method or fabrication-contour qualification.

```python
trace = csx.AddMetal('finite_trace')
trace.AddLinPoly(points=[[0, 2, 2, 0], [0, 0, .4, .4]],
                 norm_dir=2, elevation=.27178, length=.018, priority=10)
```

Coordinates use the configured drawing unit. Normals X/Y/Z, reversed vertex
order and negative extrusion length are covered. Finite planar vertices and
nonzero extrusion height are required. Transforms, non-Cartesian coordinates,
rotated-body polygon subclasses and polygonal material volumes remain rejected.
Axial ends must align with local cell boundaries; the conservative check uses
the object's transverse bounding box. Existing coarse/fine PEC trace guards
remain: an unresolved gap/contact is not silently closed.

A thin-box control exposed an adapter alignment error independently of polygon
support: undoing the translated physical coordinate could cancel most of an
8mm origin when reconstructing a 0.018mm face. The bounded roundoff estimate
now includes the cancelled origin and relative-coordinate operands. Its
`1e-8` local-cell cap remains; real axial displacements of `±1e-9mm` are still
rejected. Geometry is not snapped. Native CSXCAD membership retains the ordinary
operator's drawing tolerance; the correction does not change CSXCAD itself.

[Native regression](python/Tests/Cartesian_Subgrid_Polygons.py) covers 15 shape/
axis/unit groups on basic/reference/physical Metal, resolved field exclusion,
18µm strip/ground port contacts, complex Z/S, phase and passive power, eight
thin alignment controls and ten strict negative controls. Generic field
differences stay below `1.56e-5`; Cu18 port differences stay below `6.37e-6Ω`
and `6.26e-8` complex S, with power-balance error below `1.81e-5`. Numerical
gates remain `2e-4` fields/S, `0.01Ω` Z and `1e-3` power. These are adapter
controls, not a calibrated 24GHz manufacturing model or mesh-convergence study.

The exclusion witness keeps the **whole native probe interpolation footprint**
inside copper. A raw vector probe near a thin edge can sample neighbouring
staggered base positions outside metal; its nominal point alone does not prove
that all sampled components are interior. Output coordinates remain the base
observation grid, without an implicit global fine output lattice.

Qualification used a new isolated prefix
`output/openems-subgrid-polygons-qualified-20261115/install`, configure-time
parent `a846ebf`, core SHA256
`62ca999ad42efba6ac4b48387e0a9ab626f2d2344d72f650664581117d13ebe7` and
NF2FF SHA256
`3080a152bc7d32b2cb4d2c15b5669100f04cedfd0a0c35b597b8d18247400e4a`.
All 17 qualification stages pass: 21 CTest controls, 14 Python unit modules,
prior CPU/reference/Metal field/port/DC/PML/dipole and ordinary-engine controls,
plus the new polygon regression. All 11 earlier core/NF2FF prefixes remain
unchanged, and 11,255 prior numerical arrays exactly match the frozen address-
index runtime. This CAD increment is not physically qualified on CUDA, with
lossy conductors or a representative PCB antenna; separate pre-CAD CUDA tests
do not transfer those qualifications automatically.

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
along X/Y/Z, with and without refinement. Adjacent Mur corners still fail at
setup; Mur and PML follow the separate contracts below.
This is scoped serial CPU reference functionality, not absorbing-boundary,
antenna or GPU qualification. The six identity mirror controls agree with ordinary
Yee within 1.31e-6 of E peaks and 7.33e-7 of H peaks (gate 2e-4). Local and
graded mirror controls agree within 1.67e-15 of field peaks. The mixed-cavity
frequency error decreases from 0.2553% to 0.0640% on doubling base resolution;
all-PMC error decreases from 0.8569% to 0.2164%. The 100k-step two-region
invariant drift is 5.78e-15 (gate 1e-9). These values do not establish general
interface scattering or production performance.

### First-order Mur absorbing faces

The existing `MUR` boundary code (2), Python/Matlab `SetBoundaryCond` and
optional Mur phase-velocity setting are reused. This increment implements a
paired `Operator_Extension` / `Engine_Extension` on the serial compact engine;
it does not change ordinary Yee Mur. For example:

```python
fdtd.SetBoundaryCond(['MUR', 'MUR', 'PEC', 'PEC', 'PMC', 'PMC'])
```

One face or two opposite, parallel faces are supported. Adjacent Mur faces
and mixed Mur/PML configurations are rejected. Each Mur face needs
two unrefined base cells of homogeneous lossless material, free of internal
PEC, enabled sources and lumped resistors. The collar covers the entire face,
not just a local refinement neighbourhood. Setup rejects violations. This
strong guard keeps both boundary and first interior electric traces independent
and their mass diagonal. Interior local refinement remains permitted.

Mur lives on the outer electric base plane, with the neighbour one base line
inside. For each tangential voltage trace:

```text
a = (v*dt - d)/(v*dt + d)
V_boundary(new) = V_inside(old) + a*(V_inside(new) - V_boundary(old))
```

Here `d` is the actual outer interval, including graded meshes; the default
`v=1/sqrt(epsilon*mu)` uses the collar medium. Histories and trace weights are
allocated once. Ordered E/H half-sweeps expose the ordinary extension hooks;
the existing full-step path retains its arithmetic. Invalid half-sweep ordering
or changing dt between E and H is rejected. There is no device fallback:
device absorber extensions and SSE/new CPU subgridding still fail explicitly
at setup. The separate lossless closed-boundary device path is described below.

`BoundaryWorkOut()` records **signed discrete correction work**, not far-field
radiated power and not a guarantee of positive dissipation on each step.
The leapfrog invariant plus this work closes the source-free control ledger.
This is not a passivity proof for arbitrary media, fields or interface layouts.

Controls in `tests/test_cartesian_subgrid_mur.cpp` cover all six face normals,
a local-refinement packet, full versus split-step equivalence and ordering
errors. Identity pulse fields agree with ordinary Yee within 1.59e-7 of initial
amplitude (gate 2e-4). The reflected packet peak is 0.00522 (identity) and
0.00580 (local). After 100k local steps, remaining energy is 3.35e-5 of its
initial value, with peak energy no greater than initial and a maximum corrected
ledger error 1.19e-13 (gate 1e-10). The nonzero local residual is retained;
this is bounded evolution, not complete asymptotic absorption qualification.

`python/Tests/Cartesian_Subgrid_Mur.py` exercises the actual scripting/XML
path, including scope rejection. A matched long reference moves the remote
PEC echo beyond the six-nanosecond observation window. Frequency-domain
reflection is extracted from the complete incident and difference signals;
four- versus six-nanosecond magnitudes must differ by less than 0.01. A narrow
10 GHz pulse avoids the long near-cutoff tail observed with wider illumination.
The oblique TE mode has fixed transverse wave number, so 30/45 degrees refer
only to 10 GHz, not the whole band. TEM is used for normal incidence.

| Incidence at 10 GHz | Measured amplitude reflection, 9.5 / 10 / 10.5 GHz |
|---|---|
| Normal | 0.00243 / 0.00269 / 0.00298 |
| TE 30 degrees | 0.0744 / 0.0711 / 0.0653 |
| TE 45 degrees | 0.226 / 0.149 / 0.128 |
| TE 30 degrees, graded normal intervals | 0.0732 / 0.0697 / 0.0636 |

The four identity waveforms agree with ordinary Yee within 1.03e-6 of the
incident peak (gate 2e-4). The largest 4-to-6 ns magnitude change is 9.56e-4.
The continuum first-order estimate `(1-cos(theta))/(1+cos(theta))` is recorded
for context; finite mesh, finite pulse and repeat returns need not reproduce
it exactly. The local-refined matched-reference normal packet return is also
checked. Oblique TM and adjacent corners remain unsupported for Mur; the
separate PML implementation below adds scoped absorbing-face/corner controls.
Wider converged spectra and general refinement-interface scattering remain open.
These reflections are too large to treat Mur as a precision antenna absorber.

### Classical convolution PML

Use the existing C++/XML/Python/Matlab boundary setting, for example:

```python
fdtd.SetBoundaryCond(['PML_8'] * 6)
```

This is a separate paired `Operator_Ext_CartesianSubgrid_PML` /
`Engine_Ext_CartesianSubgrid_PML`, not the ordinary Yee UPML implementation.
The ordinary solver and its PML are unchanged. The initial scope was serial
`engine=basic`; the later [device PML qualification](#scoped-device-pml-and-dipole-radiation)
adds existing reference/Metal backend execution. Any subset of PML faces, including edges/corners
where layers overlap, may be used with the remaining PEC/PMC faces. Mur/PML
mixing is rejected.

Each PML face requires at least four base cells, plus **one inner guard cell**.
Its full collar must be unrefined, match the homogeneous lossless background,
and contain no internal PEC, enabled source or resistor. The guards reject
refinement/material/geometry/port intrusion and overlapping opposing layers.
The collar's electric mass is diagonal; an additional check rejects PML
corrections that touch any coupled interface mass block. Interior lossless
materials, PEC, ports and independent refinement regions retain their previous
scope. No per-half-step host/device copies or device fallback are introduced.

The extension groups the existing normalized curl and transpose terms by
derivative axis. It stretches each derivative independently, at the physical
electric-edge or magnetic-face location. For rate `r` in inverse seconds:

```text
b = exp(-r*dt)
a = expm1(-r*dt)
psi_new = b*psi_old + a*D(field)
D_stretched = D(field) + psi_new
```

This is the unit-kappa, zero-frequency-shift limit of recursive-convolution
PML, not a configurable CFS-PML. The cubic profile is
`r = -2*c*log(1e-8)/L * (depth/L)^3`, using the actual physical layer width
and background wave speed. `1e-8` is a continuum round-trip design parameter,
**not a measured discrete reflection**. No runtime dependency or new CSXCAD
attribute is added. The implementation was independently derived from the
stretched-derivative formulation; see
[Roden and Gedney (2000)](https://doi.org/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A).

Pre-update hooks advance histories from the old opposite field; post-update
hooks calculate targets, and apply hooks write them before the dependent
half-step. The lossless core arithmetic is unchanged. Histories/scratch are
allocated at setup. `BoundaryWorkOut()` records signed field-correction work,
not a nonnegative material loss or far-field radiation. In the source-free
controls the original leapfrog invariant plus this work closes to the stated
tolerance. This ledger is not a general PML passivity/stability proof.

`tests/test_cartesian_subgrid_pml.cpp` covers all six normals, a graded/local
packet, PEC negative control, and six-face overlap/corners with one or two
separate anisotropic refinement regions. The 100k-step local packet retains
3.42e-5 of its initial field energy, with peak energy no greater than initial
and a maximum corrected-work residual 7.22e-12 (gate 1e-9). Its nonzero residual
is retained. The separate 10k-step corner control and native evidence are
recorded with actual source/build hashes in the contributor's evidence archive.

`python/Tests/Cartesian_Subgrid_PML.py` measures boundary-return amplitude
against an otherwise matched long compact reference, with the remote echo
outside a 6 ns window. Frequencies are 9.5/10/10.5 GHz. Oblique TE 30/45 and
TM 30 degree labels apply at 10 GHz only; transverse wavenumber is fixed.
The 4-to-6 ns window check must change amplitudes by <0.01, and the reflection
gate is 0.02. The PEC negative control must reflect >0.8 in amplitude.
This separates boundary returns from any scattering shared by the matched
coarse/fine interface. It is not general interface-accuracy qualification.
The ordinary-Yee comparison is an early-wavefront check (<0.7 ns), covering
only about 6.2% of the eventual incident peak. Across the complete normal
record the two different absorber paths differ by 1.31e-3 of the incident
peak; this is recorded separately, not described as full-waveform equality
or attributed to a solver defect. Existing lossless identity controls remain
the check of the compact core against ordinary Yee.

Observed 10 GHz amplitude returns are approximately 7.41e-5 (normal),
6.20e-5/5.02e-5 (TE 30/45), 6.19e-5 (TM 30), and 1.59e-3 (graded TE 30).
The local-refined normal control gives 7.42e-5. The larger graded-grid return
remains visible; no claim of reflectionless discretization is made. These
small 10 GHz controls do not qualify a 24 GHz antenna. Broad/near-cutoff/grazing
spectra and general material traversals into PML remain open. The later
dipole increment below separately qualifies a scoped NF2FF workflow;
lossy conductors and GPU/multi-device subgridding remain open.

Qualification build: macOS/arm64 Release, source parent `eaa1a31` plus this
increment, with isolated CPU-only prefix. The actual loaded core SHA256 is
`83682d094e5eb3aa073bb77be4158acd8e5bf2569e7448b1fdb2b4e144b5d893`;
matching CSXCAD SHA256 is
`f1485025adfe48c70aa8664db2dc34a0abdf55834620b7536d63e77cdd6da9f0`.
Before/after source inventories and runtime hashes accompany the qualification
logs. Embedded version text names the configure-time parent; source inventories
identify the tested increment. The frozen earlier reference library remains unchanged.
CTest 11/11, Python unit modules 13/13 and native subgrid groups 7/7 pass.
The native ordinary CPU/cylindrical controls also pass (2/2 groups).
Matlab PML serialization test source is present, but no local interpreter was
available, so it was not executed. No device execution or performance result
is claimed by this build.

### Scoped PEC dipole radiation

`python/Tests/Cartesian_Subgrid_Radiation.py` qualifies the existing serial
`engine=basic` path for a **1 GHz solid rectangular PEC dipole in vacuum**.
This increment adds a regression and qualification evidence; it changes no
numerical engine or NF2FF algorithm. It is not the project's 24 GHz PCB,
a dielectric/lossy antenna qualification, or a speedup result.

The two 20 × 20 × 60 mm PEC arms occupy Z −65…−5 and +5…+65 mm.
One 50-ohm source/resistor occupies the independent local Z edge across the
10 mm gap. An independent native-XML audit checks both finite terminal
contacts, identical source/resistor/voltage spans, current-plane position,
direction and signs. A broken contact must fail that audit. Passive arm
currents use symmetric dual Z planes; their actual positions are recorded
per mesh and are not compared across different physical observation planes.

The refinement box is `[-30,-30,-75,30,30,75]` mm with ratio `[2,2,1]`.
The base X/Y spacing is 10 mm; the accepted interior Z spacing is 5 mm,
except for the preserved single 10 mm port gap. A matching tensor grid
extends the local 5 mm X/Y lines across the domain. A finer-interior control
also refines the non-PML transverse air. All controls preserve geometry,
port cells, timestep and physical PML inner boundary. This is a finite mesh
sensitivity check, not asymptotic convergence for arbitrary antenna models.

Inner/outer closed NF2FF boxes are ±80/±100 mm in X/Y and ±85/±105 mm in Z.
All six E/H faces and their node-interpolation neighbours remain in
homogeneous unrefined air, outside the entire local box and inside PML.
Actual HDF5 meshes, face coverage, frequencies and finite excited fields are
audited. A missing-face negative control must fail before postprocessing.
No mirrors, refinement-crossing surfaces or material-crossing surfaces are
qualified. The test compares two enclosures, native FD and TD NF2FF,
sphere-integrated far fields and surface flux, complex near fields,
directivity, realized gain, accepted/radiated power and passive arm symmetry.
The analytical thin-wire dipole pattern is a physical sanity check for these
finite rectangular bars, not their exact solution.

DC-free Gaussian excitation is sampled at 0.9/1.0/1.1 GHz with a common
5 ps step, over 20 ns and a 30 ns temporal control. The intentionally tiny
end criterion forces fixed duration; maximum-step warnings remain visible.
DFTs use the probes' actual timestamps, including the H half-step. Adding
another half-step causes about 1.57% error and must be detected. The maximum
FD/TD probe error is 0.047%; complex-vector NF2FF error is 0.0376% and
radiated-power error 0.0356%. Near-zero Ephi is normalized by the total or
co-polar field peak, not by its own tiny magnitude.

Observed maxima over the three frequencies and sampled angles:

| Comparison | Complex Z difference | Complex Gamma difference | Complex Etheta far-field error | Common-point E/H error |
|---|---:|---:|---:|---:|
| Compact tensor / matched ordinary Yee, PML10 | 0.00860 ohm | 5.57e-5 | 0.00560% | 0.00556% |
| Local / matching tensor | 0.235 ohm | 0.00232 | 0.249% | 1.43% |
| Local / finer-interior control | 0.555 ohm | 0.00292 | 0.347% | 1.90% |
| Local 20 / 30 ns | 2.73e-5 ohm | 2.29e-7 | 1.94e-5% | 3.82e-5% |

Fields are complex and normalized to incident port voltage. Both Ephi and
Etheta comparisons use the reference Etheta peak; relative errors near a
pattern null are not used. Local/tensor limits are 0.5 ohm, 0.005 in Gamma
and 2% in near/far fields, directivity and realized gain. The finer-interior
limits are 2 ohm, 0.02 and 5%; the tighter ordinary-Yee and temporal gates
are explicit in the test source. None was loosened after failed trials.
At 1 GHz the accepted local port is approximately `69.2296+j6.3180 ohm`.
Its inner-surface radiated/accepted-power ratio is 0.99697…0.99719 across
the sampled band (outer surface 0.99724…0.99795); this is numerical lossless
power closure, not a lossy antenna efficiency measurement. The realized-gain
local/finer-reference difference is below 0.466%. Pulse spectra are
consistently normalized; no absolute CW-watt claim is made.

Two failures are retained in the local evidence archive. The exploratory
near-zero Ephi normalization failure was a test error. The first formal
candidate used 10 mm interior Z: despite local/tensor agreement, its finer
reference changed Z by 3.9186 ohm and Gamma by 0.02643, failing the original
2-ohm/0.02 gates. A controlled Z-only refinement reduced this to the accepted
0.555 ohm/0.00292. The public test keeps the old coarse-Z case as an expected
negative under the same gates. This is mesh underresolution, not an
established upstream engine defect. PML6/PML10 controls keep the same inner
boundary and demonstrate absorber sensitivity of the ordinary-Yee comparison;
the two solvers use different PML discretizations.

Reproduce with bindings/libraries built from this branch:

```sh
python3 python/Tests/run_testsuite.py --keep --verbose Cartesian_Subgrid_Radiation
```

Ten cases include the expected-negative coarse-Z case. Actual qualification
uses an isolated macOS/arm64 Release CPU build from numerical source
`55fe3e2bf4be8807f9d136e66ac0e94e847a9d08`, plus this Python test.
Loaded core SHA256:
`e10e2a95b9508c64b2763c0f82a61dfaf5351293a8c025305d72dc0b16913f31`;
loaded nf2ff SHA256:
`41cfa2d3d34a91fa3e4ace6b0200fd349c88debeb5f64121235b526db9c64fda`;
CSXCAD SHA256:
`f1485025adfe48c70aa8664db2dc34a0abdf55834620b7536d63e77cdd6da9f0`.
Before/after inventories identify the actual test source and check both frozen
earlier runtimes unchanged; embedded version text identifies the configure-time
parent. Raw local evidence is retained in the contributor's
`experiments/openems-cartesian-local-subgrid-20261112/pec-radiation/` archive,
including failed inputs, numerical results, runtime hashes and figures.
CTest 11/11, Python unit modules 13/13, native subgrid groups 8/8 and ordinary
engine/field/radiation groups 4/4 pass on this actual build.
Matlab/Octave remains unexecuted. General CAD/ports/surfaces, the project array,
losses, device execution, resource savings and production adoption remain open.

### Banded interface factors

The serial reference stores each connected SPD interface factor in a lower
band instead of a full square. The existing breadth-first edge order is
unchanged. For a block of n edges, its graph determines half-bandwidth b;
Cholesky fill remains within that band. Every in-band coefficient and fill
term is retained, including zeros. `L(i,j)` is stored at
`i*(b+1)+i-j`, with the diagonal first. No mass lumping, iterative tolerance,
changed quadrature or dropped coupling is introduced.

Electric updates gather each block RHS into one preallocated local vector,
solve in place, then scatter to the same independent edges. Nonzero products
retain their previous ascending order. The workspace is the largest block,
not a global field vector. The conservative inverse-mass norm uses one signed
RHS for positive bipartite interface graphs, with the original n-column
calculation retained as the general fallback; see the next subsection.
The norm formula, 1.98 safety factor, common timestep and all existing
allocation guards remain unchanged. A step allocates no new workspace.

`InterfaceStorage()` exposes read-only C++ reference counts: blocks, largest
block/half-bandwidth, allocated factor entries, equivalent dense entries and
workspace entries. It adds no Python/Matlab input setting or XML attribute.
`StorageBytes()` includes the new workspace and band allocation. Its map-node
accounting remains approximate; it is not process RSS or peak memory.

The public [core benchmark](benchmark/README.md#compact-cartesian-subgrid-core)
compares unchanged dense parent `fb98953` with the candidate using the same
O3/C++17 compiler. Three interleaved pairs per case retain topology, stable
timestep, invariant and E/H checks. Selected medians on macOS/arm64:

| Core fixture | Constructor, dense/banded [s] | Update, dense/banded [s] | Accounted storage [MiB] | Peak process RSS [MiB] |
|---|---:|---:|---:|---:|
| Uncoupled 16³, 1000 steps |0.0116/0.0112|0.0819/0.0823|4.36/4.36|7.06/7.08|
| Local 24³, 500 steps |0.0893/0.0874|0.3091/0.2920|33.44/32.85|40.39/39.78|
| Two regions 16³, 500 steps |0.0580/0.0560|0.2007/0.1904|17.79/17.34|23.44/22.98|
| Strip 128×8×8, 1000 steps |0.1766/0.0788|0.8884/0.5245|25.72/21.57|30.58/27.16|
| Strip 512×8×8, 200 steps |9.2803/0.3775|2.2238/0.4249|151.32/86.72|176.61/112.02|

The largest fixture has 89,664 leaves, 229,044 independent E states and 537,984 local
H histories. Its largest block has 1024 edges and half-bandwidth 2; allocated
factor doubles fall 8,518,656→48,960, with 1024 workspace doubles. E/H differences
across all seven fixtures are at most 3.0e−15 of the reference peak. The
uncoupled control adds 24 accounted bytes; its timing difference is within spread.
Largest-case update spreads are 16.8% dense and 0.9% banded; constructor spreads
are 7.3%/5.4%. Raw repetitions are retained; these medians are not universal
performance bounds.

This measures the compact reference core, excluding native processing. It is
not a speed or resource comparison against optimized Yee/GPU execution or a
production antenna. Regular interiors still use sparse curl rows; topology
maps and duplicated H quadrature histories remain. A broad graph band may
give little saving; non-bipartite inverse-bound setup still needs n solves.
At this earlier milestone the 1024-edge cap was deliberately retained; see
[the later bounded extension](#bounded-larger-interface-factors). Larger/general interfaces, native
structured kernels, threading/SIMD and GPU execution remain separate work.

`test_cartesian_subgrid_factor` checks a frozen dense E/H oracle and timestep,
10,200-step invariant, block-local workspace, an accepted 1024-edge block and
rejection above that unchanged guard. Existing million-step reference,
100k-step material/PEC/port/boundary/absorber and 10k-step PML corner controls
are repeated. The ten-case native PEC dipole also reruns with identical
input hashes and is compared to its frozen dense-build complex port, near-field
and radiation arrays. Native ordinary-engine tests and the interleaved native
master comparison remain separate from the core timing table. All ten native
cases retain identical complex arrays (maximum normalized difference zero);
CTest 12/12, Python unit modules 13/13, long core invocations 9/9, native
subgrid groups 8/8 and ordinary groups 4/4 pass. The four XS/basic ordinary
master/candidate benchmark pairs have identical probe traces; their update
rate differences are classified as measurement noise, including the antenna
case with 24.5% spread.

Qualification build: isolated macOS/arm64 Release, source parent `fb98953`
plus this increment. Actual loaded core SHA256:
`8f4cef1025a9de8442389085f7c2fddfbfbec3df4f62e64e3783a2a86b629e91`;
nf2ff SHA256:
`e1ec4b39e47d4554f51a545e44a5605cfb2141e78a64368b684abf2b49e4554c`.
The source inventory identifies the tested files; embedded version text names
the configure-time parent. Prior reference/PML/radiation binaries are frozen
and checked before/after. Local source/build/test and benchmark evidence:
`experiments/openems-cartesian-local-subgrid-20261112/banded-interface/`.

### Signed inverse-norm setup

Positive trace weights assemble nonnegative interface off-diagonals. If their
undirected graph is bipartite, its two-colouring defines a diagonal sign
matrix D for which `D S D = I - B`, B entrywise nonnegative. Since S is SPD,
the Neumann series of `(I-B)^-1` is nonnegative. Thus
`abs(S^-1) 1 = D solve(S, D 1)`: one existing Cholesky solve gives exactly
the same absolute inverse row sums as all n unit columns, up to roundoff.
No coupling, quadrature, factorization or spectral-bound formula changes.
An odd cycle, or a nonpositive/nonfinite signed result, uses the old fallback.
The final `1 + 1e-12` inflation and allocation guards remain unchanged.

`InterfaceStorage()` additionally counts inverse-bound solves and blocks
using the sign transform. The focused regression compares against independent
dense inverses for graded chains/sheets and disconnected graphs, and rejects
an odd cycle before solving. The frozen dense-field factor control requires
one solve per current bipartite block. Seven standalone fixtures use three
interleaved pairs against source `571766f`: the largest strip's median setup
is 0.3690/0.3106 s (15.8% less), with at most 4.7e-14 normalized E/H difference.
Memory and update time are essentially unchanged. This is a setup change in
the small CPU reference, not a production-antenna or GPU speedup.

Isolated macOS/arm64 Release qualification passes CTest13/13, Python unit
modules 13/13, long core invocations 9/9, native subgrid groups 8/8 and ordinary
groups 4/4. The ten-case dipole's frozen complex port/field/radiation arrays
remain identical (maximum normalized difference zero). Source/runtime
inventories match across the run, including four earlier frozen prefixes.
The loaded candidate core SHA256 is
`2b8eea4f5f5b7f4bb44f44530eaa06f81997a9b586906c92548d6881a13a919b`;
nf2ff is `006d1b9834cf108cefd12e657c51109a269df902aac1550f41f7071bcdc00bb2`.
The existing Engine_Compare suite passes and a separate short control
explicitly records Metal on Apple M5 Pro with all extensions on device.
These ordinary-engine regressions do not qualify subgrid GPU execution.
No project antenna, conductor or multi-device adoption is included.

### Exact translated curl stencils

The CPU reference streams normalized curl rows into shared exact translated
patterns. Each row retains its base index and pattern ID; every coefficient
and relative column offset is compared exactly before reuse. There is no
quantization, dropped term or change in row/term order. Transpose updates and
classical PML consume the same patterns. Setup no longer materializes a global
curl CSR before compression, and releases temporary mass assembly before
allocating the fields. Topology maps and duplicated magnetic quadrature states
remain; this is an intermediate memory reduction, not structured Yee interiors.

An independent CSR oracle checks forward/transpose arithmetic, repeated and
unique rows, empty rows, one-ULP coefficient distinctions and malformed input.
Seven fixtures with three interleaved pairs against parent `cb36dde` have
identical retained E/H arrays and timestep/invariant checks. The 24-cell local
fixture uses 32.85/17.05 MiB accounted storage and 39.72/20.56 MiB peak RSS;
the 512-cell strip uses 86.72/51.18 MiB and 112.05/62.34 MiB respectively.
The memory saving has a CPU cost: updates change from 0.2896 to 0.3254 s and
0.4223 to 0.4598 s respectively. Setup is nearly unchanged. These small-core
measurements do not establish production speedup or memory requirements.

Isolated Release qualification passes CTest 14/14, Python unit modules 13/13,
long core invocations 9/9, native groups 8/8 (including the ten dipoles), and
ordinary groups 4/4. All frozen dipole port/field/radiation arrays are unchanged.
Source and four earlier frozen runtime inventories remain unchanged during
qualification. Loaded core SHA256:
`9db83252b72729aabaefc64135abae80f072cd54115b05cabf23b0b90d3a82d9`;
nf2ff: `3677d50a93a2d6d9a711fd105e7fb1fd21f0477b47c854ba848f49d4b1265835`.
Reproducible contributor evidence is in
`experiments/openems-cartesian-local-subgrid-20261112/translated-stencils/`,
including final benchmark `core-comparison-r05-final.json`. Earlier timing
attempts are retained but are not the final benchmark. This increment does
not qualify GPU, lossy conductor, PCB antenna or multi-device execution.

### Upstream 92cd92a integration

The reference is explicitly rebased by merge onto upstream
`92cd92a8d924c78a93453c25fa541384c2ac1d0e`. This adopts the common GPU cell
model/point-operation architecture and the upstream fixes for excitation
frequency initialization, rounded periodic signal length, XML timestep factor
and explicit zero end criteria. Subgrid periodic sources remain unsupported.

A separate isolated build passes CTest 14/14, Python unit modules 13/13, long
core invocations 9/9, native groups 8/8 and ordinary groups 4/4, including the
expanded CPU/reference-GPU/Metal engine comparison. Frozen runtime prefixes
and inventoried source remain unchanged during qualification.

The previous strict 2e-11 frozen-array identity check initially fails; this
failure is retained. Upstream now compiles engine extensions without fused
multiply-add. All ten dipole results remain within 6.45e-7 peak-normalized
error (existing upstream single-precision comparison gate 1e-5). Subgrid
port/near-point arrays remain unchanged; subgrid far-field differences are
below 1.89e-8. An isolated identity-dipole control restores only the prior PML
contraction mode and recovers every frozen array exactly. Its loaded library
path/hash is checked; the first control with an incorrect soname alias was
not accepted as causal evidence. The counterfactual is not an adopted build.
No frozen reference or numerical gate in the solver tests was overwritten.

Loaded merged core SHA256:
`8648669ef485b93f4502490cea44906b72a8b73f19c751adbf0c36f03845ad85`;
nf2ff: `3cee7e023b9d7387e273e7cd259ffcee56a47082e529613ce4420f9d3f6098f7`.
Source/runtime inventories, the retained failed strict comparison and causal
control are in the contributor evidence directory
`experiments/openems-cartesian-local-subgrid-20261112/master92-integration/`.
This is compatibility qualification, not a subgrid GPU or antenna speedup.

### Scoped native device execution

Compact subgrid operations use the existing backend `DeviceAlloc`, `RunPoints`
and `DeviceDownload` contract. Shared portable point operations implement curl
and its gather transpose, independent banded mass solves, electric updates and
ordered additive sources. Metal, HIP and the explicit CPU reference backend
compile the same operations. **Physical qualification here covers Metal on
Apple M5 Pro only; HIP/CUDA hardware qualification remains open.** No additional
backend or silent host half-step execution is introduced.

The native adapter accepts `engine='gpu-reference'` and `engine='gpu'` for its
lossless closed PEC/PMC field/source scope and the separately qualified
one-edge RF resistor scope below. It rejects Mur and any extension without a paired device implementation.
Classical convolution PML has the separately qualified device path below.
Device state and coefficients are FP32. The CPU reference remains double.
Electric/magnetic fields stay on device between sweeps; complete-step field
mirrors are refreshed once per processing batch for existing reconstruction,
probes and FD dumps. This is host DFT processing, not device DFT/gather parity.

The first upload expands curl/transpose rows and uploads existing banded
factors. It preserves local leaves without a global fine parent lattice, but
is not structured block indexing and does not remove current allocation guards.
Reported device allocation bytes exclude host topology/observation mirrors,
command buffers and backend overhead. Large-antenna VRAM/speed claims remain
unsupported.

The independent test-only FP32 oracle, separately implemented from the shared
point kernels, detects arithmetic/order/dispatch differences. Three seeded
uniform/local/graded source fixtures run 2000 steps; a source-free local case
runs 20000 steps. Physical Metal agrees exactly with this oracle in all retained
fields. Maximum peak-normalized difference from the double reference is
9.75e-5 (gate 2e-4); source-free invariant drift is 7.59e-8 (gate 5e-4).
Repeated backend creation/destruction is exercised; this is not a general
resource-leak or allocation-limit qualification.

Four native fixtures (local, graded unequal ratios, multiple regions, and
lossless dielectric plus resolved PEC) compare full E/H time records and
multi-frequency FD dumps. Maximum double-CPU/device difference is 5.30e-6;
reference-backend/Metal records are identical, within the declared 2e-4 gate.
Invalid device field dimensions, nonfinite fields, incomplete half-steps and
unsupported native device absorbers have rejection controls. The later PML
increment adds positive PML coverage and retains Mur as a rejection control. The first increment
also rejected ports; the subsequent scoped RF increment below replaces that
restriction with positive port tests. The native
Python binding returns setup code4 for engine-construction failure; tests check
that code, the actual native error reason and absence of propagated output.

Reproduction:

```sh
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_subgrid_device --metal
python3 python/Tests/Cartesian_Subgrid_Device.py
OPENEMS_SUBGRID_DEVICE=metal python3 python/Tests/Cartesian_Subgrid_Device.py
```

Qualification source parent is `039b128`; the exact source inventory, loaded
libraries, commands, failed development controls and numerical metrics are in
`experiments/openems-cartesian-local-subgrid-20261112/native-metal/`.
Loaded canonical candidate core SHA256:
`e3b143d182dda9ebb8d9aa1c85740476702207e04c12a3f75b54c8374f19e6bf`;
nf2ff: `8bf8d8e01be96032d94b71e3aef6e7c2645954354e7255763c575ed74438581c`.
Canonical qualification passes CTest 15/15, physical full-field controls 4/4,
Python unit modules 13/13, native device fixtures 4/4, the native CPU base and
six increment groups, two frozen CPU dipoles and three ordinary-engine groups.
The two frozen dipoles retain identical arrays against the explicitly merged
master reference. Source/runtime inventories, including five earlier frozen
prefixes, remain unchanged throughout this qualification.

Device ports, absorbers, lossy media/conductors, native device antenna/NF2FF,
structured topology scaling and multiple devices remain separate milestones.

### Scoped device resistor ports

The shared point operations now implement the CPU reference's small midpoint /
Woodbury resistor solve, sparse port response, and ordered soft sources. The
unsolved transpose curl is gathered only at source edges before mass solving;
source mass neighbours preserve the CPU work identity. Device loss/work uses
compensated float pairs, refreshed with the fields only at complete processing
boundaries. No CPU resistor updates or half-step field downloads are used.
Native geometry remains restricted to nonoverlapping positive R on independent
local edges; distributed/reactive/modal ports remain unsupported.

The default core port controls cover passive decay, a driven Norton pulse,
zero-DC soft kicks at a mass-coupled interface edge, and mixed coupling through
20k steps. Maximum double/device field error is 1.055e-4 (gate 2e-4); maximum
work/loss/balance error is 1.167e-6 (gate 2e-4). Reference backend and physical
Apple M5 Pro Metal controls agree. Native 24 GHz microstrip identity/outside/
inside/reversed cases retain exactly identical deterministic-colour input XML.
Across 20–28 GHz, max complex Z error is 1.501e-5 ohm and max complex S11/S21
error is 1.273e-7 (gates 0.01 ohm / 2e-4). Time U/I, phase, power and polarity
gates remain explicit. This is implementation precision on a deliberately
coarse control, not a calibrated manufacturing line or device antenna.

An inside-refinement fixture exposed Cholesky fill below normal FP32 range.
The device conversion now explicitly drops only subnormal off-diagonal fill,
retaining every normal entry and positive diagonal. Its summed perturbation
per block must be below `eps(float)^2 * smallestDiagonal`; unsafe/diagonal
underflow, nonfinite values and invalid dimensions fail. The observed fixture
has 136 tails, total L1 4.088561e-39 and maximum block-relative bound
1.023807e-39. Counts/bounds are logged in scientific notation. CPU double
factors remain unchanged; this rule does not silently truncate fields or curl
coefficients. Focused tests exercise preservation, bounds and rejection.

Isolated canonical qualification passes CTest 17/17, Python unit modules 13/13,
eight physical Metal core controls (fields and ports), four native field and
four native port geometries, CPU base plus six increment groups, two unchanged
merged-master dipoles and three ordinary-engine groups. Source/runtime
inventories, including six earlier frozen prefixes, remain unchanged. Parent
source is `4a7cf69`; contributor evidence:
`experiments/openems-cartesian-local-subgrid-20261112/device-ports/`.
Core SHA256: `bd64fa18e484a5393340a1c0473cf3bd0b180bffd17bca1b94afc420f3df1f49`;
nf2ff: `0ac606c88ee32284b6e8f9b9b44735c4d39155609b02b24651ea0fa858b7fad5`.

**Precision limitation of the preceding `0b255a8` increment:** a strong unmodulated Gaussian soft kick on a
mass-coupled edge misses the unchanged 2e-4 full-field target at 10k steps:
2.956e-3 with ports, identically on reference and Metal. It also fails without
ports on the immutable preceding `4a7cf69` runtime (2.774e-3). Changing only
the pulse to a zero-DC derivative reduces that control to 5.604e-5. This supports
a DC/FP32 diagnosis, not a proved cause or a correction, and is not labelled
an upstream regular-engine defect. The failure remains reproducible:

```sh
build/tests/test_cartesian_subgrid_device_ports --dc-stress
build/tests/test_cartesian_subgrid_device_ports --metal --dc-stress
```

Against the preserved preceding runtime these commands fail their field target;
they are not counted as passing controls or hidden by a higher tolerance.
The follow-up correction below retains these failures and adds a passing
regression. Default RF /
zero-DC passes do not establish arbitrary-source accuracy. Investigate curl
nullspace preservation, metric/factor representation and accumulation before
closing U05 or claiming full device parity. Device PML/Mur, losses/conductors,
antenna qualification and production scaling remain open.

### Device curl and timestep precision correction

The preceding strong-DC failure is retained against the immutable `0b255a8`
runtime. The follow-up correction keeps each normalized curl/transpose
coefficient as a high/low FP32 pair, uses explicit product-error FMA and
compensated dot sums, and preserves the supplied timestep/final field addition
through two-part FP32 arithmetic. No FP64 device type, source filtering, global
fast-math change or CPU half-step fallback is introduced. CPU double updates
and existing ordinary engine kernels remain unchanged. Non-normal coefficient
or timestep residuals reject explicitly; the separate bounded Cholesky-fill
conversion remains unchanged.

Arithmetic controls isolate coefficient/timestep rounding in this case; the
all-double diagnostic reproduces the core exactly. A coefficient-only candidate
improved DC but failed the unchanged 20k mixed-source test (2.59649e-4 > 2e-4).
Its source diff/runtime and failure remain preserved. The final candidate also
corrects timestep/final-addition precision; tolerances were not increased.
The test-only FP32 oracle explicitly models that arithmetic independently;
a separate double-core comparison remains mandatory.

| Scoped control | Final maximum normalized field difference |
|---|---:|
| Four physical Metal field/source controls, including 20k source-free | 8.97473e-06 |
| Four passive/Norton/mixed-source port controls, up to 20k | 6.48865e-06 |
| Previously failing 10k interface DC stress | 1.7602e-05 |
| Native field/probe/dump geometries | 2.7312397e-06 |

The DC stress uses the original 2e-4 field/work/balance gate on reference and
physical Metal, and is now a permanent CTest regression. Native complex Z error
is at most 1.0406748e-05 ohm, complex S11/S21 error
8.1456147e-08. CTest 18/18, Python unit modules 13/13,
physical fields/ports 4+4 and DC on both backends, native fields/ports 4+4, CPU base plus six groups,
two unchanged frozen CPU dipoles and three ordinary-engine groups pass.
Seven prior frozen runtime prefixes and the source/runtime inventory stay
unchanged. This does not establish arbitrary-source/long-time accuracy or
close the broader U05 antenna/backend matrix.

Three interleaved completed-batch Metal comparisons (24³ base / local 2×,
10,000 steps, final readback included) against `0b255a8`:

| Metric | Prior RF-port runtime | Corrected candidate |
|---|---:|---:|
| Median batch seconds | 0.714895 | 0.840039 |
| Accounted device bytes | 11764896 | 16580256 |

The corrected batch is about 17.5% slower and accounts for about 40.9% more
device memory on this fixture;
these small-core measurements do not predict antenna speed. The separately
compiled reproduction driver is [benchmark/subgrid_device_metal.cpp](benchmark/subgrid_device_metal.cpp). Added residuals
cost one FP32 value per curl/transpose nonzero. Expanded topology/coefficients
and host batch mirrors remain; structured interiors and scale guards are
independent open milestones. A smaller factorized-metric arithmetic experiment
is retained as unqualified research, not substituted for these controls.

Configure-time parent `0b255a8`; exact source/build/loaded-library inventories:
contributor archive `experiments/openems-cartesian-local-subgrid-20261112/device-curl-precision/qualification.json`.
Frozen prefix `output/openems-subgrid-device-curl-20261115/install`:
core SHA256 `35aea4bce694036437ec9c0a31c623a8eddeda2d07dc06eed3b9fc64dcc1f0e4`;
nf2ff `cc67dadd9b55c5a5d1f521d4910e8f5bc4a3cc98d60a95cbd630ec1f0b69b1c8`. Reproduction uses the existing device tests, plus
`ctest --test-dir build -R cartesian_subgrid_device --output-on-failure`.
Device Mur, loss/conductor, PCB antenna/scalability and physical HIP/CUDA
qualification remain open. The following increment adds scoped device PML/dipole evidence.

### Scoped device PML and dipole radiation

The existing shared backend now executes the paired classical convolution PML
on normalized compact fields. This is the same scoped PML described above,
not ordinary Yee UPML/CFS parity. Derivative histories, half-step corrections
and signed-work compensation stay on device. Pre/Post/Apply hooks preserve the
extension order; Post stages the correction and Apply publishes it. Complete
field mirrors/work diagnostics are read once per processing batch. No CPU
half-step fallback or per-half-step field download is introduced. Native Mur
still fails explicitly rather than silently using the CPU.

The implementation checks fixed timesteps, derivative bounds and coefficients,
unique grouped target ownership, complete-step observation and hook ordering.
Initial subnormal state tails may be zeroed only when their summed magnitude
relative to the original maximum state is below `epsilon(float)^2`; all normal
state values remain unchanged, an entirely unrepresentable state is rejected,
and counts/bounds are logged. This is separate from bounded Cholesky fill.
The double CPU reference and previously qualified runtimes are unchanged.

Qualification on Apple M5 Pro and the explicit reference backend:

- Eleven same-input PML controls cover six normals, graded/local packets,
  near-unit PEC reflection negatives and multi-region six-face corners.
  E/H samples include the PML and interior. The 100k-step packet and 10k-step
  corner controls pass on both backends. Largest sampled field/signed-work
  error across these core controls is 5.91e-7 (unchanged 2e-4 device gate).
- Six native normal/TE30/TE45/TM30/graded/local geometries yield 18 CPU/reference/
  Metal comparisons of E/H time records and FD planes. Largest normalized
  difference is 1.508e-6, below the 2e-4 gate.
- The public solid PEC vacuum dipole runs with identity and local-Z refinement
  on CPU/reference/Metal. Exact native contact, six-face surface completion,
  near fields, complex Etheta/Ephi on two surfaces, accepted/radiated power,
  gain/directivity and time/FD consistency are checked. Missing-face and broken
  contact controls discriminate failures. Device errors are at most 4.927e-6 ohm
  in Z and 4.914e-7 for the sampled complex near/far-field vectors. All retained
  double CPU arrays are exactly unchanged from the frozen merged-master fixture.
- CTest19/19, Python unit modules 14/14, physical field/port/DC controls,
  native fields/ports, CPU base plus six increment groups, two frozen dipoles
  and three ordinary-engine groups pass. Source/runtime identities and eight
  prior frozen prefixes remain unchanged throughout qualification.

Reproduction (after building/installing the isolated candidate):

```sh
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_subgrid_device_pml --metal --long --long-corners
OPENEMS_SUBGRID_DEVICE=metal python3 python/Tests/Cartesian_Subgrid_Device_PML.py
OPENEMS_SUBGRID_DEVICE=metal python3 python/Tests/Cartesian_Subgrid_Device_Radiation.py
```

Configure-time parent `ca6f755`; contributor evidence is under
`experiments/openems-cartesian-local-subgrid-20261112/device-pml/qualification.json`.
Frozen prefix `output/openems-subgrid-device-pml-qualified-20261115/install`:
core SHA256 `3bfa5d5a9311e44e1572577dbb4b3cf70c08a8fa05148f2cc35eea42ab9f23fd`;
nf2ff `552d8d5bb139cc69c1b3428072c2f956d1d10e708d07b792d9d61bb12ec04b78`.

This is a small vacuum PEC antenna/device milestone, not representative PCB,
lossy conductor, general CAD/source/interface, structured scalability or
multi-device qualification. DFT/NF2FF remains host processing. Current guards,
expanded upload metadata and host topology mirrors remain. CUDA preparation
uses the **existing HIP backend**, shared operations and NVIDIA `nvcc`, as
[documented separately](CARTESIAN_SUBGRID_CUDA.md). Host harness/selection guards
pass, including rejection of a non-CUDA build; no NVIDIA compilation or
physical HIP/CUDA evidence is claimed here.

### Exact device curl patterns

The device execution now shares exactly equal translated row patterns for curl
and its gather transpose. Each output retains a pattern number and input-base
index; coefficients remain the same FP32 high/low values, in the same term
order. No coefficient quantization, magnetic-history merging, interface-solve
change, field precision change or global fast-math setting is introduced.
The shared operation definition still serves the existing reference/Metal/HIP
backends. Physical qualification here is Metal only.

The independent oracle retains its expanded CSR export. Reconstruction and
input guards check every index/coefficient bit, empty/translated/irregular rows,
signed-zero residuals and malformed input. The new shared gather is compared
with the same expanded gather for overwrite and additive updates. All prior
FP32, double field, port, PML and dipole gates are unchanged. CTest20/20, Python
unit modules 14/14 and the full previous qualification matrix pass. In addition,
11,255 NPZ arrays, HDF5 datasets and time-probe arrays are numerically identical
to the frozen `bb2edbb` expanded device runtime, including CPU/reference/Metal
fields, ports, PML and dipole radiation. Nine prior prefixes remain unchanged.

Three interleaved physical Metal completed-batch pairs on the same local24
fixture (400 warm-up +10,000 measured steps, final readback included):

| Metric | Expanded upload/runtime | Exact patterns |
|---|---:|---:|
| Median batch seconds | 0.846029 | 0.570051 |
| Accounted device bytes | 16,580,256 | 3,781,236 |
| Process peak RSS bytes | 330,563,584 | 318,619,648 |
| Device/backend setup seconds | 0.040093 | 0.042943 |

This fixture takes 32.6% less batch time and 77.2% less accounted device memory.
Process peak RSS falls only 3.6%; it includes macOS/Metal overhead. Device setup
includes backend/cache initialization and is not a cold compiler or production
startup benchmark. All three paired final energies are exactly equal. These
measurements do not predict large-antenna or NVIDIA speed. The reproduction
driver is [benchmark/subgrid_device_metal.cpp](benchmark/subgrid_device_metal.cpp).
The public local-Z dipole separately accounts for 167,703,460 ->90,212,052
allocated device bytes, including unchanged PML histories; no native timing
speedup is claimed from the qualification runs.

Both expanded CSR tables still exist temporarily in the **host upload** before
compression, alongside the pattern dictionary. This is not a general reduction
in peak setup RAM, a structured field kernel or a reason to raise guards. Host
topology/observation mirrors, large-interface scaling and device DFT remain open.

Configure-time parent `bb2edbb`; actual source/runtime inventories, interleaved
logs and exact prior-array comparison are in contributor archive
`experiments/openems-cartesian-local-subgrid-20261112/device-patterns/`.
Frozen candidate prefix `output/openems-subgrid-device-patterns-20261115/install`:
core SHA256 `45d0681845f90409d0eb6bba13b2a81e1bbb6fdd4cd49a7173835ff7f2766f33`;
nf2ff `b4216540dc74cf21186c8042f662ebcff8445ef0c21f31b87b963971bb8d4f83`.
Reproduce the existing physical device suite plus
`ctest --test-dir build -R cartesian_subgrid_device_patterns --output-on-failure`.
No project PCB/lossy-conductor, physical CUDA or multi-device qualification is
added by this metadata representation.

### Structured independent-edge address index

The compact core now addresses independent edges through structured base/local
slot arrays and a base-cell region-owner table instead of a tree. A base-aligned
edge keeps one common address when seen from either side of a local boundary;
other edges resolve through at most four adjacent owners. Boundary/PEC/dependent
slots remain explicitly absent. Original traversal, independent-state numbering,
coordinates, constraints, mass/curl coefficients and arithmetic order are unchanged.
Positive observation lookup is const and allocates nothing. The explicit
`uniformFine` numerical control uses one internal full-domain address mask,
not only the user's local boxes.

Eight tree-oracle address fixtures cover all leaf edges, duplicate/sparse/missing
keys, forward/reverse lookup, independent ratios/regions, ratio 16, 64 regions
and global uniformFine 222/412. The standalone core and address test also pass
ASan/UBSan. CTest 21/21, Python unit modules 14/14 and the previous full native
CPU/reference/physical Metal port/PML/dipole/ordinary-engine matrix pass.
All 11,255 sampled numerical arrays remain exactly equal to frozen `f5aedf2`,
including time probes and complex port/radiation values. Ten prior runtime
prefixes remain unchanged.

Three interleaved standalone CPU pairs per fixture, with each executable built
from its matching source/header snapshot:

| Case | Setup before/after [s] | Accounted before/after [MiB] | Peak RSS before/after [MiB] |
|---|---:|---:|---:|
| Uniform16 |0.0118/0.0058|2.45/1.90|4.97/4.36|
| Local24 |0.0888/0.0433|17.05/13.31|20.59/18.14|
| Strip512 |0.3097/0.1635|51.18/39.45|62.39/50.14|

Across all seven fixtures, setup is about 43–51% shorter; update time is
essentially unchanged and sampled fields/CFL are exactly equal. Buffer
accounting uses actual structured capacities, while the old tree accounting
was approximate; process RSS is recorded separately. These are standalone-core
measurements, not ordinary-engine, native-antenna or GPU speed claims.
Reproduce with [benchmark/subgrid_compare.py](benchmark/subgrid_compare.py).
Changing a private Grid member changes its C++ object layout: rebuild C++
callers with matching headers, and never compare shared libraries by relinking
an old object compiled for the other layout.

Configure-time parent `f5aedf2`; exact source/runtime inventories, sanitizers,
interleaved benchmark and old-array comparison are in contributor archive
`experiments/openems-cartesian-local-subgrid-20261112/edge-index/`.
Frozen prefix `output/openems-subgrid-edge-index-20261115/install`:
core SHA256 `b73b6824b641ddffcd805acd32684bce0393ab9df33dea7009d86fefd0bcdc5b`;
nf2ff `63c8e91a7503371485496d9f01c5d9da4dfc6e538d4c0fd35a2f7cfe00a44908`.
The new [address regression](tests/test_cartesian_subgrid_edge_index.cpp) runs
with the existing suite. This is an address-storage increment, not structured
interior field kernels, larger allocation guards, general CAD/losses, physical
CUDA or multi-device qualification.

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
build/tests/test_cartesian_subgrid_mur --long
build/tests/test_cartesian_subgrid_pml --long
build/tests/test_cartesian_subgrid_pml --long-corners
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

The prototype has explicit allocation guards: at most 16,384 base intervals per
axis, one million base cells, one million leaves summed across all regions,
64 regions and 16,384 coupled edges per interface factor. Above 1,024 edges,
half-bandwidth must be at most four and the one-RHS inverse bound must succeed;
see [bounded larger interfaces](#bounded-larger-interface-factors). These are
implementation safeguards, not physical limits of subgridding. Exact banded
storage reduces current interface cost; increasing the guards still requires
separate memory/setup/stability qualification rather than removing assertions.

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

Beyond the scoped vacuum PEC dipole: broader material/PEC geometry and
interfaces, distributed/modal ports/sources, absorber spectra, NF2FF placements
and representative antenna spatial qualification remain necessary.
Before practical acceleration: structured interior kernels, compact interface
metadata, memory/RSS and interleaved native benchmarks. Then device kernels,
domain partitioning and communication can reuse this numerical contract.
Nested refinement levels and local subcycling are separate changes. No
running or frozen production model should adopt this prototype silently.


## Bounded larger interface factors

Configure-time parent `3a502cc`; no mass, factorization, curl, state ordering or
update arithmetic changes. The old per-axis and per-component caps blocked long
narrow regions despite their already banded factors. Base intervals and coupled
interface edges may now reach 16,384. Components above 1,024 edges additionally
require half-bandwidth<=4 and a verified positive-bipartite one-RHS inverse
bound. Reject broad factors before allocation and reject an unavailable large
one-RHS bound instead of falling back to a quadratic column sweep. The existing
one-million base/leaf caps and 64-region/non-nesting constraints remain.

[CPU/resource controls](tests/test_cartesian_subgrid_factor.cpp) preserve the
frozen dense oracle, exercise the old rejection boundary, test guard boundaries
and a 4,096-edge component over 2,000 steps. The latter invariant drift is
7.33e-15; factor/work storage stays band/block-sized. [Device controls](tests/test_cartesian_subgrid_device_large.cpp)
exercise 1,536-edge components in XYZ orientations and a driven case, retaining
processing reads without restarting fields. Physical Metal maximum error against
double is 1.08e-6 (gate 2e-4); independent FP32 oracle arrays are identical.

Additional standalone maximum-bound controls use base cells `{8194,4,4}`,
spacing `{.001,.0013,.0009}` metres and region `{{2,1,1},{8192,3,3}}` with
ratio 2/2/2. The actual 16,384-edge component (half-bandwidth 2, 360,424 leaves)
retains its invariant over 2,000 CPU steps within 1.0e-14. Its 393,168 factor
entries replace 1,074,265,984 dense-equivalent entries. Physical Metal 512-step
passive/driven controls have double error 1.56e-6 and identical FP32 oracle
arrays. Raising only the long dimension to 8195 produces an actual 16,386-edge
component, rejected before factor allocation. These are functional/resource
controls, not RF convergence or production throughput measurements.

The immutable full suite passes all 18 stages: CTest 22/22, Python 14/14, the
previous native CPU/reference/Metal, ports, DC, PML, dipole, polygon and ordinary
engine controls. Twelve older core/NF2FF prefixes remain unchanged; 11,255 prior
numerical arrays exactly equal the frozen polygon-qualified runtime. Separate
maximum-bound sources, hashes and logs are retained with
`experiments/openems-cartesian-local-subgrid-20261112/large-interface/` in the
project evidence tree; `qualification.json` distinguishes those controls from
the full suite. No new ASan or physical CUDA claim is made.

Frozen candidate prefix:
`output/openems-subgrid-large-qualified-20261115/install`.

- Core SHA256: `3eea599befbe5c3bd590e6a50e8e2f0cae8dd6ee44d9ddf669adff31da1b59cd`.
- NF2FF SHA256: `bc6961138b4d8dba3689c773567d9bf82d492a0fb5742c819a0715cee64c046e`.

General broad-interface scaling, structured interior fields, proportional
setup memory, losses/conductors, representative PCB spatial accuracy and
multi-device execution remain open. Allocation success and energy conservation
do not establish antenna accuracy; this increment changes no active or
production runtime.

## Compact domain implementation

The subsequent [C++ domain increment](CARTESIAN_SUBGRID_DOMAINS.md) adds fixed
ownership, worker/coordinator contracts and boundary-only host transport, with
reference, physical Metal and single-RTX3080 logical-domain controls. It retains
the existing one-device native paths. Physical multi-GPU, domain ports/PML/loss
histories and distributed processing remain open; see the canonical checklist.


## Paired compact advance extensions

The experimental compact engine can use one explicitly paired
`Operator_Ext_CartesianSubgrid_Advance` /
`Engine_Ext_CartesianSubgrid_Advance` to replace its actual electric and
magnetic half-step. This is an extension boundary, not conductivity support
or post-step damping. Existing pre/post/apply ordering remains around the
CPU recurrence. Device providers must return an actual implementation through
the selected backend; an unsupported device cannot silently fall back to CPU.
Duplicate or missing providers reject. A failed recurrence or its diagnostics
poisons further continuation and field observation on that engine.

The transparent test recurrence in
[the C++ regression](tests/test_cartesian_subgrid_device_advance.cpp) checks
call counts, CPU hook order, complete-step observations, device failures and
unchanged paired PML fields. It rebuilds a nonzero excitation after operator
setup and separately tests zero initial fields; an earlier local test revision
had only initial-field coverage because setup reset its excitation. That gap
and the stronger supplement are retained in the evidence archive.

The 21-stage compatibility suite passes on CPU/reference/physical Metal:
CTest25/25, Python14/14, native fields/ports/PML/dipole/polygons, six CPU
increment groups, frozen dipoles and three ordinary-engine groups. Sixteen
older runtime prefixes and11,255 numerical arrays remain unchanged. The
stronger canonical source-only supplement has maximum Metal/double field
error6.164e-7, below the unchanged2e-4 gate. Transparent PML fields remain
exact. A separate physical RTX4090/RTX3080 suite retains both domain orders,
ports and single-device controls; its full CTest25/25 passes. The same stronger
source-only supplement also passes on each physical NVIDIA device. These are
bounded extension/field controls, not lossy/native PCB or scaling acceptance.

Frozen local core SHA256:
`75e8b7c09f7645132623d0753010f80c5f1d4c6cc7c2101514f785fc763d1658`.
Frozen CUDA core SHA256:
`1c7e917010b6c57d1f3247a66bd6c05a8cdc2ba3069162fab855899af1fa3bfe`.
Evidence is under the contributor's
`experiments/openems-cartesian-local-subgrid-20261112/advance-provider/`;
source-only supplements use unchanged cores and separately identified tests.
No production or antenna runtime is changed. Conductivity, absorber/material
history combinations, distributed processing and representative antenna
qualification remain separate milestones.


## Native source waveform tails

The experimental native adapter now follows the ordinary source's sample-zero
fallback outside a stored nonperiodic waveform. The two-sample Step retains
one, rather than silently dropping to zero after step two. Gaussian/Dirac
zero tails are unchanged. C++ periodic sample indices use the global timestep;
normal native periodic/steady-state scripting remains explicitly rejected.

A frozen negative control reproduced the old Step error (about 0.9074 maximum
field discrepancy). Direct Step/custom/Dirac/periodic controls retain the
unchanged 2e-4 CPU/device field gate and a negative legacy-tail control. Actual
Python/native Step uses the real two stored unit samples; CPU/reference/Metal
and HIP/NVIDIA comparisons pass. This is not H-source, delayed-source, TFSF or
steady-state parity.

Local source/configure parent is b48e819, plus this focused source-tail patch:
core SHA-256 237a52de2af1d0d1be7173fe83239e659c9804d28fc9b8d0786bb570fb849709.
All 21 compatibility stages pass, including26/26 CTest, 14/14 Python,
fields/ports/PML/dipole/polygons and two frozen CPU dipoles. Seventeen older
runtime prefixes and 11,255 numerical arrays remain unchanged. Earlier b48e819
ordinary-engine controls are retained with unchanged ordinary implementation.

The separately built HIP/NVIDIA runtime passes 26/26 CTest, eight waveform cases
on each physical RTX4090 and RTX3080, actual native Step, and retained two-device
field/port controls in both orders. Largest waveform field error 2.454e-6.
Runtime/source/terminal hashes and compact logs are in the contributor's
`native-waveforms/{qualification,cuda-qualification,prior-array-comparison}.json`
archive. The direct periodic tests do not override the scripting restriction.

```sh
build/tests/test_cartesian_subgrid_device_waveforms --metal
# Use a separately built HIP/NVIDIA runtime:
build/tests/test_cartesian_subgrid_device_waveforms --cuda
```
