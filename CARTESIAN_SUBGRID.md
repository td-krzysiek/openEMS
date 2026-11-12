# Experimental compact Cartesian XYZ subgridding

`CartesianSubgridBox` replaces selected base cells with `rx * ry * rz`
children. `CartesianSubgridRatio` independently selects the refinement on
X, Y and Z: each factor is 1, 2, 4, 8 or 16. The default `[2,2,2]` retains
the original eight-child behaviour; `[1,2,8]` makes 16 children while leaving
X unchanged, and `[2,4,8]` makes 64 children. Spacing on each axis is divided
by its selected factor inside one aligned box. Outside the box
the coarse cells remain. This is a separate Cartesian numerical core, not
angular interpolation or a slab copied from cylindrical multi-grid.

It is a **serial CPU reference implementation**, suitable for numerical
review and small validation models. It is not yet a replacement for the
production antenna engine.

## Geometry, storage and execution

- One refinement level with configurable axis ratios and one rectangular box, strictly inside the domain.
  Each box face coincides with a base-grid line; at least one base cell remains
  between the box and each outer boundary.
- A Cartesian base mesh, uniform on each axis; the three axis spacings may
  differ. Arbitrary nonuniform lines, multiple/nested blocks and moving
  refinement are not implemented.
- Coarse cells covered by the box are removed. Only leaf cells and independent
  electric edges are stored. There is **no global fine parent E/H lattice**.
- One common timestep bounded by the coupled operator, not local subcycling.
  The smallest scale still limits the timestep.
- The reference uses double fields, sparse curl rows and interface mass
  factors. This favours inspectability, not production speed. Regular
  interiors need native structured kernels before a fair performance claim.

For base count N, M selected base cells and `R=rx*ry*rz`, the leaf count is
`N + (R-1) M`. A fully fine grid has `R N` cells. In the test with 8x8x8 base cells and a
4x4x4 selected block, the counts are 960 versus 4096. Accounted persistent
storage is about 1.09 versus 4.36 MB **for this same reference core**. These
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
# Supply a uniform Cartesian CSX grid; e.g. -8:2:8 on every axis.
# Add a zero-delay soft E box source and E/H probes or dumps.
fdtd.Run('simulation', engine='basic', numThreads=1)
```

Matlab/Octave uses its generic attribute serialization:

```matlab
FDTD = InitFDTD('CartesianSubgridBox', '-4,-4,-4,4,4,4', ...
                'CartesianSubgridRatio', '1,2,8', 'TimeStepFactor', 0.8);
FDTD = SetBoundaryCond(FDTD, repmat({'PEC'}, 1, 6));
```

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

Currently supported: homogeneous positive isotropic, lossless background;
outer PEC; zero-delay nonperiodic soft E box sources; E/H observations.
Unsupported requests fail at setup: internal material/metal objects, lumped
ports, conducting sheets, dispersion/loss, PML/Mur, periodic sources, transformed
primitives, other field types and GPU/SSE/new CPU engine selection. Only
`engine=basic` selects this reference adapter; its execution is serial.

## Tests and observed accuracy

```sh
cmake -S . -B build -DOPENEMS_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_subgrid --long
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
axis, 20 million base cells, one million leaves and 1024 coupled edges per
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

Before antenna use: material and PEC geometry mapping, nonuniform base lines,
absorbing boundaries, physical port/source qualification and scattering tests.
Before practical acceleration: structured interior kernels, compact interface
metadata, memory/RSS and interleaved native benchmarks. Then device kernels,
domain partitioning and communication can reuse this numerical contract.
Multiple refinement levels and local subcycling are separate changes. No
running or frozen production model should adopt this prototype silently.
