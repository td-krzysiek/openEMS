# Experimental planar Cartesian multi-grid

This prototype adapts the existing cylindrical multi-grid operator/half-step
scheduler to a **single planar Cartesian interface**. It is a numerical
development step, not an antenna-ready embedded refinement box.

## Supported domain

With `CartesianMultiGridX = s`, the grid at `x < s` uses every second Y line;
the grid at `x >= s` retains every Y line. X and Z lines and the timestep are
shared. Thus the refinement ratio is 2:1 **in Y only**, not in all three axes.

```text
                   exact X mesh line s
  coarse Y spacing       |       fine Y spacing
  +-------+-------+      |      +---+---+---+---+
  |       |       |      |      +---+---+---+---+
  +-------+-------+      |      +---+---+---+---+
  |       |       |      |      +---+---+---+---+
  +-------+-------+      |      +---+---+---+---+
             same X/Z coordinates; one shared timestep
```

The operator requires:

- Cartesian coordinates and the basic, SSE or CPU host engine.
- Uniform Y spacing, with an odd number of Y lines (at least nine).
- An exact interior X split line, with three collar cells on either side.
- PEC on all six outer faces.
- No material, metal, source, probe or dump primitive intersecting that collar.
  The conservative bounding-box test includes primitives outside the domain;
  it does not infer clearance from a few Yee samples.
- Only the excitation extension. TFSF, dispersive/sheet/lumped extensions,
  absorbing boundaries and GPU engines are rejected, without fallback.
- Untransformed primitives with verifiable Cartesian bounding boxes. Other
  property types (including base-operator lumped RC elements) are rejected.

The collar's half-width is three local X cells, **not** three Y cells. Do not
place a port or an NF2FF box across it. These diagnostics and material
interfaces need separate qualification before antenna use. Homogeneous
background and structures wholly away from the collar can be represented,
but the numerical regression here qualifies vacuum cavities only.

## Coupling and timestep

The implementation reuses `Engine::SetSubGrid()` and its half-step barriers:
update both E fields, restrict the interface voltages, update both H fields,
then transfer the interface currents. No engine scheduler change is needed.

Voltage and current unknowns are edge integrals. At the split, the restriction
`Vc = R Vf` samples Ez at even fine Y nodes and sums each pair of Ey edges:

```text
Vc,z[j] = Vf,z[2j]
Vc,y[j] = Vf,y[2j] + Vf,y[2j+1]
```

The conjugate current transfer is `If = transpose(R) Ic`, with component signs
supplied by the existing curl. This makes the signed interface power pairing
`Vf^T If = (R Vf)^T Ic`. Point-interpolating the currents instead is not this
adjoint: an initial implementation developed a late-time instability and was
discarded. This observation concerns the new Cartesian port, not a finding of
a defect in cylindrical multi-grid.

The Ez interface capacitance uses the fine-side half dual cell and the
coarse-side half dual cell at even Y nodes. Relative to the native full fine
capacitance its factor is
`(right + (even_y ? 2*left : 0)) / (left + right)`, where left/right are half
the adjacent X intervals. Ey uses row-sum lumping of the restricted coarse
capacitance. This is a low-order, diagonal-mass interface, not a higher-order
hanging-node solve.

Both operators first calculate their native timestep bounds. The common step
is also bounded by the coupled normalized curl A: each incidence has magnitude
`sqrt(VI_e * IV_h)`, including dt; the code enforces
`sqrt(norm1(A) * normInf(A)) < 2`. It rebuilds with a smaller step when needed.
An explicitly forced step that violates either bound is rejected. This extra
bound is conservative; a subgrid does not necessarily allow a larger timestep.

Inactive parent fields are linearly reconstructed for processing at batch
boundaries. The current ghost plane used for interface flux is deliberately
not overwritten with point-interpolated H. It is not a physical H sample,
which is one reason observations are prohibited near the interface.

## Usage

Python:

```python
from openEMS import openEMS
fdtd = openEMS(CartesianMultiGridX=0.0, TimeStepFactor=0.8)
fdtd.SetBoundaryCond(['PEC'] * 6)
# Configure CSXCAD, an exact split X line and uniform Y lines.
# Keep sources/observations/structures outside the interface collar.
fdtd.Run('simulation', engine='basic', numThreads=2)
```

Matlab/Octave uses the existing generic XML-attribute path:

```matlab
FDTD = InitFDTD('CartesianMultiGridX', 0, 'TimeStepFactor', 0.8);
FDTD = SetBoundaryCond(FDTD, repmat({'PEC'}, 1, 6));
```

The optional attribute is on the FDTD element:
`<FDTD CartesianMultiGridX="0" ...>`. Existing inputs without it retain their
normal operator. Coordinates are in the mesh's drawing units. See
`python/Tests/Cartesian_MultiGrid.py` for a complete small source/probe model.

## Reproducing the tests

Build with `OPENEMS_BUILD_TESTS=ON` and the dependencies described in README:

```sh
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_multigrid --long
python3 python/Tests/run_testsuite.py --unittests Cartesian
python3 python/Tests/run_testsuite.py --engine=cpu Cylindrical_Cavity Engine_Compare
```

The C++ test checks the transfer's power adjoint, basic/SSE/CPU E and H agreement
within 2e-4 of each field's peak, random three-component field stability, a
30,000-step closed cavity and a transmitted pulse against a uniform fine grid
at the same timestep. `--long` extends the mode test to one million steps.
The mode must be within 2% of its analytic frequency and 1% of the uniform fine
reference. The sampled field norm must remain bounded; it is **not** an exact
discrete Hamiltonian or an antenna energy-decay qualification.

For the fixed physical Gaussian-pulse test, halving the mesh twice reduced the
Ez peak-normalized time-record discrepancy from 6.53% to 3.04% to 1.43%.
The finest result must be below 2% and the error must decrease. The coarse
result is retained as evidence of limited accuracy, not counted as a 2% pass.
The orthogonal, constant-Y Ey pulse is a transparent-interface control.
This does not establish arbitrary-angle scattering accuracy.

## Remaining work

The parent still allocates full fine fields and coefficients for processing,
and a coarse child is added. **No memory saving or speedup is claimed.** The
displayed cell count is active update work, not total allocated storage. Setup
currently builds coefficients multiple times and temporarily allocates row
sums to qualify the coupled timestep.

Before a useful antenna solver: compact the inactive storage; extend the
interface accuracy and its angular tests; qualify materials, PEC crossings,
absorbers and diagnostics; implement bounded local boxes and additional
refinement axes. GPU kernels and multi-device transport come after those
numerical contracts. Local subcycling is a separate numerical change.
