# Scoped passive surface states for Cartesian subgrids

Generated with Codex at the contributor's request. This is an experimental
C++ component in PR29, not native conducting-sheet/skin/CAD admission.

## Numerical contract

`Operator_Ext_CartesianSubgrid_Conductivity::SetSurfaceAdmittances` accepts
physical oriented voltage paths with a finite nonnegative conductance G and
positive R/L branches. The component uses the same paired advance provider as
constant bulk conductivity. Independent native providers remain forbidden.

For normalized electric unknown q, magnetic unknown h and x = sqrt(L) I,
the component integrates

```
M dq/dt = C^T h - G q - K x + source/port terms
    dx/dt = K^T q - rho x,             rho = R/L.
```

Eliminating the branch midpoint gives D = 1 + dt rho/2 and
G_eff = G + dt K D^-1 K^T/2. Use M + dt G_eff/2 for the electric solve and
its coupled resistor response; update x from that same electric midpoint
before ordered soft sources. The CPU and shared device operations use the
same equations. Zero-pole paths retain the existing bulk-only path.

Stored branch energy is x^T x/2. Ohmic loss uses the physical G and rho at
the midpoint; G_eff is an elimination coefficient, not extra dissipation.
Source work, resistor loss, bulk loss, PML signed work and surface loss retain
separate ledgers. Branch energy contributes to the compact engine's energy
observation. Midpoint passivity does not establish temporal accuracy: very
stiff modes can ring slowly because this scheme is not L-stable.

Histories stay on the selected device. Preparation builds edge/branch incidence
once; updates require no atomics or CPU branch-history download. Reductions
occur at complete-step observations. Matching rebuilt C++ callers and backend
plugins are required by the added virtual surface diagnostics; plugin mode
itself is not covered by the built-in backend tests. Dynamic CPU/device history
migration, restart and repeated attachment of warm state are not qualified.

Surface/PML combination retains the prepared disjoint-support contract in
[CARTESIAN_SUBGRID_CONDUCTIVITY.md](CARTESIAN_SUBGRID_CONDUCTIVITY.md#disjoint-conductivity-and-pml).
Overlap, Mur and incomplete observations reject. This does not implement
surface histories across multiple devices.

## Boundary-roundoff correction

A reproduced PR29-only reconstruction defect selected the positive-side leaf
within the existing FindLeaf roundoff envelope but retained a tiny negative
basis weight from the coordinate just below that boundary. A physical local
surface path acquired an unrelated root at relative weight about -1.39e-14
and correctly failed the mass-graph guard. The correction applies that same
bounded envelope to local basis fractions and integration-cut support. It
does not move mesh lines, geometry or refinement corners and retains genuine
interior interpolation. The mass-graph guard is unchanged.

The negative input/source/binary/log remains in the contributor's evidence
archive. `test_cartesian_subgrid_integral_roundoff` covers all three normals,
translations, two-ULP ties, resolved interior paths and refined magnetic cuts.
This is not an established defect in the ordinary upstream solver.

## Qualification and reproduction

Build with `OPENEMS_BUILD_TESTS=ON`, matching dependencies and the desired
existing Metal or HIP backend. No new backend or dependency is introduced.

```sh
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_subgrid_surface --reference
build/tests/test_cartesian_subgrid_surface_guards
build/tests/test_cartesian_subgrid_conductivity_pml --surface-only --metal
# NVIDIA HIP build: use --cuda in place of --metal.
```

The surface executable accepts the device selectors supported by
`tests/subgrid_device_target.h`. The independent full monolithic double oracle
compares fields and physical energy/loss against the eliminated recurrence.
Seven cases exercise constrained traces, a pure inductor, bulk E/H losses,
coupled ports, ordered soft sources, stiff relaxation and seven poles from the
unchanged public smooth-Cu surface fitter (Rq=0). Two disjoint PML cases add a
real source/port packet and a six-face corner. Fourteen negative cases retain
path/scalar/rate/mass-graph guards. No private supplier data is in these inputs.

All 30 local stages passed: CTest 33/33, Python 14/14, physical Metal, native
fields/ports/PML/dipole/polygon controls and prior CPU dipoles. Maximum sampled
surface Metal field error is 4.408e-6 and normalized energy balance error is
2.078e-7, with the original 2e-4 gates. Ten-thousand-step pure-L and fitted-state
controls also passed CPU/reference/Metal (maximum field error 5.624e-6,
balance 7.820e-7). All 11,255 compared prior numerical arrays remain bit-identical.

All 19 physical CUDA gates passed, including CTest 33/33. Surface histories and
joint PML run separately on RTX4090 and RTX3080. Retained **lossless** domain and
port controls run on the actual pair in both orders; this is not distributed
surface/bulk/PML qualification. The CUDA surface maxima are 4.40722e-06 for fields
and 2.07763e-07 for normalized balance, with unchanged 2e-4 gates.

Pinned numerical parent: `08e22748daf3d5b775e3e58b6e2786adcfdbfab3`.
Local core SHA256: `3644721f4789b13006eba727f539dd6f794bc38bfbfb1669b4d8c30512619720`.
CUDA core SHA256: `4b9a4afc469c30518e6a5b49a593fa81af4dceb85510805123d2e29c46420a03`.
Exact numerical files were matched independently to both immutable qualification
source inventories before transfer into PR29. The contributor's local archive
`experiments/openems-cartesian-local-subgrid-20261112/lossy-materials/surface-device/`
contains receipts, source/runtime hashes, command logs and comparison arrays;
`output/openems-subgrid-surface-device-20261116/FROZEN-r01.json` protects the
qualified prefix. These local archives are not public PR attachments.

## Remaining scope

A surface state attached to a voltage path does not automatically represent
physical metal. Native volume exclusion, exposed-face quadrature, reference
planes, finite-thickness conventions, conductor contacts and short-line phase/
attenuation need separate implementation and comparison. General ADE media,
roughness, overlapping absorbers, arbitrary CAD/ports, distributed histories,
DFT/pause, large-interface scaling and the PCB array remain open. Frozen
production/antenna runtimes are unchanged.
