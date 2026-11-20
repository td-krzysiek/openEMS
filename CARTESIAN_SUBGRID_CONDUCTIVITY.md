# Experimental Cartesian-subgrid conductivity

**Qualified scope:** constant native isotropic conductivity boxes and the
paired single-device recurrence on CPU/reference/Metal and HIP/NVIDIA. The
combined compatibility gates below pass. This is not distributed-loss,
conductor/skin or production antenna qualification. A separately qualified
disjoint-support conductivity/PML combination is described below. Retain the
allocation and geometry guards from `CARTESIAN_SUBGRID.md`.

The later [cold C++ domain increment](CARTESIAN_SUBGRID_DOMAINS.md#prepared-conductivity-domains)
has a separate implementation/qualification record; the original single-device
evidence below does not qualify that path or native distributed processing.

## Numerical contract

Use constant isotropic electric conductivity in S/m and magnetic conductivity
in ohm/m. Native boxes use the existing CSXCAD `kappa` and `sigma` properties;
permittivity/permeability, priorities, drawing units and local-cell alignment
retain their existing meaning. Negative/nonfinite coefficients reject. This is
not a broadband loss-tangent, surface-impedance, skin-depth or roughness model.
Background conductivity, spatially weighted loss and unsupported material
primitives remain rejected. Constant diagonal loss is implemented by the later
R80 candidate below, with its own qualification scope; earlier isotropic test
evidence does not qualify that increment automatically.

The paired `Operator_Ext_CartesianSubgrid_Conductivity` builds the electric loss
form `G` from local physical-volume quadrature and the same restricted traces as
the original electric mass `M`. Its coupled implicit factor is
`A = M + dt G/2`; the electric increment obeys
`A (E_new - E_old) = dt (C^T H - G E_old)`.
Resistor-port responses and the small coupled port solve use that same factor,
instead of damping the result of an otherwise lossless timestep.

Soft-E sources remain ordered kicks after the implicit/port update. Their work
uses the original mass and transpose curl. Electric dissipation is the positive
midpoint quadratic form before those kicks. Magnetic conductivity uses its own
constitutive midpoint update and positive loss ledger. The current CFL bound
is retained; midpoint passivity does not imply L-stability or adequate spatial
resolution for a highly conducting solid.

One shared implementation uses the existing reference, Metal and HIP backends.
Device factors preserve the original trace/summation order, including full XYZ
refinement. Histories remain on device. Compensated per-trace/per-face ledgers
are reduced hierarchically only when an observation requests them. Allocation,
phase, fixed-timestep and incomplete-observation guards remain explicit.

## Evidence and reproduction

The independent dense oracle in
[`test_cartesian_subgrid_conductivity.cpp`](tests/test_cartesian_subgrid_conductivity.cpp)
constructs physical edge/volume forms rather than reading prepared private
factors. Twelve bounded cases cover zero/stiff loss, pure electric/magnetic
loss, full XYZ refinement, noncommuting mass/conductance, real pulse/DC sources,
source-only operation, coupled ports and repeated observations.
[`test_cartesian_subgrid_conductivity_guards.cpp`](tests/test_cartesian_subgrid_conductivity_guards.cpp)
checks invalid coefficients, construction, phases, rates and unsupported absorber combinations.

The earlier isolated XYZ prototype passed 27/27 CTest checks and twelve physical cases
on each of Metal, RTX4090 and RTX3080. The largest sampled FP32 field error
against the independent recurrence was `1.733e-6`; the largest normalized
energy-balance discrepancy was `2.329e-7`. The CUDA attempt also retained the
separate lossless domain/port controls on the physical two-card pair, in both
orders. **That does not run conductivity on multiple devices.**

Isolated XYZ library SHA-256:

- Metal: `7f2e8809dc0441dbb9f952cd971c9754c0e44c4463c87dff9c7ece1bae5739a0`.
- CUDA: recorded in the contributor's identity-checked `xyz-cuda-r01.json`.

Native controls in
[`test_cartesian_subgrid_native_conductivity.cpp`](tests/test_cartesian_subgrid_native_conductivity.cpp)
check analytic uniform-electric-field decay, second-order time convergence,
translated drawing units, priority/insertion-order changes, a higher-priority
lossless pocket and repeated preparation. This is a constitutive zero-curl
control, not a radiating cavity. Actual scripting controls in
[`Cartesian_Subgrid_Conductivity.py`](python/Tests/Cartesian_Subgrid_Conductivity.py)
compare native CPU/reference/device fields and require a measurable material
effect. Both native groups pass in the combined CPU/reference/Metal qualification and
on each physical RTX4090/RTX3080. Native Python electric/magnetic/combined
material controls pass on Metal and CUDA; their largest sampled CPU/device
probe discrepancy is 1.425e-6.

A separate native homogeneous PEC-cube TE110 pilot now adds continuum wave
evidence. At 6/12/24 base cells per axis, frequency errors are approximately
1.162%/0.2955%/0.07451%; doubling resolution reduces that error about fourfold.
CPU/Metal six-cycle traces differ by at most `1.416e-6`. The finest six-cycle
waveform still differs from the continuum trace by 2.70%, explicitly retained
as a spatial-error observation. The maximum normalized balance discrepancy is
`1.702e-7` on Metal. This is a homogeneous closed-mode pilot, not a lossy
transmission line, heterogeneous scattering or antenna qualification. Its
standalone source/runtime/trace evidence is `damped-cavity-pilot-r03.json` in
the contributor's native-admission archive; it does not change the running
combined-source qualification.

```sh
ctest --test-dir build --output-on-failure
build/tests/test_cartesian_subgrid_conductivity --metal
build/tests/test_cartesian_subgrid_native_conductivity --metal
# On a separately built and qualified HIP/NVIDIA runtime:
build/tests/test_cartesian_subgrid_conductivity --cuda
build/tests/test_cartesian_subgrid_native_conductivity --cuda
OPENEMS_SUBGRID_DEVICE=metal python3 python/Tests/Cartesian_Subgrid_Conductivity.py
```

Preserve failed attempts: the early XYZ upload incorrectly required sorted
master-edge traces; validation now permits their actual order without sorting
arithmetic. An earlier CUDA CTest attempt exhausted its wall-time budget on the
unchanged large reference test; a separate retry passed. The later isolated XYZ
attempt completed all 27 tests in one attempt. These are different source/build
identities, not interchangeable green checks.

## Combined source/build qualification

The source was compiled from parent b48e819 with the recorded waveform and
conductivity patches. The later source-tail commit has the same three qualified
numerical waveform blobs; no production/frozen runtime was rebuilt.

- Local: all 25 stages pass, CTest 29/29 and Python 14/14, native fields/ports/
  PML/dipole/polygons and two frozen CPU dipoles. Twenty-one earlier library
  prefixes are unchanged; 11,255 prior numerical arrays remain bit-identical.
  Two harness-only preparation/import failures are preserved; continuation
  checked unchanged source/runtime identities before retaining prior green gates.
- Physical CUDA: all 13 gates pass, CTest 29/29, twelve independent conductivity
  cases per card, five native analytic and four priority/unit controls per card,
  actual native material/Step scripts and retained lossless pair controls in
  both RTX4090/RTX3080 orders. Conductivity runs on one GPU at a time;
  this does not implement distributed material histories.
- Largest sampled conductivity field error 1.733e-6; unchanged 2e-4 field gate.
  Exact numerical tolerances, actual loaded libraries, source inventories,
  terminal states and failed attempts are recorded in the contributor's
  `native-admission/{qualification,cuda-qualification,prior-array-comparison}.json`.
  The ordinary solver implementation is unchanged; its qualified b48e819
  ordinary-engine suite is retained, not a new run.

Frozen combined core SHA-256:

- CPU/Metal: `51b9bd40dc81e3f50c9cdde47dd387a4ec60e02e8a52790ccde0ec1d051ed27e`.
- HIP/NVIDIA: `951ad4afd0d82e6e70e8cc52bdf3282ad98b12e1664115f68b9b9d72e2bf017d`.

## Disjoint conductivity and PML

The joint increment permits loss only outside the actual prepared PML targets
and protected collar/guard supports. Affected electric supports are closed over
each coupled mass block. Magnetic supports use the actual corrected face set.
A bounding-box test or extension insertion order is insufficient. Validation
runs after all operators have prepared their tables and again before engine
creation; direct construction with unprepared PML metadata also rejects.
Mur, overlapping loss/PML and loss inside the protected homogeneous collar
remain rejected. No field-update or PML kernel formula changes in this increment.

The new `test_cartesian_subgrid_conductivity_pml` covers XYZ packets, reverse Y,
both declaration orders and native mapping, exact zero-loss controls, a six-face
corner, real built soft-Dirac sources, coupled 50/75-ohm ports, pure electric or
magnetic loss, and actual-target/collar overlap guards. There are 24 ledger
cases and eight CPU/device field comparisons per physical GPU. Native scripting
also exercises separated loss with PML.

- All 26 local compatibility stages pass, CTest 30/30 and Python 14/14.
  Twenty-two previous runtime prefixes and 11,255 numerical arrays are unchanged.
- All 15 physical CUDA gates pass, including CTest 30/30, the new cases on each
  RTX4090 and RTX3080, and retained lossless two-card controls in both orders.
  Joint conductivity/PML runs on **one device at a time**, not across domains.
- Largest joint sampled CPU/device field error: `7.431e-7`; normalized balance
  error: `1.805e-7`, both below the unchanged `2e-4` gate.
- A single-mesh matched 8 mm lossy slab absorbs about 0.26108 of initial energy;
  analytic attenuation predicts 0.26021, a relative difference below 0.338%.
  This is a packet sanity check, not spatial convergence or broadband reflection
  acceptance. Joint long-time/near-cutoff/antenna coverage remains open.

Frozen joint core SHA-256:

- CPU/Metal: `8d6011bece45f32582964b721d1b6dee8dadd8c81b9a9454b7001245d2fd62bf`.
- HIP/NVIDIA: `da4e5f836e917f8a655bb0f82764c3eaf6948b49eaad95cf5d8fce9a0cb9ec08`.

Contributor evidence: `lossy-materials/disjoint-pml/` qualification, CUDA,
prior-array and source-publication receipts. Earlier failed source-amplitude,
port-setup ordering and driven-energy-envelope fixtures are retained separately;
corrected controls retain the numerical gates and record actual source/builds.

```sh
build/tests/test_cartesian_subgrid_conductivity_pml
build/tests/test_cartesian_subgrid_conductivity_pml --metal
# Separately built/qualified HIP/NVIDIA runtime:
build/tests/test_cartesian_subgrid_conductivity_pml --cuda
```

## Remaining admission gates

Native damped-wave/cavity and short-line propagation need broader analytic and
regular-Yee comparisons, including spatial/time sensitivity and power.
Conducting sheets, physical conductor skin/ADE histories, distributed loss/PML,
DFT/pause and representative lossy PCB/array behavior remain open. These scoped
checks do not establish scalability, end-to-end speedup or antenna accuracy.

## Constant diagonal conductivity — R80 candidate

The midpoint extension also accepts constant Cartesian diagonal electric and
magnetic SI conductivities. This candidate extends R79 diagonal epsilon/mu;
its scoped CPU/reference and single-Metal controls pass independent review.
It is not a rotated/full tensor,
spatial weight, dispersive loss-tangent or roughness implementation.

For C++ samplers, `Coefficients{electric, magnetic}` still assigns scalar values
to every direction. Optional `electricDiagonal` and `magneticDiagonal` arrays
assign three finite nonnegative components. An explicit array's first component
must equal its scalar alias. A zero component is valid, including an explicit
all-zero array with zero alias. Rebuild callers using the changed C++ type.

Native CSXCAD uses its existing material fields, without new XML attributes:

```cpp
material->SetIsotropy(false);
for (unsigned int n = 0; n < 3; ++n)
{
    material->SetEpsilon(epsilonRelative[n], n);
    material->SetMue(muRelative[n], n);
    material->SetKappa(electricConductivitySI[n], n);
    material->SetSigma(magneticConductivitySI[n], n);
}
```

Native admission checks every component, including loss present only in Y/Z.
The existing material-priority/volume quadrature is applied component by
component when needed. Isotropic conductivity keeps its single material lookup.
Each electric-axis trace uses `sqrt(sigma_e[n] * volume / 4)` and the existing
restriction/mass normalization; magnetic faces use `sigma_m[n] / mu[n]`.
The existing midpoint factors, loss ledgers and normalized device recurrence
are reused. No backend kernel or field-update branch is added.

Focused tests preserve the original Maxwell, invariant and isotropic-loss
controls and add analytic per-axis decay, zero-loss directions, isotropic
epsilon/mu with directional loss and invalid sampler/native values. Shaped-loss
tests preserve eighteen isotropic cases and add six directional cases with
independently calculated geometry/priority fractions at refinement interfaces.
CPU82, sixteen bit-exact scalar-parent files, eighteen physical Apple M5 Pro
Metal controls and sixteen exact native field/probe arrays pass independent
review. Maximum native decay field error is 5.911e-7 and maximum shaped-loss
device field error is 1.674e-6 under the retained 2e-4 gates. R81 additionally
qualifies actual CSXCAD Python material/XML roundtrip with 40 native CPU/FP32
field/energy controls and three intended admission rejections. It does not
qualify the full openEMS Python FDTD wrapper. The same XML additionally passes
20 physical Metal field/energy controls and three admission rejections with
the exact accepted CPU/reference timestep. Reproduce the public synthetic
controls using [the XML diagnostic](benchmark/subgrid_material_xml.md).

Broader convergence, physical CUDA/domain and antenna qualification must be
recorded separately; these implementations do not qualify anisotropic absorber
collars or remove the existing geometry/contact/stability guards.
