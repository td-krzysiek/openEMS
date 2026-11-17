# Native fixed-domain execution

Experimental local integration of the existing domain coordinator into normal
CLI/Python runs. It does not qualify the project array or production scaling.
Keep physical geometry in XML and execution placement outside it. An external
planner chooses fixed partitions; the solver does not migrate histories or
automatically rebalance devices.

## Configuration

For a model already containing Cartesian refinement regions:

```sh
openEMS model.xml --engine=gpu --subgrid-domain-axis=x \
  --subgrid-domain-cuts=0 --subgrid-domain-workers=auto@0,auto@1
```

```python
fdtd.Run(sim_path, engine='gpu', subgrid_domain_axis='x',
         subgrid_domain_cuts='0', subgrid_domain_workers='auto@0,auto@1')
```

The existing Matlab option-string path carries the same flags:

```matlab
RunOpenEMS(Sim_Path, 'model.xml', ...
  '--engine=gpu --subgrid-domain-axis=x --subgrid-domain-cuts=0 --subgrid-domain-workers=auto@0,auto@1');
```

The actual CLI and Python paths have executed controls. The Matlab interpreter
has not been executed; this example uses its existing CLI pass-through.

Provide all three options together. Cuts are strictly increasing, finite,
**global drawing coordinates** in the mesh's units; the adapter subtracts the
axis origin before converting to SI. They assign whole cells/constraint
groups to owners, rather than adding mesh planes. Refinement-region corners
and XYZ subdivision factors retain their separate geometry API.

The ordered worker list contains one through64 entries and one fewer cuts:

- `auto@0,auto@1`: explicitly selected device-list indices. Inspect the device
  list on the actual host; indices are not persistent physical identities.
- `reference,reference`: two CPU point-operation reference workers, used with
  `engine=gpu-reference`. This is not the optimized CPU engine.
- `auto,reference`: intentional GPU/reference mixing with `engine=gpu`.
  There is no silent device-to-CPU fallback.
- One worker uses `subgrid-domain-cuts=none`. Boost rejects an empty adjacent
  CLI value; the explicit sentinel works through both scripting paths.

Repeated `auto` entries may select the same best device. They do not request
automatic enumeration of distinct GPUs. Preserve the actual ordered device
identities in the external run manifest. Backend creation restores the prior
global device preference. Invalid/incomplete options, wrong engine/model,
unknown workers and external/unordered assignment cuts reject.

## Native execution and observations

The engine constructs the checked cold plan from the grid, conductivity
extension, PML extension or their explicitly disjoint combination. Unknown
companion extensions reject. Existing factories own the PML hooks; the native
engine does not apply them twice. E/H ghosts carry fields, not material/PML
histories. All workers retain a common fixed timestep.

The coordinator advances complete steps. At existing processing boundaries,
the engine reads the global host field mirror and complete-epoch material and
signed PML-work diagnostics. Normal TD/FD probes, physical-path U/I integrals
and homogeneous unrefined NF2FF surfaces then use that mirror. DFT/NF2FF remain
on the host; this is not a device-gather or device-DFT implementation.

An advance/observation failure invalidates further native stepping, energy
and physical-path U/I sampling. The previous native adapter rejected voltage
and current reads but could return stale fast energy; the strengthened
injected-failure regression fails against the actual R17 parent and passes
on the new build. This is an experimental adapter defect, not an established
ordinary master-engine bug. Low-level `Grid()` access remains a numerical/debug
interface; `ObservationGrid()` enforces the native observation health guard.

No external collective pause/resume/snapshot API is enabled by this adapter.
Normal complete-batch processing is not such an API. PML is the checked
classical convolution extension, not UPML/CFS parity. The global host topology,
field mirror and per-state ownership maps still require scaling work.
The separately qualified [global row stream](CARTESIAN_SUBGRID_DOMAINS.md#streamed-global-domain-preparation)
removes the temporary pair of global expanded operators, without changing
local upload arrays or runtime kernels in that R19 increment. The later
[local row package](CARTESIAN_SUBGRID_DOMAINS.md#streamed-local-domain-preparation)
also removes expanded local uploads, with unchanged runtime kernels and
separately recorded qualification.

## Qualification

R18-r01's full local replay passes20 stages, CTest44/44, Python14/14 and11,255
exact retained numerical arrays. Focused C++ reference/Metal native tests cover
XYZ assignment, nonzero coordinate origins, mm/um units, disjoint bulk loss/PML
and repeated32-step processing batches. Six cases per backend retain exact
one-device FP32 fields and integrated energy; signed-work tolerances remain
unchanged.

Actual Python runs compare identical native XML, TD/FD fields and two PEC
two-port fixtures. Reference domains are exact; Metal complex-S11 differences
are below6.58e-8 and input-Z differences below9.44e-6ohm. A public PEC dipole
with two homogeneous unrefined surfaces also passes on two logical Metal
workers against the frozen same-input single-reference run: deltaZ3.75e-6ohm,
deltaGamma2.67e-8 and complex far-field difference2.24e-7 relative to peak.
Port/power and TD/FD timing gates are retained. This is not a PCB-array test.

R18-r02 changes only the option parser and its tests to add `none`. Its own
reference/Metal native fields, actual CLI/Python ports, intentional Metal/CPU
reference mixing, six CLI rejection controls, one-worker CLI control and all44
CTest cases pass. The numerical adapter/physics source is identical to R01;
its complete replay and dipole evidence remain separately identified. All
source/runtime identities are checked, rather than silently mixing prefixes.

Local R18-r01 core SHA256:
`97b15556d75964497a9d9a0c3927ad545fdf80a27d0313b578c65ce989e96778`.
Local R18-r02 core SHA256:
`b0b157b27840bde3debef17ff1f914c8b9a049905dbed9b7716059fd8a404de4`.

Physical R18-r01 passes48 CUDA gates and CTest44/44 on the authorized
RTX4090/RTX3080 pair, including native C++ domains in both device orders and
16 actual Python TD/FD/U/I cases. Its terminal exit is0, without OOM; source
and runtime identities are rechecked. CUDA core SHA256:
`8cc678373160b2433b178f6827771f712ce96ef596612f2920b1e38b08c0f6a3`.
The parser-only R02 follow-up also passes seven physical CUDA gates and
CTest44/44, including both pair orders,16 native Python cases and actual
one-worker CLI sentinel validation. Earlier low-level physical-domain
evidence does not substitute for native processing tests. Two logical Metal
workers are not two physical GPUs.

Evidence is under the separately frozen `native-runtime-r18/` and
`native-runtime-r18-r02/` studies. No frozen antenna, production runtime or
allocation guard is changed. Required PCB/conductor accuracy, interface/host
scaling, actual end-to-end savings and representative array acceptance remain
open, as do broader master features and collective external controls.


## Exact native child coordinates

`SetCartesianSubgridChildLines([x, y, z])` optionally retains unequal physical
child positions instead of equal subdivisions. Coordinates use mesh drawing
units and the same origin/unit conversion as the base axes. A populated axis
contains `max_ratio * base_intervals + 1` finite, strictly increasing entries,
including every base anchor **exactly**. Empty axes retain equal subdivision;
three empty axes clear the setting. The selected region ratios remain powers
of two and non-nested/separated; child coordinates do not change that topology
or admit touching covers. Invalid requests fail rather than snapping anchors.
Configure before operator preparation. Rebuild matching C++/Python bindings.

C++ uses `SetupCartesianSubgridChildLines`; Python also accepts the constructor
keyword `CartesianSubgridChildLines`. Optional XML stores comma-separated X/Y/Z
attributes in one `CartesianSubgridChildLines` child of FDTD. The Matlab helper
`SetCartesianSubgridChildLines(FDTD, {x,y,z})` is supplied; actual Octave/Matlab
execution remains unqualified.

R27-r02's corrected core fixture passes twelve owning-geometry/256-step-energy
cases and twelve single/domain reference cases; six no-port cases also match
an independent FP32 oracle. The first degenerate symmetric-field diagnostic
fails on both old and new meshes and is retained separately, not interpreted
as a solver defect. R28 adds XML/API/units and matching Cython bindings:53/53
Linux CPU tests,11 Python API tests, three native E/H runs, two invalid mesh
setups, sixteen exact frozen-parent operator/field files and eight native
reference cases pass. The numerical grid source remains byte-identical to
R27-r02; the actual R28 Linux `libopenEMS.so` (including the new adapter) has
SHA256 `78a0ee87b0a61bc3e1460af1c2cad6261c29ae431fc0e6773f5cd4e2c6517ea9`.

Eighteen physical CUDA gates pass on the RTX4090/RTX3080 pair, including48 new
child-coordinate cases across each device and both orders,24 exact domain
metadata/field comparisons, retained loss/PML controls and16 native TD/FD/U/I
cases. An independently hashed physical-only driver changes only backend
selection, label and scope wording in the native child E/H script. It runs
three cases without changing geometry, timestep or thresholds. Source,
runtime, test driver and terminal exit0/no-OOM were independently rechecked.
These are scoped controls, not Metal, full-array, RF-convergence or production
qualification. Native unequal-coordinate loss composition is separately
tested by R30; it is not inferred from equal-subdivision loss controls.

`tests/test_cartesian_subgrid_device_child_loss.cpp` supplies that separate
composition fixture. Fourteen CPU/reference and58 physical CUDA cases pass
against the same R28 library, covering each GPU and both pair orders. It
requires actual unequal child lengths, preserves synthetic surface-path
anchors and aligns native18/35um volume-copper XY limits to the supplied
knots. All original field/ledger/passivity thresholds remain. Single/domain
fields are exact; maximum double-field error is2.033e-6 and normalized energy
balance error1.543e-7. Source/runtime/fixture/binary identities and terminal
exit0/no-OOM are independently rechecked. This unchanged standalone fixture
is now registered in CTest; the combined registration passes in R31's
56-test superset. Physical Metal remains pending. It is not a PCB RF-convergence result.
Allocation/interface guards remain unchanged.


## Native rectangular resistive ports

Native lumped resistors may span a resolved transverse sheet or volume rather
than only a one-dimensional line. The experimental adapter retains the native
bounds and positive R, checks complete independent local-edge support, and
normalizes dual-area/edge-length conductances to the requested total resistance.
Hanging traces, unresolved endpoints, missing/opaque support and overlapping
ports reject. `Caps=true` adds tangential PEC electrodes at both longitudinal
ends. The existing line-port path and allocation guards remain unchanged.
This is geometry/preparation support, not an RF-calibrated antenna port.

R34-r03 passes Linux CPU CTest59/59, sixteen exact frozen-R33 operator/field
files and retained native reference controls. Its new fixture covers48 cases:
three directions, uniform/graded axes, sheet/volume footprints, matching-grid
ordinary-Yee fields or actual internal refinement, and CPU/requested backend.
The geometry-normalized resistor coefficients are compared with the actual
ordinary operator's VV/VI coefficients, not a duplicated formula alone.
Shared timesteps remain below both operators' CFL limits.

All26 physical gates pass on RTX4090/RTX3080, including48 new port cases on
each GPU and both pair orders. Maximum normalized coefficient, field and
balance errors are2.027e-7,5.304e-7 and9.723e-8. The retained24 exact domain
field cases, child-coordinate/loss controls and16 native Python TD/FD/U/I
cases pass. Actual source/runtime hashes and terminal exit0/no-OOM were
rechecked; CUDA libopenEMS SHA256 is
`970f273bd64fe96721ea86c6fb9613b3ae015d3f59525f6ceb729b460bebbf59`.
Evidence: the separately frozen `native-box-ports-r34-r03/` study.

Earlier test-only compile, margin and collector-marker failures remain recorded;
the accepted run uses a fresh source-before-build prefix. Metal, driven
microstrip/Cu/PML composition, RF convergence, full-array memory and antenna
qualification remain open. General path quadrature can retain roundoff-sized
extra support at some graded cuts; the separately prepared R35 descriptor
follow-up diagnoses that strict-admission limitation without changing R34's
frozen inputs or relaxing numerical thresholds.


## Retained native line normalization diagnosis

The earlier native one-dimensional mapper's equalR/N assignment is not
ordinary-Yee parity on graded longitudinal intervals. A minimal GPU-free
comparison against actual ordinary `(1-VV)/VI` coefficients reproduces a
maximum relative per-cell conductance difference0.53846154 in six graded
cases (three axes, both orientations), while six uniform controls differ
only1.355e-8 after FP32 metadata conversion. Both operators use identical
knots and a valid shared timestep; no FDTD or RF accuracy follows from this
coefficient-only test. This is an experimental adapter discrepancy, not an
upstream engine defect or a53.85% antenna-impedance prediction.

The native finite-width rectangular path above is separate. Frozen curved
R23 line ports span eight uniformZ intervals inside constantXY refinement;
read-only XML checks therefore do not support assigning their spatial RF
failure to this graded-longitudinal counterexample. Preserve those failures.
A separately frozen R36 candidate distributes lineR by physical segment
length, retaining local voltage functionals, reference-port limits and the
explicit C++ circuit API. It is being built/tested, not promoted. Required
checks include actual-old-library failure, ordinary coefficients, retained
interface/passivity/partial-edge controls and physical backends. Evidence:
`native-line-normalization-r36-diagnosis/`; the correction is qualified below.

## Physical-length native line normalization

R36 changes only the native one-dimensional adapter: segment resistance is
`R * abs(segmentLength) / totalLength`. The ordinary matching tensor grid has
constant transverse dual area, so its normalized series conductance implies
this physical-length weighting. Equal R/N is equivalent only on equal
longitudinal intervals. The local constrained voltage functionals, rectangular
ports and explicit C++ circuit semantics remain unchanged.

The new12-case actual ordinary-VV/VI coefficient comparison covers all axes,
uniform/graded lines and reversal. It fails the first graded case against the
actual frozen R35 library after two uniform controls. The independent manual
path fixture retains fields, losses, passivity and partial-edge guards.
CPU61/61 and16 exact frozen-parent files pass. Physical28-gate qualification
on RTX4090/RTX3080 in both orders retains24 exact domain cases, losses/PML,
native boxes and16 native Python cases. The actual qualified Linux core is
`8a290446557f83301e9268a5325639747a319a38ed737c523e2c377164865a92`.
Terminal exit0/no OOM and380 archived files/758 source hashes were rechecked
before promoting exact numerical/test bytes. Project evidence is in
`experiments/openems-cartesian-local-subgrid-20261112/lossy-materials/native-line-normalization-r36/`.
Metal, driven copper RF, array scaling and production adoption remain open.


## Reusing resolved native port support

R35-r02 retains each already verified native rectangular-port edge ID, one
incident leaf and its local edge coordinates. Overlap checks and admittance
preparation use this immutable support rather than repeating the general
full-leaf voltage-path integration for every edge. Resistance, order, physical
endpoints, source geometry and numerical update kernels remain unchanged;
the older owning `DistributedResistor` API remains available. Rebuild C++
callers against the matching headers. Descriptors require unchanged topology.

The original test incorrectly required the general integrator to have exactly
one support term. An unchanged-library diagnosis finds an extra quadrature
term3.6266e-16 relative to the identical main coefficient at a rounded cut.
Fresh r02 changes only the test oracle: keep the original2e-14 bound on the
complete functional's relativeL1 difference, with wrong-main/wrong-support
negative controls. Exact direct topology and all hanging/missing-edge guards
remain. This avoids a strict-cardinality false rejection, rather than silently
admitting unsupported geometry. Original failure and diagnostic are retained.

CPU60/60,48 new descriptor controls,16 exact actual-R34 frozen files and all26
physical CUDA gates pass on4090/3080 in both orders. These include48 new
native-box cases per hardware lane,24 exact domain metadata/field cases,
retained child loss/PML and16 native Python observations. The actual library
SHA256 is `a94dd254ddffb78c38988e4f290cd679f04539378ed7334786376465871894be`.
Terminal exit0/no-OOM and all377archived files/source identities are rechecked.
Metal, driven RF, full-array memory and antenna acceptance remain open.

A GPU-free native294-trace port benchmark uses identical binary-exact geometry
on the actual frozen R34/R35 libraries. Three interleaved repetitions time
complete CalcECOperator, before untimed prepared-loss fingerprint extraction:

| Leaves | R34 median | R35 median | Preparation ratio |
|---|---:|---:|---:|
|16,192|2.803563s|0.218491s|12.83x|
|118,720|20.918891s|1.580225s|13.24x|

Within each size, prepared row/rate/timestep fingerprints and counts are exact
across libraries/repetitions; source/runtime/binary hashes are bound separately.
PeakRSS is nearly unchanged. This shared-load small native setup comparison
is not a master-versus-branch FDTD benchmark, full-array wall-time prediction
or evidence that all preparation is13x faster. The inherited no-port domain
benchmark is not used to establish this gain. Evidence: the separately frozen
`resolved-port-edges-r35-r02/` study and `port-benchmark-summary-r01.json`.
