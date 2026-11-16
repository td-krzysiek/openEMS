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
field mirror and expanded domain preparation still require scaling work.

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

Physical R18-r01 CUDA48-gate/CTest44 qualification is running on the authorized
RTX4090/RTX3080 pair. The separately built R02 focused follow-up waits for that
terminal result. These live matrices are **not accepted physical native-domain
qualification yet**. Earlier low-level physical-domain evidence does not
substitute for native processing tests. Two logical Metal workers are not two
physical GPUs.

Evidence is under the separately frozen `native-runtime-r18/` and
`native-runtime-r18-r02/` studies. No frozen antenna, production runtime or
allocation guard is changed. Required PCB/conductor accuracy, interface/host
scaling, actual end-to-end savings and representative array acceptance remain
open, as do broader master features and collective external controls.
