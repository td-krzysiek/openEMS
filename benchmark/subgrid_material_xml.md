# Constant-diagonal material XML diagnostic

`subgrid_material_xml.py` uses the actual CSXCAD Python API to construct five
constant-diagonal material controls, write and reload XML, check every material
getter and write a second XML. Cases include directional epsilon/mu,
electric/magnetic loss, zero-loss axes and a completely lossless control.
Three additional files contain deliberately invalid directional coefficients.
Only public, synthetic geometry is used.

`subgrid_material_xml.cpp` reads those files through native CSXCAD and prepares
the existing Cartesian subgrid operator. It checks exact SI constitutive
values and 100-step electric/magnetic midpoint decay against an independent
analytic expression, including the loss ledger. The box has 120 leaves,
nonzero origin and either millimetre or micrometre drawing units. No antenna,
port matching or radiation model is involved.

Use the matching branch's native libraries and CSXCAD Python bindings; the
recorded wrapper source is CSXCAD `cc37cb4ee0ef20c6794acc3361b5f56953eb5189`.
This adds no solver dependency. Build this standalone diagnostic like the C++
tests, with the checkout root and CSXCAD/fparser headers on the include path,
`-std=c++17 -ffp-contract=off -DTIXML_USE_STL`, and the matching openEMS, CSXCAD,
fparser and TinyXML libraries/rpaths. Do not combine different CSXCAD headers
from auxiliary prefixes.

```sh
python benchmark/subgrid_material_xml.py /tmp/subgrid-material-xml-cases
# After compiling subgrid_material_xml.cpp as /tmp/material-xml-control:
/tmp/material-xml-control /tmp/subgrid-material-xml-cases/case-0.xml 0 0 --cpu none
/tmp/material-xml-control /tmp/subgrid-material-xml-cases/case-0.xml 0 1 --reference none
```

Arguments are XML path, case index (0–4), magnetic seed (0/1), execution flag
and expected setup-rejection text (`none` for valid inputs). Run both original
and `-roundtrip` XML for every case, both seeds and `--cpu`/`--reference`:
40 positive controls. The field/energy limits are 2e-10 for CPU and 2e-4 for
FP32 reference. Invalid files use case index 0 and `--cpu`; require exit 1,
`CartesianSubgrid: positive finite diagonal` for invalid epsilon or
`CartesianSubgrid: finite nonnegative diagonal` for invalid kappa/sigma, and
no prepared-grid/device witness. Native geometry admission rejects these
before the helper's later operator-build catch. An exit 1 for a valid case
is always a failure.

Independent R81 review accepts these 40 positive cases and three admission
rejections, every XML component/flag/grid/bound and actual source/runtime/tool
identity. Maximum field and energy-balance errors are 9.013e-7 and 5.118e-7.
CPU/reference and original/reloaded inputs preserve their separate witnesses.
The same 13 XML files also pass 20 actual Apple M5 Pro Metal field/energy controls
and 3 native admission rejections, with the exact CPU/reference timestep and
maximum field/balance errors 9.013e-7/5.118e-7. `--cuda` remains staged, not
physical CUDA/multi-device qualification. This does not qualify the full
openEMS Python FDTD wrapper,
rotated tensors, property weights, broader spatial convergence or an antenna.
