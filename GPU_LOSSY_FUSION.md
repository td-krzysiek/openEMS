# Experimental GPU LossyMetal branch fusion

`OPENEMS_GPU_LOSSY_FUSION=off|apply|all` selects an opt-in HIP/Metal path.
The default is `off`. Set the environment before starting the native process;
Python and Matlab callers can use their normal environment forwarding. No XML
or physics parameter changes. Unknown values fail explicitly.

`apply` keeps the existing eight pre-update launches and combines the eight
field-application launches. `all` also combines the pre-update launches.
The main E update remains between these hooks. All eight branches, coefficients,
histories, FP32 operations and sequential branch-order field subtractions remain.
No reduced-order model, summed/reassociated subtraction, PML fusion, field
migration or compiler/FMA-policy change is included.

Eligibility requires an actual LossyMetal extension, exactly eight voltage-only
non-Lorentz branches and identical, strictly ordered unique position lists.
Unsupported generic Drude/Lorentz configurations retain their existing path.
Histories remain owned by the original domain; a domain with no local surface
positions does no fused work. Memory packing is unchanged.

At verbose level 2 each eligible device reports the mode and owned position
count. Absence of this diagnostic is not proof that fusion was exercised.

Build matching core/backends, run CTest, then compare public fixed-step fields
and probe traces against a pinned parent, for example:

```sh
python3 tests/manual/gpu_lossy_fusion.py --binary /candidate/bin/openEMS \
  --baseline /parent/bin/openEMS --output /new/qualification
```

This needs NumPy/h5py, but no CSXCAD Python bindings. Multiple Metal logical
domains on one GPU test ownership, not physical multi-device performance.
Run actual HIP devices before adoption. Use interleaved uninstrumented A/B
benchmarks and separately retain numerical output comparisons; public controls
do not qualify a production antenna's RF/mesh/material accuracy.
