# Experimental main-update coefficient storage policy

> Historical runtime documentation. The current cell-model integration replaces
> these kernel controls; see [migration and qualification](UPSTREAM_GPU_CELL_MODEL.md).

`OPENEMS_GPU_MAIN_COEFF=auto|raw` applies to HIP and Metal. The default `auto`
retains the existing bit-exact coefficient-set compression heuristic and raw
fallback. `raw` bypasses table construction for the main E/H update and uploads
the existing full coefficient arrays. PML/material tables and all equations
are unchanged. Unknown values fail explicitly. Set it before native startup,
including when launching through Python/Matlab; no new XML attribute is used.

At verbose level 2, raw mode reports its main-coefficient byte count for each
device, including ghost storage. Raw mode can consume substantially more VRAM;
preserve actual allocation/headroom gates. This is a measurement option, not a
claim that compression is slow or that raw is universally preferable.

Build matching core/plugins. CTest covers policy rejection and bit-exact 16/32bit
table reconstruction. Fixed-step native controls compare both representations,
also combined with optional LossyMetal fusion:

```sh
python3 tests/manual/gpu_coeff_policy.py --binary /candidate/bin/openEMS \
  --baseline /parent/bin/openEMS --output /new/qualification
```

Benchmark the actual workload and devices interleaved. Keep setup, FDTD and
observation timings distinct; default-off fusion permits independent A/B tests.
No production runtime, active job or mesh is changed by exposing this option.
