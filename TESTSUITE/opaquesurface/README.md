# Opaque-surface timestep regression

This is a numerical regression for the experimental opaque RF boundary, not
material calibration. The standalone XML describes a uniform smooth-copper
microstrip with35µm physical strip/ground height, a graded Cartesian mesh,
passive isotropic dielectric and passive surface R-L coefficients. No antenna
layout, supplier data or external Python RF libraries are needed.

The previous extension replaces boundary capacitance after timestep selection.
On this grid its default1.59336e-14s step produces exponential energy growth
on both standard and tiled engines; PEC on the same grid decays normally.
The correction inserts the same physical exterior capacitance into the EC
arrays before the existing timestep calculation. It adds no stabilizing
capacitance and does not change the R-L response or material inputs.

```sh
OPENEMS_BIN=/absolute/path/to/patched/openEMS python3 test_stability.py
```

Optional `OPENEMS_BASELINE_BIN` runs the unfixed experimental boundary as a
negative control. Every run has a1200-second wall-time limit and an energy
watchdog. The threshold1e-5 is far above the normal driven energy (~1e-12).
Successful cases must reach−70dB, avoid the timestep cap, and retain finite,
nonzero probes. The no-surface control must retain the original timestep.
The standard and tiled engines must select the same timestep to log precision.

The checks are opt-in because they perform full field decay on about1.6million
grid nodes. Boundary geometry remains restricted to disjoint conforming boxes;
these checks do not qualify arbitrary connected copper, roughness, DC transfer
or a full antenna. Frozen simulation results remain unchanged.

## GPU comparisons

`test_gpu.py` uses small coarse versions of the same line, rotated through all
three axes, with unequal fit orders (including 128 branches) and a PEC control.
It compares the complete probe waveforms and complex DFTs at 22, 24.125 and
26 GHz against the basic CPU engine. Peak waveform error must be below 1e-4;
relative complex error must be below 5e-4. The CPU GPU-reference backend also
runs, and is expected to reproduce the CPU boundary without FP32 device casts.
These fixed-length cases test indexing and state updates, not field decay or
spatial mesh convergence. No CSXCAD Python bindings are required.

```sh
OPENEMS_BIN=/absolute/path/to/openEMS OPENEMS_TEST_GPU=1 \
  OPENEMS_GPU_BACKEND=Metal python3 test_gpu.py -v
OPENEMS_BIN=/absolute/path/to/openEMS OPENEMS_TEST_GPU=1 \
  OPENEMS_GPU_BACKEND=Metal python3 test_stability.py -v
```

The optional backend string is checked against the actual solver log; use the
reported HIP backend name when checking HIP/CUDA. With GPU enabled, both tests
require all extensions on the device and reject host fallback. The stability
test additionally compares the entire pulse with a 1e-4 relative peak tolerance
and the 24.125-GHz complex probe voltages within 3e-5. It requires full decay to
-70 dB for each of standard CPU, tiled CPU and GPU.

`OPENEMS_TEST_OUTPUT` retains XML, logs, probes and JSON comparison records.
Choose a fresh directory for each run; existing output is not overwritten.
GPU states/coefficient buffers use FP32, while the CPU boundary state uses
FP64. A successful short-line comparison is not proof of antenna-scale loss
accuracy. Metal and HIP must be tested on their respective hardware; the
GPU-reference backend alone does not exercise their kernels.
