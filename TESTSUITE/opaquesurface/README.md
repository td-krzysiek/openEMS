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
