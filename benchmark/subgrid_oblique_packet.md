# Experimental oblique packet observations

`subgrid_oblique_packet.cpp` is the exact R72-R02 diagnostic, linked against
an already built matching openEMS library. It changes neither the solver nor
the interpolation API. No GPU is used. A normal exit checks finite evolution,
nonzero transmission and the existing modified-energy invariant; it is **not**
an accuracy certificate. Keep and compare the CSV histories separately.

Compile against a matching source/install pair, with normal build dependencies:

```sh
c++ -std=c++17 -O3 -ffp-contract=off -I"$SUBGRID_SOLVER_SRC" \
  benchmark/subgrid_oblique_packet.cpp -L"$SUBGRID_SOLVER_PREFIX/lib" \
  -Wl,-rpath,"$SUBGRID_SOLVER_PREFIX/lib" -lopenEMS \
  -o /tmp/subgrid-oblique-packet
/tmp/subgrid-oblique-packet --sampling-controls
PACKET_OUTPUT=$(mktemp -d)
/tmp/subgrid-oblique-packet 48 30 TE 2 2 2 .25 1.5 .7 "$PACKET_OUTPUT/te.csv"
/tmp/subgrid-oblique-packet 48 30 TM 2 2 2 .25 1.5 .7 "$PACKET_OUTPUT/tm.csv"
```

Arguments are base cells/axis, XY angle in degrees, TE/TM, three ratios,
longitudinal envelope width, wavelength, timestep fraction and CSV destination.
Envelope width and wavelength are normalized by 10 mm. Admitted
diagnostic resolutions are 12/24/36/48, ratios 1/1/1, 2/2/2 and 2/4/1, widths
0.18–0.35, wavelength 1.5 and step fractions 0.7/0.35. The physical domain is a
30 mm cube, with a centred 10 mm refinement box and a shared timestep no larger
than either grid's bound. The largest globally fine reference has 884,736 leaves.
The driver does not enforce an external RAM/time budget; admit it separately.
The fixed observation window is about 60 ps and still contains a finite tail.
It is not a completed-pulse or asymptotic reflection measurement.

TE uses a divergence-free E_z initial packet and a half-step Faraday H seed.
TM uses an H_z packet and PMC Z caps; its initial time-stagger residual is
O(dt²). XY boundaries are PEC. These are closed-box wave packets, not TFSF
plane-wave scattering or a 24 GHz antenna calibration. Initial envelopes,
propagation, interpolation and interface effects all contribute to differences.

## Distinguish the observables

The CSV retains point values and adds the exact mean of the existing local
reconstruction over the same physical 1.6 mm cube at each receiver. Its midpoint
integral is exact for that multilinear local basis; it is not an exact integral
of the unknown physical field. Affine/constant controls check the integration.

H_z is constant in the two transverse directions inside each leaf. A common
point therefore samples different XY cell centres on coarse and fine grids.
R72 point-TM30 errors are 38.47/9.69/17.70% at 24/36/48. Those failed 5% gates
remain failed; a volume mean is a different observable, not a repaired point H.
An unfitted geometrical time-offset diagnostic explains much of the difference.

The added observations leave all 16 original point histories exactly unchanged.
Three-resolution mean-field results, peak-normalized at the transmitted probe:

| XY packet |24 base cells |36 |48 |
|---|---:|---:|---:|
| TE30°, ratio2/2/2 |9.324%|3.802%|2.182%|
| TM30°, ratio2/2/2 |11.080%|6.660%|3.964%|

The 48-cell TE 45° and TM 45° ratio 2/4/1 errors are 2.293% and 3.885%.
Broader-spectrum controls reach 2.580% and 4.548%. Predeclared 5% peak/L2
mean-field gates pass at 48; interpolation of the half-timestep controls changes
TE/TM mean histories by 0.3684%/0.2654% of the fine-time peak, below 1%.
Ratio 1 controls are exact. These are scoped CPU comparisons against the same
core on a globally fine grid, not an independent ordinary-engine replay.

Qualification retains source/library/input and raw-history hashes, external
1.5 GiB child/12.5 GiB reserve guards and separate resource telemetry. Maximum
observed child RSS is 413,220,864 bytes on Apple M5 Pro. No performance claim
or antenna/device admission follows. Other planes/polarizations, broader ratio
families, separated reflection, independent-reference and physical GPU replay,
boundary/window sensitivity and representative antenna checks remain open.
