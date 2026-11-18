# Native magnetic point-readout diagnostic

The existing `Magnetic`/`MagneticSample` reconstruction is linear in a
component's normal direction and constant in its two transverse directions
inside a leaf. For the Z-invariant packet H_z this means that requesting an
XY point uses the transverse leaf centre. A common requested point can
therefore describe different effective positions on two meshes. This is the
documented reconstruction, not a newly established propagation defect.

`subgrid_point_readout.cpp` records the actual native requested/centre/leaf
coordinates. It runs no FDTD steps. Eighteen affine controls independently
verify both direct and weighted H readout, including linear normal-Z variation.
Six native geometries cover N48/N72/N96 with ratios 2/2/2 and 2/4/1. The
1,200,000-leaf allowance is explicit in these fixtures only.

For the transmitted 45-degree receiver the longitudinal centre offset is
+129.126 / -18.188 / +129.126 micrometres on N48/N72/N96. The two dyadically
related meshes can have the same projected offset even though their cells
shrink. This explains why point-error improvement need not be monotonic.

The independent continuum cavity modes are evaluated again at these actual
centres, using the original H timestamps. There is no fitted delay, phase or
amplitude, no altered source, no replacement point result and no change to
fields or libraries. With one fixed original requested-point normalization,

```text
observed - exact(requested)
  = [exact(native centre) - exact(requested)]
  + [observed - exact(native centre)].
```

The residual on the last line still combines propagation, discrete initial
state/time staggering and any remaining observation approximation. It is not
a pure interface-reflection measurement. The original point gates remain
failed. The N96 TM30 L2 error is 5.0009%, even though its peak error is
4.9402%; a sub-5% peak alone does not satisfy both gates. At N96 the
largest requested-point peak error is 13.276%, while the
largest native-centre residual is 1.379%. This diagnoses the sampling effect;
it does not qualify physical point reconstruction.

| TM packet | N48 requested / centre residual | N72 | N96 |
|---|---:|---:|---:|
| 30 degrees | 16.344% / 4.610% | 2.940% / 1.993% | 4.940% / 1.124% |
| 45 degrees, anisotropic ratio | 13.068% / 3.757% | 3.631% / 1.559% | 10.419% / 0.909% |
| 30 degrees, half step | 16.454% / 4.914% | 2.949% / 2.129% | 4.860% / 1.171% |
| 45 degrees, broader spectrum | 17.348% / 5.624% | 5.697% / 2.386% | 13.276% / 1.379% |

N48 retains a 5.624% centre residual in the broader-spectrum case; geometric
explanation is not a universal accuracy certificate. The original data and
normalizations remain unchanged, and each row uses the maximum across three
receivers. Both modal/quadrature convergence gates pass at every new centre.

## Reproduction

Compile against a matching source/install pair and its normal dependencies:

```sh
c++ -std=c++17 -O3 -ffp-contract=off -I"$SUBGRID_SOLVER_SRC" \
  benchmark/subgrid_point_readout.cpp -L"$SUBGRID_SOLVER_PREFIX/lib" \
  -Wl,-rpath,"$SUBGRID_SOLVER_PREFIX/lib" -lopenEMS -o /tmp/subgrid-point-readout
READOUT_OUTPUT=$(mktemp -d)
/tmp/subgrid-point-readout 96 241 "$READOUT_OUTPUT/geometry.csv"
```

The first argument is N48/N72/N96; the second is 222 (30 degrees) or 241
(45 degrees). The destination must be fresh. Columns retain all three
requested coordinates, actual centre/start/length, leaf identity/count and
affine direct/weighted checks. Native setup resources must be admitted
externally; a successful exit is not a propagation or antenna qualification.

`subgrid_local_packet.cpp` is the exact separately tested local-only packet
driver used for the three-level convergence study. Build it in the same way.
It evolves only the local field lattice for about 80 ps and retains all
original N48 local columns exactly, before N72/N96 refinement. For example:

```sh
/tmp/subgrid-local-packet 48 45 TM 2 4 1 .25 1.5 .7 \
  "$READOUT_OUTPUT/packet48.csv" 2.0248077747204605e-13
/tmp/subgrid-local-packet 96 45 TM 2 4 1 .25 1.5 .7 \
  "$READOUT_OUTPUT/packet96.csv" 2.0248077747204605e-13
```

The final argument is the frozen N48 timestep; the driver scales it by 48/N
and independently rejects any coupled-CFL violation. Other retained N48
steps are 3.443786467501156e-13 s (TE, ratio222, fraction0.7),
3.0372116620806914e-13 s (TM, ratio222, fraction0.7); half-step controls
use half these values. The local-only driver is limited to the documented
N/angle/ratio/source combinations, with explicit allocation and energy gates.
It does not silently raise the library's global guard or claim accuracy from
energy conservation.

To inspect the fixed-coordinate continuum diagnostic, use `Modes`, `profile`
and `evaluate` from `subgrid_packet_cavity_reference.py`: requested XY and
native-centre XY are passed as distinct positions at the same raw `time_h`
values. Use 192 modes/512 quadrature points, checking against 128/512 and
192/384 as documented there. Preserve the requested-point error as well as
the diagnostic residual. The original standalone CLI intentionally retains
its original fixed receivers.

These sources are byte-exact copies of executed CPU fixtures; they change no
numerical solver, interpolation default, frozen input or production runtime.
The 24 local packet cases and 18 affine geometry checks are standalone
evidence, not additions to the existing 81-entry CTest suite. Physical GPU,
broader scattering and representative antenna qualification remain separate.
