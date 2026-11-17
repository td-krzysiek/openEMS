# Independent reference for the closed-box XY packet

`subgrid_packet_cavity_reference.py` supplies a continuum modal reference for
the specific vacuum packet in `subgrid_oblique_packet.cpp`. It imports NumPy,
not openEMS. The optional CLI consumes only immutable complete-epoch E/H
timestamps and refuses to overwrite an output or reinterpret collocated H data.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 benchmark/subgrid_packet_cavity_reference.py --controls
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 benchmark/subgrid_packet_cavity_reference.py /tmp/te30.csv \
  --angle 30 --polarization TE --width .25 --output /tmp/te30-analytic.csv
```

Use the matching angle/polarization/envelope from the numerical packet. The
scope is a 30 mm vacuum square, observation times up to81ps and the original
three receivers/common1.6mm cube means. The Z-invariant physical mode makes
the cube's Z integral constant. This is not a generic material/scattering or
antenna reference. Admit RAM/time externally; no device or solver is launched.

## Physical solution and observations

The benchmark labels TE/TM relative to the XY incidence plane: TE has E_z and
Dirichlet XY boundaries at PEC walls; TM has H_z and Neumann XY boundaries.
These labels are not TE_z/TM_z waveguide terminology. PEC/PMC Z caps admit the
corresponding Z-invariant mode. The scalar quantity is E_z or eta*H_z.

Expand the initial rotated Gaussian/carrier in sine/sine modes for Dirichlet,
or cosine/cosine modes for Neumann, using coordinates shifted by L/2. Their
frequencies are `c*sqrt((m*pi/L)^2+(n*pi/L)^2)`. Each initial coefficient
evolves with `cos(omega*t)` for zero initial conjugate field. Neumann includes
the static zero mode: its normalization is1/L per zero index, versus2/L for
positive indices. Dirichlet uses2/L per index. Tensor Gauss–Legendre projection
uses the unchanged continuous seed, with no fitted phase, time shift or gain.

Volume means integrate each sine/cosine analytically, introducing one
`sin(k*a)/(k*a)` factor per transverse axis. They do not reuse the solver's
reconstruction weights or substitute a point sample. Actual numerical E/H
timestamps remain separate. The numerical half-step initialization residuals,
state projection, spatial/time error and observation approximation all remain
in the comparison; especially, the TM code assigns its seed at H(+dt/2).
The reference is the continuous initial-value problem, not a fitted discrete
initial state. Closed-box tails include the physical boundary response and
are not a pure interface-reflection measurement or a completed outgoing pulse.

## Qualification and retained failures

Basis/normalization, a known standing mode, boundary values, exact volume
integrals and the constant Neumann mode pass. Independent coefficient checks
use640-point quadrature; direct scalar-time modal sums agree within6.67e-15.
Resolved cases compare256/384/512-point quadrature and64/96/128/192 maximum
mode indices. Final trace changes are below the original1e-5 quadrature and
1e-4 modal gates. The CLI checks these gates again at the supplied timestamps.

Eight actual CLI runs reproduce16,488 independently reviewed analytical field
values bit-for-bit, with unchanged metadata. Existing-output and incorrect
H/E-stagger controls reject before export. Qualification used NumPy2.5.3 on
Apple M5 Pro; this is a mathematical utility test, not physical GPU execution.

Against this reference, the unchanged48-base-cell numerical results have the
following maximum physical-volume peak discrepancies over three receivers:

| Packet |Local grid |Globally fine grid |
|---|---:|---:|
| TE30°, ratio2/2/2 |3.757%|2.433%|
| TM30°, ratio2/2/2 |**5.855%**|2.499%|
| TE45°, ratio2/2/2 |3.355%|2.002%|
| TM45°, ratio2/4/1 |4.752%|1.369%|
| TE30°, half timestep |3.925%|2.803%|
| TM30°, half timestep |**6.106%**|2.782%|
| TE30°, broader spectrum |**5.555%**|4.465%|
| TM45°, broader spectrum |**6.230%**|2.531%|

Four local volume peak/L2 gates exceed5%; point-H failures also remain. The
earlier local-versus-fine5% agreement is still a valid relative comparison but
does not establish5% continuum accuracy. Neither backend agreement nor energy
conservation removes these failures. Sampling, initialization and propagation
must be separated before diagnosing a solver defect. Separate three-level local CPU controls now reduce the worst volume peak
error from 6.230% to 2.670% and 1.446%, with decreasing peak/L2 errors in
every case and accepted nominal/half-step controls. Physical source, domain,
refinement corners and observation cubes remain fixed; space/time refine
together. Eight original histories remain exact. Three finest-grid point
gates still fail (maximum peak 13.276%), so the volume result does not qualify
point reconstruction. Broad R11 and representative antenna qualification remain open.

The utility changes no numerical library, default guard, frozen model or
production runtime. It is a standalone diagnostic, not an additional CTest
entry or a claim that the existing81-test suite establishes these accuracy gates.
