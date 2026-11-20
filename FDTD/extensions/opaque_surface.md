# Experimental opaque RF surface: connected boxes

`ExperimentalOpaqueSurfaceV1` on an existing CSXCAD Metal property opts into
the passive opaque RF boundary. The attribute format is unchanged:
`g;w1,tau1;w2,tau2;...`, with nonnegative conductance and positive R-L weights
and time constants (up to128 branches). This is not a supplied copper model.

Positive-height, untransformed Cartesian boxes must conform to mesh lines.
Touching or overlapping boxes with identical fits form a single conductor:
internal faces are excluded and exterior half-faces are counted once.
Identical fits may belong to different Metal properties. Different-fit
contacts are rejected, including face/edge/corner contacts. Cylindrical grids,
smooth polygon boundaries and DC transmission through metal are unsupported.
Staircase approximations of diagonal outlines need their own refinement check.

The geometry/EC construction is shared by CPU, tiling and GPU. Device update
equations and FP32 histories are unchanged; CPU boundary histories use FP64.
The physical exterior capacitance participates in timestep selection.
Conforming volume occupancy requires one integer per primary grid cell during
operator construction; it is not retained during time stepping.

Python uses existing CSXCAD calls:

```python
metal = CSX.AddMetal('conductor')
metal.AddAttribute('ExperimentalOpaqueSurfaceV1', passive_fit_string)
metal.AddBox(start1, stop1)
metal.AddBox(start2, stop2)
```

Octave/Matlab can use the same existing property attributes:

```matlab
CSX = AddMetal(CSX, 'conductor');
CSX = AddPropAttribute(CSX, 'conductor', ...
    'ExperimentalOpaqueSurfaceV1', passive_fit_string);
CSX = AddBox(CSX, 'conductor', 10, start1, stop1);
CSX = AddBox(CSX, 'conductor', 10, start2, stop2);
```

The full-simulation regression `python/Tests/Opaque_Surface_Union.py` compares
whole, segmented, overlapping/duplicate and separate-property representations.
Its explicitly fixed4096-step window checks numerical equivalence, not energy
convergence or material calibration. Use `--gpu` to additionally check a device.
The geometry unit test separately covers a concave contact and a cavity.
