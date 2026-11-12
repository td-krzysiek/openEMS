# SPDX-License-Identifier: GPL-3.0-or-later
"""Native PEC objects, field exclusion, closed-wall isolation and scope guards.

The C++ companion checks finite Cu18 contacts/hollow barrels, interface patterns,
lossless energy, analytic cavity convergence and native Yee equivalence.
These are numerical controls, not a physical-port or antenna qualification.
"""
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_PEC')
os.makedirs(Sim_Path, exist_ok=True)


def fixture(shape='box', offset=(0, 0, 0), scale=1, reverse=False):
    offset = np.asarray(offset)

    def point(values):
        return ((np.asarray(values)+offset)*scale).tolist()

    fdtd = openEMS(NrTS=2048, EndCriteria=0, TimeStepFactor=.8,
                  CartesianSubgridBox=point([-4]*3)+point([4]*3))
    fdtd.SetGaussExcite(13e9, 9e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    csx.GetGrid().SetDeltaUnit(1e-3/scale)
    for n, axis in enumerate('xyz'):
        csx.GetGrid().SetLines(axis, ((np.arange(-8, 9, 2)+offset[n])*scale).tolist())
    csx.AddExcitation('source', exc_type=0, exc_val=[0, 0, 1]).AddBox(point([-6, 0, -2]), point([-6, 0, 2]))
    # Keep base-grid E/H observation stencils clear of the curved PEC surface.
    csx.AddProbe('received', p_type=2).AddPoint(point([6, 4, 2]))
    csx.AddProbe('magnetic', p_type=3).AddPoint(point([6, 4, 2]))
    csx.AddProbe('inside_e', p_type=2).AddPoint(point([0, 0, -2 if shape == 'box' else 0]))
    csx.AddProbe('inside_h', p_type=3).AddPoint(point([0, 0, -2 if shape == 'box' else 0]))

    def dielectric():
        csx.AddMaterial('dielectric', epsilon=2.25).AddBox(point([-8]*3), point([8, 8, 0]), priority=1)

    def conductor():
        metal = csx.AddMetal('conductor')
        if shape == 'box':
            metal.AddBox(point([-2, -2, -4]), point([2, 2, 0]), priority=10)
        elif shape in ['sheet', 'aperture']:
            metal.AddBox(point([0, -8, -8]), point([0, 8, 8]), priority=10)
        elif shape == 'cylinder':
            metal.AddCylinder(point([0, 0, -3]), point([0, 0, 3]), 2.5*scale, priority=10)
        elif shape == 'shell':
            metal.AddCylindricalShell(point([0, 0, -3]), point([0, 0, 3]), 2.5*scale, 1*scale, priority=10)
        elif shape != 'none':
            raise ValueError(shape)

    for build in ([conductor, dielectric] if reverse else [dielectric, conductor]):
        build()
    if shape == 'aperture':
        csx.AddMaterial('air-aperture', epsilon=1).AddBox(point([-2]*3), point([2]*3), priority=20)
    return fdtd, csx


def run(fdtd, label, setup_only=False):
    previous = os.getcwd()
    try:
        fdtd.Run(os.path.join(Sim_Path, label), engine='basic', numThreads=1,
                 exact_endcriteria=True, setup_only=setup_only)
    finally:
        os.chdir(previous)


def data(label, name):
    values = np.loadtxt(os.path.join(Sim_Path, label, name), comments='%')
    assert values.shape[0] > 10 and np.isfinite(values).all()
    return values


for label, options in [('reference', {}), ('reverse', {'reverse': True}),
                       ('translated', {'offset': [31, -23, 17]}), ('micrometre', {'scale': 1000}),
                       ('sheet', {'shape': 'sheet'}), ('cylinder', {'shape': 'cylinder'}),
                       ('shell', {'shape': 'shell'}), ('aperture', {'shape': 'aperture'})]:
    print('PEC case:', label, flush=True)
    fdtd, csx = fixture(**options)
    run(fdtd, label)
    for name in ['received', 'magnetic']:
        values = data(label, name)
        if label in ['reverse', 'translated', 'micrometre']:
            reference = data('reference', name)
            np.testing.assert_allclose(values[:, 0], reference[:, 0], rtol=1e-11, atol=0)
            peak = np.max(np.abs(reference[:, 1:]))
            error = np.max(np.abs(values[:, 1:]-reference[:, 1:]))/peak
            assert peak > 1e-7 and error < 1e-8, (label, name, error)
            print(label, name, 'peak-relative error', error)
        if label == 'sheet':
            assert np.max(np.abs(values[:, 1:])) < 1e-12, 'closed PEC plane leaked to the other cavity'
        else:
            assert np.max(np.abs(values[:, 1:])) > 1e-7, (label, name, 'no exterior response', np.max(np.abs(values[:, 1:])))
    if label in ['reference', 'reverse', 'translated', 'micrometre', 'cylinder']:
        for name in ['inside_e', 'inside_h']:
            assert np.max(np.abs(data(label, name)[:, 1:])) < 1e-12, (label, name, 'field inside PEC')


def reject(label, change, reason=None, base_shape='box'):
    fdtd, csx = fixture(shape=base_shape)
    change(csx)
    try:
        run(fdtd, label, setup_only=True)
    except (ValueError, RuntimeError) as error:
        assert 'CartesianSubgrid' in str(error), str(error)
        if reason is not None:
            assert reason in str(error), (label, str(error), reason)
        print(label, 'correctly rejected:', error)
        return
    assert False, label+': unsupported PEC setup accepted'


reject('subcell-box', lambda c: c.AddMetal('bad').AddBox([-1.1]*3, [1]*3, priority=20))
reject('subcell-sheet', lambda c: c.AddMetal('bad').AddBox([.1, -2, -2], [.1, 2, 2], priority=20))
reject('partial-interface', lambda c: c.AddMetal('bad').AddBox([-4, -3, -3], [-4, -1, -1], priority=20), 'boundary cuts a leaf')
reject('curved-interface', lambda c: c.AddMetal('bad').AddCylinder([0, 0, -4], [0, 0, 4], 2.5, priority=20), 'coarse interface trace')
reject('sphere', lambda c: c.AddMetal('bad').AddSphere([0, 0, 0], 2))
reject('wire', lambda c: c.AddMetal('bad').AddBox([-2, 0, 0], [2, 0, 0]))
reject('axial-gap-unresolved', lambda c: c.AddMetal('bad').AddCylinder([0, 0, -3.1], [0, 0, 3], 2.5, priority=20),
       'cylinder end cuts a leaf')
reject('diagonal-cylinder', lambda c: c.AddMetal('bad').AddCylinder([-2, -2, -2], [2, 2, 2], 1))
reject('source-in-metal', lambda c: c.AddMetal('cover-source').AddBox([-6, -8, -8], [-6, 8, 8], priority=20), 'source touches no unconstrained master edge')
reject('transformed-metal', lambda c: c.AddMetal('bad').AddBox([-2]*3, [2]*3).AddTransform('Translate', [1, 0, 0]))
reject('unresolved-cylinder', lambda c: c.AddMetal('tiny').AddCylinder([.25, .25, -2], [.25, .25, 2], .01, priority=20),
       'no PEC edge resolved', base_shape='none')
reject('all-metal', lambda c: c.AddMetal('filled').AddBox([-8]*3, [8]*3, priority=20),
       'no unconstrained electric edges')
print('PASS: native PEC volumes/sheets/cylinders/shells, priorities, units, field exclusion and rejection guards')
