function FDTD = SetCartesianSubgridMaterialSamples(FDTD, samples)
% Explicit scalar leaf-volume quadrature for shaped isotropic materials.
% 0 (default): strict aligned boxes; 2/4/8: samples per axis in intersecting leaves.
% This is not conformal integration. Check quadrature and mesh convergence.
if ~isnumeric(samples) || ~isreal(samples) || ~isscalar(samples) || ...
        ~isfinite(samples) || ~ismember(samples, [0, 2, 4, 8])
    error('openEMS:CartesianSubgrid', 'Material samples must be 0, 2, 4 or 8 per axis');
end
if samples == 0
    if isfield(FDTD.ATTRIBUTE, 'CartesianSubgridMaterialSamples')
        FDTD.ATTRIBUTE = rmfield(FDTD.ATTRIBUTE, 'CartesianSubgridMaterialSamples');
    end
else
    FDTD.ATTRIBUTE.CartesianSubgridMaterialSamples = samples;
end
end
