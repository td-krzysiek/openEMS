function FDTD = SetCartesianSubgridMaxLeaves(FDTD, count)
% Explicit experimental leaf-count budget; default remains 1,000,000.
% This is not a RAM estimate or guarantee. Bound the process externally and
% qualify larger meshes before production. Interface-factor limits still apply.
if ~isnumeric(count) || ~isreal(count) || ~isscalar(count) || ...
        ~isfinite(count) || count ~= floor(count) || count < 1 || count > 357913941
    error('openEMS:CartesianSubgrid', 'Leaf budget must be an integer from 1 through 357913941');
end
if count == 1000000
    if isfield(FDTD.ATTRIBUTE, 'CartesianSubgridMaxLeaves')
        FDTD.ATTRIBUTE = rmfield(FDTD.ATTRIBUTE, 'CartesianSubgridMaxLeaves');
    end
else
    FDTD.ATTRIBUTE.CartesianSubgridMaxLeaves = count;
end
end
