function FDTD = SetCartesianSubgrids(FDTD, boxes, ratios)
% FDTD = SetCartesianSubgrids(FDTD, boxes, ratios)
% Replace all experimental Cartesian refinement regions.
% boxes: N-by-6 [x0 y0 z0 x1 y1 z1] in mesh drawing units.
% ratios: N-by-3 powers of two from 1 through 16 (default: all 2).
% Each box must align with base-grid lines, remain inside the domain and
% have at least one base-cell gap from every other box. The solver checks
% alignment/margins against the mesh. Nesting, overlap and touching fail.
% Serial engine=basic, homogeneous lossless medium and outer PEC only.
% See CARTESIAN_SUBGRID.md and InitFDTD.

if ~isnumeric(boxes) || ~isreal(boxes) || ndims(boxes) ~= 2 || ...
        size(boxes, 2) ~= 6 || size(boxes, 1) < 1 || size(boxes, 1) > 64 || ...
        any(~isfinite(boxes(:))) || any(any(boxes(:, 1:3) >= boxes(:, 4:6)))
    error('openEMS:CartesianSubgrid', 'Expected 1 through 64 finite ordered boxes (N-by-6)');
end
if nargin < 3
    ratios = 2 * ones(size(boxes, 1), 3);
end
if ~isnumeric(ratios) || ~isreal(ratios) || ...
        ~isequal(size(ratios), [size(boxes, 1), 3]) || ...
        any(~ismember(ratios(:), [1 2 4 8 16]))
    error('openEMS:CartesianSubgrid', 'Expected N-by-3 axis ratios from 1,2,4,8,16');
end
for i = 1:size(boxes, 1)
    for j = 1:i-1
        if ~any(boxes(i, 4:6) < boxes(j, 1:3) | boxes(j, 4:6) < boxes(i, 1:3))
            error('openEMS:CartesianSubgrid', 'Regions must be separated; nesting, overlap and touching are unsupported');
        end
    end
end

regions = cell(1, size(boxes, 1));
for i = 1:size(boxes, 1)
    box = sprintf('%.17g,', boxes(i, :));
    ratio = sprintf('%u,', ratios(i, :));
    regions{i}.ATTRIBUTE.Box = box(1:end-1);
    regions{i}.ATTRIBUTE.Ratio = ratio(1:end-1);
end
if isfield(FDTD, 'ATTRIBUTE')
    for key = {'CartesianSubgridBox', 'CartesianSubgridRatio'}
        if isfield(FDTD.ATTRIBUTE, key{1})
            FDTD.ATTRIBUTE = rmfield(FDTD.ATTRIBUTE, key{1});
        end
    end
end
FDTD.CartesianSubgrids.Region = regions;
end
