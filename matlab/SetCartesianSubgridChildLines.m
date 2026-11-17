function FDTD = SetCartesianSubgridChildLines(FDTD, axes)
% FDTD = SetCartesianSubgridChildLines(FDTD, {x, y, z})
% Preserve supplied nonuniform finest-lattice coordinates in drawing units.
% Empty axes retain equal subdivision; {[],[],[]} clears this option.
% Every populated axis must contain max_ratio*(base_intervals)+1 coordinates
% including all exact base anchors. The solver checks this against the mesh.

if ~iscell(axes) || numel(axes) ~= 3
    error('openEMS:CartesianSubgrid', 'Expected three child-coordinate axes in a cell array');
end
children = struct();
names = {'X', 'Y', 'Z'};
for n = 1:3
    value = axes{n};
    if ~isnumeric(value) || ~isreal(value) || ...
            (~isempty(value) && ~isvector(value)) || numel(value) == 1 || ...
            numel(value) > 16384*16+1 || any(~isfinite(value(:))) || ...
            any(diff(value(:)) <= 0)
        error('openEMS:CartesianSubgrid', 'Child coordinates must be finite, strictly increasing vectors');
    end
    if ~isempty(value)
        text = sprintf('%.17g,', value(:));
        children.ATTRIBUTE.(names{n}) = text(1:end-1);
    end
end
if isfield(FDTD, 'CartesianSubgridChildLines')
    FDTD = rmfield(FDTD, 'CartesianSubgridChildLines');
end
if isfield(children, 'ATTRIBUTE')
    FDTD.CartesianSubgridChildLines = children;
end
end
