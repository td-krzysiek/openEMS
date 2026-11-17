function pass = cartesian_subgrids(varargin)
% Check multi-region validation and its native XML structure, without a run.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helperscripts'));
ts_options(varargin{:});
pass = 1;
boxes = [1 2 2 4 6 6; 7 2 2 10 6 6];
ratios = [2 2 4; 4 1 2];
FDTD = InitFDTD('CartesianSubgridBox', '1,2,2,4,6,6', 'CartesianSubgridRatio', '2,2,4');
FDTD = SetCartesianSubgrids(FDTD, boxes, ratios);
regions = FDTD.CartesianSubgrids.Region;
pass = ts_check('region serialization', numel(regions) == 2 && ...
    strcmp(regions{1}.ATTRIBUTE.Box, '1,2,2,4,6,6') && ...
    strcmp(regions{2}.ATTRIBUTE.Ratio, '4,1,2'), 'two independent regions') && pass;
pass = ts_check('single-box attributes removed', ...
    ~isfield(FDTD.ATTRIBUTE, 'CartesianSubgridBox') && ...
    ~isfield(FDTD.ATTRIBUTE, 'CartesianSubgridRatio'), 'no mixed XML definitions') && pass;
FDTD = SetCartesianSubgrids(FDTD, boxes);
pass = ts_check('default ratios', strcmp(FDTD.CartesianSubgrids.Region{2}.ATTRIBUTE.Ratio, '2,2,2'), ...
    'each region defaults independently') && pass;
bad = {boxes(1, 1:5), [1 2 2 4 6 6; 4 2 2 7 6 6], ...
       [1 2 2 4 6 6; 2 3 3 3 4 4], [1 2 2 NaN 6 6]};
for i = 1:numel(bad)
    rejected = false;
    try
        SetCartesianSubgrids(FDTD, bad{i});
    catch err
        rejected = strcmp(err.identifier, 'openEMS:CartesianSubgrid');
    end
    pass = ts_check('invalid region rejected', rejected, 'malformed, touching or nested region') && pass;
end
rejected = false;
try
    SetCartesianSubgrids(FDTD, boxes, [2 3 4; 4 1 2]);
catch err
    rejected = strcmp(err.identifier, 'openEMS:CartesianSubgrid');
end
pass = ts_check('invalid ratio rejected', rejected, 'non-power-of-two ratio') && pass;
% The compact solver uses the existing boundary XML/API. Native PML physics
% is exercised by Cartesian_Subgrid_PML.py; no separate attribute is needed.
FDTD = SetBoundaryCond(FDTD, {'PML_8', 'PML_8', 'PEC', 'PEC', 'PMC', 'PMC'});
pass = ts_check('PML boundary serialization', ...
    strcmp(FDTD.BoundaryCond.ATTRIBUTE.xmin, 'PML_8') && ...
    strcmp(FDTD.BoundaryCond.ATTRIBUTE.xmax, 'PML_8') && ...
    strcmp(FDTD.BoundaryCond.ATTRIBUTE.zmin, 'PMC'), 'existing boundary API') && pass;
FDTD = SetCartesianSubgridChildLines(FDTD, {[0 .13 1], [], [-1 .27 1]});
pass = ts_check('child-coordinate serialization', ...
    strcmp(FDTD.CartesianSubgridChildLines.ATTRIBUTE.X, '0,0.13,1') && ...
    ~isfield(FDTD.CartesianSubgridChildLines.ATTRIBUTE, 'Y'), 'exact explicit and empty axes') && pass;
invalid_axes = {{[0], [], []}, {[0 1 1], [], []}, {[0 NaN], [], []}, {[0 1], []}};
for i = 1:numel(invalid_axes)
    rejected = false;
    try
        SetCartesianSubgridChildLines(FDTD, invalid_axes{i});
    catch err
        rejected = strcmp(err.identifier, 'openEMS:CartesianSubgrid');
    end
    pass = ts_check('invalid child coordinates rejected', rejected, 'malformed or nonincreasing axes') && pass;
end
FDTD = SetCartesianSubgridChildLines(FDTD, {[], [], []});
pass = ts_check('child-coordinate clearing', ~isfield(FDTD, 'CartesianSubgridChildLines'), 'no stale XML element') && pass;
end
