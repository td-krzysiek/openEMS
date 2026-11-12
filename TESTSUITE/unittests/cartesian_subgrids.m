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
end
