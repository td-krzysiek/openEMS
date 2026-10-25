function ControlRun(sim_path, action, request_id)
% ControlRun: request a pause, continue or graceful stop of an opt-in live run.
% Enable --runtime-control and --decay-gate-id=RUN_TOKEN at launch. Manual
% request IDs start with RUN_TOKEN-. Read runtime-gate.xml before continuing.
% State remains in memory; this is not a process-restart checkpoint.
%
% SPDX-License-Identifier: GPL-3.0-or-later
if ~any(strcmp(action, {'pause', 'continue', 'stop'}))
    error('ControlRun: action must be pause, continue or stop');
end
if isempty(regexp(request_id, '^[A-Za-z0-9_-]{1,64}$', 'once'))
    error('ControlRun: invalid request ID');
end
target = fullfile(sim_path, ['RUNTIME_' upper(action)]);
if exist(target, 'file')
    error('ControlRun: another request is pending');
end
temporary = tempname(sim_path);
fid = fopen(temporary, 'w');
if fid < 0
    error('ControlRun: cannot write request');
end
fprintf(fid, '%s\n', request_id);
fclose(fid);
% As with RequestDFTSnapshot, use a single controller for each simulation.
[ok, message] = movefile(temporary, target);
if ~ok
    delete(temporary);
    error('ControlRun: cannot publish request: %s', message);
end
end
