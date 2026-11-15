% SPDX-License-Identifier: GPL-3.0-or-later
function RequestStartupDecision(sim_path, gate_id, accept)
% Accept/reject an already paused --startup-gate run after external checks.
% Use a single controller; this is not a physics or convergence certificate.
if isempty(regexp(gate_id, '^[A-Za-z0-9_-]{1,64}$', 'once'))
    error('Invalid startup gate ID');
end
if ~islogical(accept) || ~isscalar(accept)
    error('accept must be a scalar logical');
end
witness = fileread(fullfile(sim_path, 'startup-gate.xml'));
if isempty(regexp(witness, ['<StartupGate[^>]*id="' gate_id '"'], 'once'))
    error('Startup witness ID mismatch');
end
if exist(fullfile(sim_path, 'STARTUP_CONTINUE'), 'file') || ...
        exist(fullfile(sim_path, 'STARTUP_REJECT'), 'file')
    error('A startup decision already exists');
end
name = 'STARTUP_REJECT';
if accept
    name = 'STARTUP_CONTINUE';
end
temporary = tempname(sim_path);
fid = fopen(temporary, 'w');
if fid < 0
    error('Cannot create startup decision');
end
fprintf(fid, '%s\n', gate_id);
fclose(fid);
[ok, message] = movefile(temporary, fullfile(sim_path, name));
if ~ok
    delete(temporary);
    error('Cannot publish startup decision: %s', message);
end
end
