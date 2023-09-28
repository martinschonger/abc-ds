% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function res = process_recorded_trajectory(input_fullpath, mode, yz)

arguments
    input_fullpath;
    mode = "record"; % "eval"
    yz = false;
end

opts = detectImportOptions(input_fullpath);
opts.DataLines = [3 Inf];
opts.VariableNamesLine = 2;
full_recording = readtable(input_fullpath, opts);


if ~yz
    rec = full_recording(:, ["time", "O_T_EE_12_", "O_T_EE_13_", "O_T_EE_14_", "x", "y", "dx", "dy"]);
else
    rec = full_recording(:, ["time", "O_T_EE_13_", "O_T_EE_14_", "O_T_EE_12_", "x", "y", "dx", "dy"]);
end
rec_arr = rec{:,:};


timediffs = rec_arr(2:end, 1) - rec_arr(1:end-1, 1);
vels = [];
velsx = rec_arr(2:end, 2) - rec_arr(1:end-1, 2);
velsx = velsx ./ timediffs;
velsy = rec_arr(2:end, 3) - rec_arr(1:end-1, 3);
velsy = velsy ./ timediffs;
vels = [velsx, velsy];
vels = lowpass(vels, 100, 1e3);

% cut off beginning and end


switch mode
    case "record"
        vel_norms = sqrt(sum(vels(:, 1).^2 + vels(:, 2).^2, 2));
        vel_norm_th = 3e-6;
        
        safety_margin_beginning = 1000; % 1
        start_idx = safety_margin_beginning + find(vel_norms(safety_margin_beginning:end) > vel_norm_th, 1);
        safety_margin = 650; % 0
        start_idx = find(vel_norms(1:start_idx+safety_margin) <= vel_norm_th, 1, 'last');
        tmp_end_idx = find(vel_norms(start_idx+1:end) < vel_norm_th, 1);
        if isempty(tmp_end_idx)
            tmp_end_idx = length(vel_norms(start_idx+1:end)) - 1;
        end
        end_idx = start_idx + 1 + tmp_end_idx;
        rec_arr_proc = rec_arr(start_idx:end_idx, :);
        vels_proc = vels(start_idx:end_idx, :);
    case "eval"
        rec_arr_proc = rec_arr(1:end-1, :);
        vels_proc = vels;
end

seconds_factor = 1e3;
rec_arr_proc(:, 1) = (rec_arr_proc(:, 1) - rec_arr_proc(1, 1));
rec_arr_proc(:, 1) = rec_arr_proc(:, 1) ./ seconds_factor;
vels_proc = vels_proc .* seconds_factor;

res = [rec_arr_proc, vels_proc];

end