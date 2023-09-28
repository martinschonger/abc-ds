% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [Data, Target, indivTrajStartIndices, timestamps, stuff] = recorded_trajectories_to_refdata(input_fullpaths, num_subsamples, mode, yz)

arguments
    input_fullpaths;
    num_subsamples = 100;
    mode = "record"; % "eval"
    yz = false;
end

Data = [];
stuff = [];
indivTrajStartIndices = [1];
timestamps = {};
for i = 1:length(input_fullpaths)
    res_tmp = process_recorded_trajectory(input_fullpaths{i}, mode, yz);
    disp(size(res_tmp, 1));
    res_tmp = res_tmp';
    if num_subsamples > 0
        idx = 1:size(res_tmp, 2);
        didx = linspace(min(idx), max(idx), num_subsamples);
        res = zeros(size(res_tmp, 1), num_subsamples);
        for j = 1:size(res_tmp, 1)
            res(j, :) = interp1(idx, res_tmp(j, :), didx);
        end
    else
        res = res_tmp;
    end

    Data = [Data, [res([2, 3, 9, 10], :)]];
    stuff = [stuff, [res([5, 6, 7, 8], :)]];
    indivTrajStartIndices = [indivTrajStartIndices, 1 + size(Data, 2)];
    timestamps{end+1} = res(1, :);
end

M = 2;
Target = mean(Data(1:M, [indivTrajStartIndices(2:end) - 1]), 2);

end