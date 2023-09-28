% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [Data, Target, indivTrajStartIndices, Timestamps] = generateRefData(fun, initial_set_sampled)

timestamps = {};
simulated_trajectories = {};
Opt = odeset('events', @isAtOrigin, 'reltol', 1e-7, 'abstol', 1e-7); % stop integrating at some boundary


[M, N] = size(initial_set_sampled);

for n = 1:N
    [t, y] = ode45(@(t, xi) fun(xi), [0, 60], initial_set_sampled(:, n), Opt);
    timestamps{end+1} = t';
    simulated_trajectories{end+1} = y';
end


Data = [];
indivTrajStartIndices = [0];
for n = 1:N
    xdiff = simulated_trajectories{n}(:, 2:end) - simulated_trajectories{n}(:, 1:(end -1));
    tdiff = timestamps{n}(:, 2:end) - timestamps{n}(:, 1:(end -1));
    tdiff_extended = repmat(tdiff, M, 1);
    vel = xdiff ./ tdiff_extended;
    vel = [vel, zeros(M, 1)];

    Data = [Data, [simulated_trajectories{n}; vel]];

    indivTrajStartIndices(end+1) = size(Data, 2);
end
indivTrajStartIndices = indivTrajStartIndices + 1;

Target = zeros(M, 1);

Timestamps = timestamps;

end