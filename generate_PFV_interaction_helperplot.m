% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [fig] = generate_PFV_interaction_helperplot(output_root_path, exp_id, sim_id, plot_id, yz)

arguments
    output_root_path;
    exp_id;
    sim_id;
    plot_id;
    yz = false;
end


frames = 1000;
linewidth_global = 1;
linewidth_global_major = 2;
linewidth_global_minor = 0.5;

[options, result, f, V, B] = read_exp(exp_id);

exp_list = options.dataset_opts.exp_list;
exp_list_cellarr = fileread(exp_list + "reference_trajectories.txt");
exp_list_cellarr = regexp(exp_list_cellarr, '\r\n|\r|\n', 'split');

[Data, Target, indivTrajStartIndices, timestamps] = recorded_trajectories_to_refdata(exp_list_cellarr, frames, "record", yz);
M = 2;
shift = Target;
Data(1:M, :) = Data(1:M, :) - shift;
Target = [0; 0];
rd = RefData;
rd.directInit(Data, Target, indivTrajStartIndices, timestamps, true);


x_min = min(rd.Data(1, :));
x_max = max(rd.Data(1, :));
x_range = x_max - x_min;
y_min = min(rd.Data(2, :));
y_max = max(rd.Data(2, :));
y_range = y_max - y_min;
xy_range = max(x_range, y_range);
xy_range = xy_range * 1.75;
x_lowerlim = x_min - 0.5 * (xy_range - x_range);
x_upperlim = x_max + 0.5 * (xy_range - x_range);
y_lowerlim = y_min - 0.5 * (xy_range - y_range);
y_upperlim = y_max + 0.5 * (xy_range - y_range);
limits = [x_lowerlim, x_upperlim, y_lowerlim, y_upperlim];


initial_set_center = rd.xi0_mean;
initial_set_radius = 0.05;
unsafe_set_center = [1; 1];
unsafe_set_radius = 0.1;

workspace = [[-1; -0.5], [0.5; 1]];


xi = sdpvar(rd.M, 1);

% initial/safe set
r1 = initial_set_radius;
r21 = r1 * r1;
center1 = initial_set_center;
initial_set = {};
initial_set{1} = r21 - sum((xi-center1).^2, 1);

% unsafe set
unsafe_set = {};
[unsafe_set_coefs, unsafe_set_monomials] = file2poly('datasets/robot_recordings/two_obstacle/unsafe_set_poly.json', xi);
unsafe_set{1} = dot(unsafe_set_coefs, unsafe_set_monomials);


plot_dim1 = 1;
plot_dim2 = 2;

[fig, ax] = setup_figure(2*3*3.41275152778, 3*3.41275152778, false, false, false);


set(ax,'ztick',[],'xcolor','black','ycolor','black');
set(ax,'color',[0 0 0]);
fontsize(fig, 10, "points");
box on;


axis_limits = axis;


resolution = 0.01;
[X, Y] = meshgrid(axis_limits(1)-resolution:resolution:axis_limits(2)+resolution, axis_limits(3)-resolution:resolution:axis_limits(4)+resolution);
XY = [X(:)'; Y(:)'];



% actual trajectories based on executing the DS on the robot
sim_path = sim_id;
[Data3, Target3, indivTrajStartIndices3, timestamps3] = recorded_trajectories_to_refdata({sim_path}, frames, "record", yz);
M3 = 2;
shift3 = shift;
Data3(1:M3, :) = Data3(1:M3, :) - shift3;
Target3 = Target3 - shift3;
rd3 = RefData;
rd3.directInit(Data3, Target3, indivTrajStartIndices3, timestamps3, true, rd.state_maxnorm, rd.vel_maxnorm);

% generate sample trajectories
initial_set_center_est = rd3.xi0_mean;

rd_est = RefData;
[Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est] = generateRefData(f, initial_set_center_est);
rd_est.directInit(Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est, false);


tmp_xiref = rd_est.getCellArrs();
tmp_xiref_subsampled = {};
num_subsamples = frames;
for n = 1:rd_est.N
    idx = 1:size(tmp_xiref{n}(1:2, :), 2);
    didx = linspace(min(idx), max(idx), num_subsamples);
    res = zeros(size(tmp_xiref{n}(1:2, :), 1), num_subsamples);
    for j = 1:size(tmp_xiref{n}(1:2, :), 1)
        tmp_xiref_subsampled{n}(j, :) = interp1(idx, tmp_xiref{n}(j, :), didx);
    end
end


tmp_xiref3 = rd3.getCellArrs();

for n = 1:rd3.N
    plt_objs3.trajectories{n} = plot(1:frames, tmp_xiref3{n}(1, :), 'color', '#3392ff', 'linewidth', linewidth_global_major, 'displayname', '$\xi^{\mathrm{ref}}$', 'color', '#3392ff', 'displayname', '\ x');
    plt_objs3.trajectories{n} = plot(1:frames, tmp_xiref3{n}(2, :), 'color', '#3392ff', 'linewidth', linewidth_global_major, 'displayname', '$\xi^{\mathrm{ref}}$', 'color', '#3392ff', 'displayname', '\ y');
    if n > 1
        set(plt_objs3.trajectories{n}, 'handlevisibility', 'off');
    end
    hold on;
end


leg = legend();
leg.ItemTokenSize(1) = 10;
filename = strcat(output_root_path, plot_id, '.pdf');
print(fig, filename, '-dpdf');
% close(fig);

end