% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [fig] = generate_PFV_trajectory_following(output_root_path, exp_id, sim_id, plot_id, obstacle_poly, obstacle_polygon_cellarr, yz)

arguments
    output_root_path;
    exp_id;
    sim_id;
    plot_id;
    obstacle_poly;
    obstacle_polygon_cellarr;
    yz = false;
end


frames = 200;
linewidth_global = 5;
linewidth_global_major = 10;
linewidth_global_minor = 1.5;

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
[unsafe_set_coefs, unsafe_set_monomials] = file2poly(obstacle_poly, xi);
unsafe_set{1} = dot(unsafe_set_coefs, unsafe_set_monomials);


plot_dim1 = 1;
plot_dim2 = 2;

[fig, ax] = setup_figure(3*3.41275152778, 3*3.41275152778, false);

myVideo = VideoWriter(output_root_path + plot_id, 'MPEG-4');  % open video file
myVideo.FrameRate = 10;  % can adjust this, 5 - 10 works well for me
open(myVideo);

set(ax,'xtick',[],'ytick',[],'ztick',[],'xcolor','black','ycolor','black');
set(ax,'color',[0 0 0]);
fontsize(fig, 40, "points");
box on;


% reference trajectories and equilibrium
scatter([], [], 0.001, 'o', 'sizedata', 600, 'markeredgecolor', 'white', 'markerfacecolor', 'white', 'linewidth', linewidth_global_major, 'displayname', '\ Equilibrium $\xi^*$');
plot_objs = rd.plotLines(plot_dim1, plot_dim2, {'sizedata', 600, 'handlevisibility', 'off', 'markeredgecolor', 'white', 'markerfacecolor', 'white'}, {'linewidth', linewidth_global, 'color', 'none', 'handlevisibility', 'off'});
quiver(100,100,1,0, 'white', 'linewidth', linewidth_global_minor, 'displayname', '\ Dyn. sys. $f(\xi)$');

% axis limits (depending on ref data)
xlim([limits(1), limits(2)]);
ylim([limits(3), limits(4)]);
axis_limits = axis;


% plot DS
streamlines_plt = plot_streamlines_for_f(f, axis_limits, 200, 1);
set(streamlines_plt, 'linewidth', linewidth_global_minor);
set(streamlines_plt, 'color', 'white');
set(streamlines_plt(1), 'displayname', '$f(\xi)$');
set(streamlines_plt(1), 'handlevisibility', 'off');
set(streamlines_plt(2:end), 'handlevisibility', 'off');

resolution = 0.01;
[X, Y] = meshgrid(axis_limits(1)-resolution:resolution:axis_limits(2)+resolution, axis_limits(3)-resolution:resolution:axis_limits(4)+resolution);
XY = [X(:)'; Y(:)'];


% scale obstacles to [0,1]-space
pgon = {};
for pi = 1:length(obstacle_polygon_cellarr)
    pgon_tmp= obstacle_polygon_cellarr{pi};
    pgon_tmp = pgon_tmp ./ rd.state_maxnorm;
    pgon{pi} = polyshape(pgon_tmp(1,:), pgon_tmp(2,:));
end


% plot Barrier zero-level curve
if options.enable_barrier
    % plot initial set
    for p = 1:length(initial_set)
        fgptilde = sdpvar2fun(initial_set{p}, xi);
        Fgptilde = reshape(fgptilde(XY), size(X));
        [cctilde, hhtilde] = contourf(X, Y, Fgptilde, [0, 0], 'linewidth', linewidth_global, 'color', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{X}_{0}$'), 'handlevisibility', 'off');
        hold on;
    end
    fill([10 10 100 100],[10 100 100 10], 'cyan', 'linewidth', linewidth_global, 'edgecolor', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', '\ Initial set $\mathcal{X}_{0}$');
    
    % plot unsafe set
    for m = 1:length(unsafe_set)
        fgm = sdpvar2fun(unsafe_set{m}, xi);
        Fgm = reshape(fgm(XY), size(X));
        [cc, hh] = contourf(X, Y, Fgm, [0, 0], 'linewidth', linewidth_global, 'color', 'green', 'facecolor', 'green', 'facealpha', 0.25, 'displayname', strcat('$\mathcal{X}_{u}$'), 'handlevisibility', 'off');
        hold on;
    end
    fill([10 10 100 100],[10 100 100 10], 'green', 'linewidth', linewidth_global, 'edgecolor', 'green', 'facecolor', 'green', 'facealpha', 0.25, 'displayname', '\ Unsafe set $\mathcal{X}_{u}$');

    for pi = 1:length(pgon)
        plt_pgon1 = plot(pgon{pi}, 'linewidth', linewidth_global, 'edgecolor', 'white', 'facecolor', 'white', 'facealpha', 0.20);
        set(plt_pgon1, 'handlevisibility', 'off');
    end

    Bfun_eval = reshape(B(XY), size(X));

    % plot barrier 0-level set
    contour(X, Y, Bfun_eval, [0, 0], 'linewidth', linewidth_global, 'color', '#d41919', 'displayname', '\ Barrier $B^{-1}(0)$');
    hold on;
end


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

for i = 1:frames

    for n = 1:rd_est.N
        plt_objs.trajectories{n} = plot(tmp_xiref_subsampled{n}(1, 1:i), tmp_xiref_subsampled{n}(2, 1:i), 'color', '#3392ff', 'linewidth', linewidth_global_major, 'displayname', '$\xi^{\mathrm{ref}}$', 'color', '#ff00ff', 'displayname', '\ Target traj. $\xi^{\mathrm{target}}$');
        if n > 1 || i > 1
            set(plt_objs.trajectories{n}, 'handlevisibility', 'off');
        end
        hold on;
    end

    for n = 1:rd3.N
        plt_objs3.trajectories{n} = plot(tmp_xiref3{n}(1, 1:i), tmp_xiref3{n}(2, 1:i), 'color', '#3392ff', 'linewidth', linewidth_global_major, 'displayname', '$\xi^{\mathrm{ref}}$', 'color', '#3392ff', 'displayname', '\ Actual traj. $\xi^{\mathrm{actual}}$');
        if n > 1 || i > 1
            set(plt_objs3.trajectories{n}, 'handlevisibility', 'off');
        end
        hold on;
    end


    pause(0.01)  % Pause and grab frame
    frame = getframe(gcf);  % get frame
    writeVideo(myVideo, frame);
end


close(myVideo);


leg = legend();
leg.ItemTokenSize(1) = 50;
set(leg, 'Interpreter', 'latex');
set(leg, 'TextColor', 'white');
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf, 'Color','k');
filename = strcat(output_root_path, plot_id, '.pdf');
print(fig, filename, '-dpdf');
close(fig);

end