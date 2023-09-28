% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [fig] = generate_sim_plot_lyapunov_v2(output_root_path, exp_id, sim_id, plot_id)

arguments
    output_root_path;
    exp_id;
    sim_id;
    plot_id;
end

[options, result, f, V, B] = read_exp(exp_id);

exp_list = options.dataset_opts.exp_list;
exp_list_cellarr = fileread(exp_list + "reference_trajectories.txt");
exp_list_cellarr = regexp(exp_list_cellarr, '\r\n|\r|\n', 'split');

[Data, Target, indivTrajStartIndices, timestamps] = recorded_trajectories_to_refdata(exp_list_cellarr);
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
xy_range = xy_range * 2.5;
x_lowerlim = x_min - 0.5 * (xy_range - x_range);
x_upperlim = x_max + 0.5 * (xy_range - x_range);
y_lowerlim = y_min - 0.5 * (xy_range - y_range);
y_upperlim = y_max + 0.5 * (xy_range - y_range);
limits = [x_lowerlim, x_upperlim, y_lowerlim, y_upperlim];
limits = [-1.4, 0.8, -0.55, 1.65];


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

[fig, ax] = setup_figure(3.41275152778, 3.41275152778, true);
fontsize(fig, 7, "points");
box on;
xlh = xlabel(strcat('\xi_', int2str(plot_dim2)));
ylh = ylabel(strcat('\xi_', int2str(plot_dim1)));


% reference trajectories and equilibrium
scatter([], [], 0.001, 'o', 'sizedata', 0.001, 'markeredgecolor', 'black', 'markerfacecolor', 'black', 'linewidth', 0.001, 'displayname', 'Equilibrium $\xi^*$');
plot_objs = rd.plotLines(plot_dim1, plot_dim2, {'sizedata', 20, 'handlevisibility', 'off'}, {'linewidth', 0.5, 'color', 'none', 'handlevisibility', 'off'});
quiver(100,100,1,0, 'black', 'linewidth', 0.5, 'displayname', 'Dyn. sys. $f(\xi)$');

% axis limits (depending on ref data)
xlim([limits(1), limits(2)]);
ylim([limits(3), limits(4)]);
axis_limits = axis;


% plot DS
streamlines_plt = plot_streamlines_for_f(f, axis_limits, 200, 2);
% set(streamlines_plt, 'linewidth', 1);
set(streamlines_plt(1), 'displayname', '$f(\xi)$');
set(streamlines_plt(1), 'handlevisibility', 'off');
set(streamlines_plt(2:end), 'handlevisibility', 'off');

resolution = 0.01;
[X, Y] = meshgrid(axis_limits(1)-resolution:resolution:axis_limits(2)+resolution, axis_limits(3)-resolution:resolution:axis_limits(4)+resolution);
XY = [X(:)'; Y(:)'];


% plot Barrier zero-level curve
if options.enable_barrier
    % plot initial set
    for p = 1:length(initial_set)
        fgptilde = sdpvar2fun(initial_set{p}, xi);
        Fgptilde = reshape(fgptilde(XY), size(X));
        [cctilde, hhtilde] = contourf(X, Y, Fgptilde, [0, 0], 'linewidth', 1, 'color', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{X}_{0}$'), 'handlevisibility', 'off');
        hold on;
    end
    fill([10 10 100 100],[10 100 100 10], 'cyan', 'edgecolor', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', 'Initial set $\mathcal{X}_{0}$');
    
    % plot unsafe set
    for m = 1:length(unsafe_set)
        fgm = sdpvar2fun(unsafe_set{m}, xi);
        Fgm = reshape(fgm(XY), size(X));
        [cc, hh] = contourf(X, Y, Fgm, [0, 0], 'linewidth', 1, 'color', 'black', 'facecolor', 'black', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{X}_{u}$'), 'handlevisibility', 'off');
        hold on;
    end
    fill([10 10 100 100],[10 100 100 10], 'black', 'edgecolor', 'black', 'facecolor', 'black', 'facealpha', 0.35, 'displayname', 'Unsafe set $\mathcal{X}_{u}$');
    

    Bfun_eval = reshape(B(XY), size(X));

    % plot certified safe region
    color_tmp = '#4CBB17';
    facealpha_tmp = 0.25;
    [~, plt_saferegion] = contourf(X, Y, -Bfun_eval, [0, 0], 'linewidth', 1, 'color', 'none', 'facecolor', color_tmp, 'facealpha', facealpha_tmp, 'handlevisibility', 'off');
    uistack(plt_saferegion, 'bottom');
    hold on;
    fill([10 10 100 100],[10 100 100 10], 'g', 'edgecolor', 'none', 'facecolor', color_tmp, 'facealpha', facealpha_tmp, 'displayname', 'Certified safe set $\mathcal{X}_s$');
    
    % plot barrier 0-level set
    contour(X, Y, Bfun_eval, [0, 0], 'linewidth', 1, 'color', '#d41919', 'displayname', 'Barrier $B^{-1}(0)$');
    hold on;
end


% actual trajectories based on executing the DS on the robot
sim_path = sim_id;
[Data3, Target3, indivTrajStartIndices3, timestamps3] = recorded_trajectories_to_refdata({sim_path}, 0, "eval");
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
plt_sampled_traj = rd_est.plotLines(plot_dim1, plot_dim2, {'markeredgecolor', 'none', 'markerfacecolor', 'none', 'handlevisibility', 'off'}, {'color', '#ff00ff', 'displayname', 'Target traj. $\xi^{\mathrm{target}}$'});

rd3.plotLines(1, 2, {'markeredgecolor', 'none', 'markerfacecolor', 'none', 'handlevisibility', 'off'}, {'color', '#3392ff', 'linewidth', 1, 'displayname', 'Actual traj. $\xi^{\mathrm{actual}}$'});


leg = findobj(fig, 'Type', 'Legend');
leg.ItemTokenSize(1) = 5;
% flushLegend(leg, ax, 'northeast');  % [OPTIONAL] https://de.mathworks.com/matlabcentral/fileexchange/57962-flush-legend?s_tid=blogs_rc_6
leg.Position = leg.Position + [0.0035, 0, 0, 0];

% rotate 90 deg
view(-90, 90);
set(gca, 'xdir', 'reverse');
set(gca, 'ydir', 'reverse');
xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'1.5', '1.0','0.5','0','-0.5', '-1', '-1.5'})
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
yticklabels({'-1.5', '-1.0','-0.5','0','0.5','1', '1.5'})

set(xlh,'rotation',0,'VerticalAlignment','middle');
xlh.Position(1) = limits(1) + 0.05;
xlh.Position(2) = xlh.Position(2) + 0.08;
ylh.Position(1) = ylh.Position(1) - 0.1;
ylh.Position(2) = limits(4) - 0.03;


% flushLegend(leg, ax, 'northeast');  % [OPTIONAL] https://de.mathworks.com/matlabcentral/fileexchange/57962-flush-legend?s_tid=blogs_rc_6
leg.Position = leg.Position + [0.014, 0.001, 0, 0];


filename = strcat(output_root_path, plot_id, '.pdf');
print(fig, filename, '-dpdf');

close(fig);

end