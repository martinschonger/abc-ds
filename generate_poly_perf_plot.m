% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [fig] = generate_poly_perf_plot(output_root_path, exp_id, plot_id)

[options, result, f, V, B] = read_exp(exp_id);
lasa_idx = options.dataset_opts.idx;
rd = RefData;
rd.loadLasa(lasa_idx);


x_min = min(rd.Data(1, :));
x_max = max(rd.Data(1, :));
x_range = x_max - x_min;
y_min = min(rd.Data(2, :));
y_max = max(rd.Data(2, :));
y_range = y_max - y_min;
xy_range = max(x_range, y_range);
xy_range = xy_range * 1.8;
x_lowerlim = x_min - 0.5 * (xy_range - x_range);
x_upperlim = x_max + 0.5 * (xy_range - x_range);
y_lowerlim = y_min - 0.5 * (xy_range - y_range);
y_upperlim = y_max + 0.5 * (xy_range - y_range);
limits = [x_lowerlim, x_upperlim, y_lowerlim, y_upperlim];


plot_dim1 = 1;
plot_dim2 = 2;

[fig, ax] = setup_figure(1.058843685114173, 1.058843685114173, false);
set(ax,'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);


% reference trajectories and equilibrium
plot_objs = rd.plotLines(plot_dim1, plot_dim2, {'sizedata', 20}, {'linewidth', 0.5});

% axis limits (depending on ref data)
xlim([limits(1), limits(2)]);
ylim([limits(3), limits(4)]);
axis_limits = axis;


% plot DS
streamlines_plt = plot_streamlines_for_f(f, axis_limits, 200, 1);
set(streamlines_plt(1), 'displayname', '$f(\xi)$');
set(streamlines_plt(2:end), 'handlevisibility', 'off');

resolution = 0.01;
[X, Y] = meshgrid(axis_limits(1)-resolution:resolution:axis_limits(2)+resolution, axis_limits(3)-resolution:resolution:axis_limits(4)+resolution);
XY = [X(:)'; Y(:)'];

% plot Lyapunov function
Vfun_eval = reshape(V(XY), size(X));
[C, plt_V] = contourf(X, Y, Vfun_eval, 'linecolor', 'white', 'linewidth', 0.5, 'displayname', '$V(\xi)$');
cmap = crameri('lajolla');
cmap = cmap(1:160, :);
cmap = cmap(round(linspace(1, 160, 10)), :);
colormap(cmap);
uistack(plt_V, 'bottom');

% generate sample trajectories
initial_set_center_est = rd.xi0_mean;

rd_est = RefData;
[Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est] = generateRefData(f, initial_set_center_est);
rd_est.directInit(Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est, false);
plt_sampled_traj = rd_est.plotLines(plot_dim1, plot_dim2, {'markeredgecolor', 'none', 'markerfacecolor', 'none', 'handlevisibility', 'off'}, {'color', '#ff00ff'});


filename = strcat(output_root_path, plot_id, '.pdf');
print(fig, filename, '-dpdf');

close(fig);

end