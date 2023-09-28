% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [fig] = generate_sim_plot_ushape_v2(output_root_path, exp_id, plot_id)

arguments
    output_root_path;
    exp_id;
    plot_id;
end

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
xy_range = xy_range * 2.5;
x_lowerlim = x_min - 0.5 * (xy_range - x_range);
x_upperlim = x_max + 0.5 * (xy_range - x_range);
y_lowerlim = y_min - 0.6 * (xy_range*0.68292518739 - y_range);
y_upperlim = y_max + 0.4 * (xy_range*0.68292518739 - y_range);
limits = [x_lowerlim, x_upperlim, y_lowerlim, y_upperlim];

initial_set_center = rd.xi0_mean;
initial_set_radius = 0.1;

switch lasa_idx
    case 1
        unsafe_set_center = [-0.51; 0.3];
        unsafe_set_radius = 0.15;
    case 3
        unsafe_set_center = [-0.31; 0.33];
        unsafe_set_radius = 0.15;
        ellipse_center_position = [-0.1234, 0.3208];
        ellipse_axes_length = [0.1746*2*2, 0.0824*2];
        ellipse_orientation = -5 / 180 * pi;
    case 4
        unsafe_set_center = [-0.81; -0.58];
        unsafe_set_radius = 0.125;
    case 5
        unsafe_set_center = [-0.2; 0.25];
        unsafe_set_radius = 0.1;
    case 11
        ellipse_center_position = [-0.3925, 0.0253];
        ellipse_axes_length = [0.3593*2, 0.1441*2];
        ellipse_orientation = 78 / 180 * pi;
    case 14
        unsafe_set_center = [-0.2; -0.09];
        unsafe_set_radius = 0.075;
    case 15
        limits = [-1.62, 1.62, -1.63, 1.63];
        unsafe_set_center = [-0.21; 0.37];
        unsafe_set_radius = 0.175;
        ellipse_center_position = [-0.2631, -0.3640];
        ellipse_axes_length = [0.5181*2, 0.1203*2];
        ellipse_orientation = 99 / 180 * pi;
        initial_set_radius = 0.15;
    case 19
        unsafe_set_center = [-0.43; 0.14];
        unsafe_set_radius = 0.1;
        initial_set_radius = 0.07;
    case 22
        limits = [-0.62, 1.12, -0.63, 1.13];
        ellipse_center_position = [0.4054, 0.5528];
        ellipse_axes_length = [0.1978*2*1.5, 0.0877*2];
        ellipse_orientation = -5 / 180 * pi;
    case 24
        limits = [-1.62, 0.62, -1.13, 1.13];
        unsafe_set_center = [-0.5; 0.07];
        unsafe_set_radius = 0.1;
        ellipse_center_position = [-0.495, 0];
        ellipse_axes_length = [0.07*2, 0.2*2];
        ellipse_orientation = 5 / 180 * pi;
    case 25
        ellipse_center_position = [-0.255, 0.09];
        ellipse_axes_length = [0.4*2, 0.08*2];
        ellipse_orientation = 96 / 180 * pi;
    otherwise
        warning('Unknown unsafe set parameters for LASA %s.', num2str(lasa_idx));
        unsafe_set_center = [1; 1];
        unsafe_set_radius = 1e-4;
        ellipse_center_position = [1,1];
        ellipse_axes_length = [0.01, 0.01];
        ellipse_orientation = 0;
end


xi = sdpvar(rd.M, 1);

% initial/safe set
r1 = initial_set_radius;
r21 = r1 * r1;
center1 = initial_set_center;
initial_set = {};
initial_set{1} = r21 - sum((xi-center1).^2, 1);

% unsafe set
unsafe_set = {};
[unsafe_set_coefs, unsafe_set_monomials] = file2poly('datasets/lasa_sine_ushape_unsafe_set_poly.json', xi);
unsafe_set{1} = dot(unsafe_set_coefs, unsafe_set_monomials);


plot_dim1 = 1;
plot_dim2 = 2;

[fig, ax] = setup_figure(3.41275152778, 0.68292518739*3.41275152778, true, true);
set(ax,'xtick',[],'ytick',[],'ztick',[],'xcolor','black','ycolor','black');
box on;


% reference trajectories and equilibrium
scatter([], [], 0.001, 'o', 'sizedata', 0.001, 'markeredgecolor', 'black', 'markerfacecolor', 'black', 'linewidth', 0.001, 'displayname', '$\xi^*$');
plot_objs = rd.plotLines(plot_dim1, plot_dim2, {'sizedata', 20, 'handlevisibility', 'off'}, {'linewidth', 0.5});
quiver(100,100,1,0, 'black', 'linewidth', 0.5, 'displayname', '$f(\xi)$');

% axis limits (depending on ref data)
xlim([limits(1), limits(2)]);
ylim([limits(3), limits(4)]);
axis_limits = axis;


% plot DS
streamlines_plt = plot_streamlines_for_f(f, axis_limits, 200, 2);
set(streamlines_plt(1), 'displayname', '$f(\xi)$');
set(streamlines_plt(1), 'handlevisibility', 'off');
set(streamlines_plt(2:end), 'handlevisibility', 'off');

resolution = 0.01;
[X, Y] = meshgrid(axis_limits(1)-resolution:resolution:axis_limits(2)+resolution, axis_limits(3)-resolution:resolution:axis_limits(4)+resolution);
XY = [X(:)'; Y(:)'];


if options.enable_barrier
    % plot initial set
    for p = 1:length(initial_set)
        fgptilde = sdpvar2fun(initial_set{p}, xi);
        Fgptilde = reshape(fgptilde(XY), size(X));
        [cctilde, hhtilde] = contourf(X, Y, Fgptilde, [0, 0], 'linewidth', 1, 'color', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{X}_{0}$'), 'handlevisibility', 'off');
        hold on;
    end
    fill([10 10 100 100],[10 100 100 10], 'cyan', 'edgecolor', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', '$\mathcal{X}_{0}$');
    
    % plot unsafe set
    for m = 1:length(unsafe_set)
        fgm = sdpvar2fun(unsafe_set{m}, xi);
        Fgm = reshape(fgm(XY), size(X));
        [cc, hh] = contourf(X, Y, Fgm, [0, 0], 'linewidth', 1, 'color', 'black', 'facecolor', 'black', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{X}_{u}$'), 'handlevisibility', 'off');
        hold on;
    end
    fill([10 10 100 100],[10 100 100 10], 'black', 'edgecolor', 'black', 'facecolor', 'black', 'facealpha', 0.35, 'displayname', '$\mathcal{X}_{u}$');
    

    Bfun_eval = reshape(B(XY), size(X));

    % plot barrier 0-level set
    contour(X, Y, Bfun_eval, [0, 0], 'linewidth', 1, 'color', '#d41919', 'displayname', '$B^{-1}(0)$');
    hold on;
end

% generate sample trajectories
initial_set_center_est = rd.Data(1:rd.M, rd.indivTrajStartIndices(1:end-1));

rd_est = RefData;
[Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est] = generateRefData(f, initial_set_center_est);
rd_est.directInit(Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est, false);
plt_sampled_traj = rd_est.plotLines(plot_dim1, plot_dim2, {'markeredgecolor', 'none', 'markerfacecolor', 'none', 'handlevisibility', 'off'}, {'color', '#ff00ff', 'displayname', '$\xi^{\mathrm{sim}}$'});


leg = findobj(fig, 'Type', 'Legend');
leg.ItemTokenSize(1) = 5;
% flushLegend(leg, ax, 'northwest');  % [OPTIONAL] https://de.mathworks.com/matlabcentral/fileexchange/57962-flush-legend?s_tid=blogs_rc_6
leg.Position = leg.Position + [+0.0005, -0.001, 0, 0];


filename = strcat(output_root_path, plot_id, '.pdf');
print(fig, filename, '-dpdf', '-r2400');

close(fig);

end