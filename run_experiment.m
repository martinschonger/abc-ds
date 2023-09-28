% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [result, options] = run_experiment(global_options, options, pre)
% Wrapper to run a single experiment to compute f, V, and B with the specified options.
% 
% Loads and pre-processes the reference data, sets up the scene with initial and unsafe sets,
% handles any optional steps of the init pipeline, and post-processes the obtained solution
% including initial plotting.
%
% :param global_options: Options struct for the wrapper code. Supported options are
%   `output_root_path` (where to put any produced artifacts like logs and figures), `generate_plots`
%   (boolean flag whether to produce figures or not), `image_resolution` (resolution of the produced
%   figures), `generator` (the random number generator to use).
% :param options: Options struct for the core optimization problem (solved in :func:`fvb`). See
%   :func:`fvbsettings` for a subset of supported options.
% :param pre: Optional init pipeline steps.
% :returns: [result, options], where result contains relevant data computed in the main optimization
%   problem, as well as f, V, and B as strings; and options is the provided options struct extended
%   by any potential default values.

yalmip('clear');

global_options = add_default_option(global_options, 'image_resolution', 600);
output_root_path = global_options.output_root_path;

timestamp = string(datetime('now', 'format', 'yyyy-MM-dd_HHmmss'));
diary_filename = strcat(output_root_path, timestamp, '_log.txt');
diary(diary_filename);

rng(options.seed, global_options.generator);

result = [];


%% load reference trajectories
switch options.dataset
    case 'lasa'
        lasa_idx = options.dataset_opts.idx;
        rd = RefData;
        rd.loadLasa(lasa_idx);
        
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
                unsafe_set_center = [-0.21; 0.37];
                unsafe_set_radius = 0.175;
                ellipse_center_position = [-0.27, -0.35];
                ellipse_axes_length = [0.5181*2, 0.1203*2];
                ellipse_orientation = 99 / 180 * pi;
                initial_set_radius = 0.15;
            case 19
                unsafe_set_center = [-0.43; 0.14];
                unsafe_set_radius = 0.1;
                initial_set_radius = 0.07;
            case 22
                ellipse_center_position = [0.4054, 0.5528];
                ellipse_axes_length = [0.1978*2*1.5, 0.0877*2];
                ellipse_orientation = -5 / 180 * pi;
            case 24
                unsafe_set_center = [-0.5; 0.07];
                unsafe_set_radius = 0.1;
                ellipse_center_position = [-0.495, 0];
                ellipse_axes_length = [0.07*2, 0.2*2];
                ellipse_orientation = 5 / 180 * pi;
            otherwise
                warning('Unknown unsafe set parameters for LASA %s.', num2str(lasa_idx));
                unsafe_set_center = [1; 1];
                unsafe_set_radius = 1e-4;
                ellipse_center_position = [1,1];
                ellipse_axes_length = [0.01, 0.01];
                ellipse_orientation = 0;
        end
    case 'robot'
        exp_list = options.dataset_opts.exp_list;
        exp_list_cellarr = fileread(exp_list + "reference_trajectories.txt");
        exp_list_cellarr = regexp(exp_list_cellarr, '\r\n|\r|\n', 'split');

        [Data, Target, indivTrajStartIndices, timestamps] = recorded_trajectories_to_refdata(exp_list_cellarr, 100, "record", options.unmatched.yz);
        M = 2;
        shift = Target;
        Data(1:M, :) = Data(1:M, :) - shift;
        Target = [0; 0];
        rd = RefData;
        rd.directInit(Data, Target, indivTrajStartIndices, timestamps, true);

        initial_set_center = rd.xi0_mean;
        initial_set_radius = 0.05;
end


x_min = min(rd.Data(1, :));
x_max = max(rd.Data(1, :));
x_range = x_max - x_min;
y_min = min(rd.Data(2, :));
y_max = max(rd.Data(2, :));
y_range = y_max - y_min;
xy_range = max(x_range, y_range);
x_limit_fact = 1.0;
y_limit_fact = 1.0;
x_lowerlim = x_min - xy_range * x_limit_fact;
x_upperlim = x_max + xy_range * x_limit_fact;
y_lowerlim = y_min - xy_range * y_limit_fact;
y_upperlim = y_max + xy_range * y_limit_fact;
limits = [x_lowerlim, x_upperlim, y_lowerlim, y_upperlim];


%% compute DS
xi = sdpvar(rd.M, 1);

% initial/safe set
r1 = initial_set_radius;
r21 = r1 * r1;
center1 = initial_set_center;
initial_set = {};
initial_set{1} = r21 - sum((xi-center1).^2, 1);

% unsafe set
unsafe_set = {};
if options.enable_barrier
    switch options.dataset
        case 'lasa'
            switch options.dataset_opts.obstacle_repr
                case "circle"
                    r = unsafe_set_radius;
                    r2 = r * r;
                    center = unsafe_set_center;
                    unsafe_set{1} = r2 - sum((xi-center).^2, 1);
                case "ellipse_eq"
                    unsafe_set{1} = 1 - (xi-C)' * A * (xi-C);
                case "ellipse_axis"
                    ell = ellipse(ellipse_center_position, ellipse_axes_length, ellipse_orientation);
                    unsafe_set{1} = ell(xi);
                case "poly"
                    [unsafe_set_coefs, unsafe_set_monomials] = file2poly("datasets/lasa_sine_ushape_unsafe_set_poly.json", xi);  % TODO: pass via option
                    unsafe_set{1} = dot(unsafe_set_coefs, unsafe_set_monomials);
            end
        case 'robot'
            [unsafe_set_coefs, unsafe_set_monomials] = file2poly(options.dataset_opts.exp_list + "unsafe_set_poly.json", xi);
            unsafe_set{1} = dot(unsafe_set_coefs, unsafe_set_monomials);
    end
end

% constrain obstacle definition to workspace
% workspace = [[x_min-0.3; y_min-0.3], [x_max+0.3; y_max+0.3]];
% unsafe_set{2} = xi(1) - workspace(1,1);
% unsafe_set{3} = workspace(1,2) - xi(1);
% unsafe_set{4} = xi(2) - workspace(2,1);
% unsafe_set{5} = workspace(2,2) - xi(2);


restrict_to_convex = 1;


history = {};

for pre_idx = 1:length(pre)
    pre_options = pre{pre_idx};
    pre_restrict_to_convex = restrict_to_convex;
    if isfield(pre_options.unmatched, 'restrict_to_convex') && pre_options.unmatched.restrict_to_convex == 0
        pre_restrict_to_convex = 0;
    end
    [~, ~, ~, ~, ~, pre_debug_output, pre_fc, pre_Vc, pre_Bc] = fvb(rd, pre_restrict_to_convex, xi, initial_set, unsafe_set, pre_options);
    history{end+1} = struct('pre_debug_output', pre_debug_output, 'pre_fc', pre_fc, 'pre_Vc', pre_Vc, 'pre_Bc', pre_Bc);
    if pre_idx < length(pre)
        if isfield(pre_options.unmatched, 'keep_fc') && pre_options.unmatched.keep_fc ~= 0
            pre{pre_idx+1}.unmatched.fc_init = history{end+1+pre_options.unmatched.keep_fc}.pre_fc;
            pre{pre_idx+1}.unmatched.fc_init_monomials = history{end+1+pre_options.unmatched.keep_fc}.pre_debug_output.f_monomials;
            pre{pre_idx+1}.unmatched.fc_init_sdpvar = xi;
        end
        if isfield(pre_options.unmatched, 'keep_Vc') && pre_options.unmatched.keep_Vc ~= 0
            pre{pre_idx+1}.unmatched.Vc_init = history{end+1+pre_options.unmatched.keep_Vc}.pre_Vc;
            pre{pre_idx+1}.unmatched.Vc_init_monomials = history{end+1+pre_options.unmatched.keep_Vc}.pre_debug_output.V_monomials;
            pre{pre_idx+1}.unmatched.Vc_init_sdpvar = xi;
        end
        if isfield(pre_options.unmatched, 'keep_Bc') && pre_options.unmatched.keep_Bc ~= 0
            pre{pre_idx+1}.unmatched.Bc_init = history{end+1+pre_options.unmatched.keep_Bc}.pre_Bc;
            pre{pre_idx+1}.unmatched.Bc_init_monomials = history{end+1+pre_options.unmatched.keep_Bc}.pre_debug_output.B_monomials;
            pre{pre_idx+1}.unmatched.Bc_init_sdpvar = xi;
        end
    else
        if isfield(pre_options.unmatched, 'keep_fc') && pre_options.unmatched.keep_fc ~= 0
            options.unmatched.fc_init = history{end+1+pre_options.unmatched.keep_fc}.pre_fc;
            options.unmatched.fc_init_monomials = history{end+1+pre_options.unmatched.keep_fc}.pre_debug_output.f_monomials;
            options.unmatched.fc_init_sdpvar = xi;
        end
        if isfield(pre_options.unmatched, 'keep_Vc') && pre_options.unmatched.keep_Vc ~= 0
            options.unmatched.Vc_init = history{end+1+pre_options.unmatched.keep_Vc}.pre_Vc;
            options.unmatched.Vc_init_monomials = history{end+1+pre_options.unmatched.keep_Vc}.pre_debug_output.V_monomials;
            options.unmatched.Vc_init_sdpvar = xi;
        end
        if isfield(pre_options.unmatched, 'keep_Bc') && pre_options.unmatched.keep_Bc ~= 0
            options.unmatched.Bc_init = history{end+1+pre_options.unmatched.keep_Bc}.pre_Bc;
            options.unmatched.Bc_init_monomials = history{end+1+pre_options.unmatched.keep_Bc}.pre_debug_output.B_monomials;
            options.unmatched.Bc_init_sdpvar = xi;
        end
    end
    
end

[fhat, V, dVdx, B, dBdx, debug_output, fc, Vc, Bc] = fvb(rd, restrict_to_convex, xi, initial_set, unsafe_set, options);

f = fhat;

result.id = timestamp;
result.mse = debug_output.mse;
result.actual_mse = result.mse * 2;
result.min_q_eigvals = debug_output.min_q_eigvals;
result.primal_residuals = debug_output.primal_residuals;
result.relative_factors = debug_output.relative_factors;
result.comp_time = debug_output.sol.solvertime;
result.fc = fc;
result.Vc = Vc;
result.Bc = Bc;

result.timestamp = timestamp;


[~, f_fh_str_arr] = mvfun2str(f);
result.f_fh_str_arr = f_fh_str_arr;
if strcmp(options.dataset, 'robot')
    f_scaled_back = @(xi) f((xi - shift) ./ rd.state_maxnorm) .* rd.vel_maxnorm;
    [~, f_fh_str_arr_scaled] = mvfun2str(f_scaled_back);
    result.f_fh_str_arr_scaled = f_fh_str_arr_scaled;
end

[~, V_fh_str_arr] = mvfun2str(V);
result.V_fh_str_arr = V_fh_str_arr;
if strcmp(options.dataset, 'robot')
    V_scaled_back = @(xi) V((xi - shift) ./ rd.state_maxnorm);
    [~, V_fh_str_arr_scaled] = mvfun2str(V_scaled_back);
    result.V_fh_str_arr_scaled = V_fh_str_arr_scaled;
end

[~, B_fh_str_arr] = mvfun2str(B);
result.B_fh_str_arr = B_fh_str_arr;
if strcmp(options.dataset, 'robot')
    B_scaled_back = @(xi) B((xi - shift) ./ rd.state_maxnorm);
    [~, B_fh_str_arr_scaled] = mvfun2str(B_scaled_back);
    result.B_fh_str_arr_scaled = B_fh_str_arr_scaled;
end


% result_json = struct;
result_json = result;
result_json_filename = strcat(output_root_path, timestamp, '_result.json');
result_json_content = jsonencode(result_json);
result_json_fid = fopen(result_json_filename, 'w');
fprintf(result_json_fid , '%s', result_json_content);
fclose(result_json_fid);


options_log_filename = strcat(output_root_path, timestamp, '_options.json');
options_copy = options;
if isfield(options_copy, 'unmatched')
    options_copy = rmfield(options_copy, 'unmatched'); % FIXME
end
options_log_content = jsonencode(options_copy);
options_log_fid = fopen(options_log_filename, 'w');
fprintf(options_log_fid , '%s', options_log_content);
fclose(options_log_fid);


if global_options.generate_plots
    %% plotting
    % plot setup
    plot_dim1 = 1;
    plot_dim2 = 2;
    
    fig = setup_figure();
    
    xlabel(strcat('\xi_', int2str(plot_dim1)));
    ylabel(strcat('\xi_', int2str(plot_dim2)));
    
    % reference trajectories and equilibrium
    plot_objs = rd.plotLines(plot_dim1, plot_dim2);
    
    % axis limits (depending on ref data)
    xlim([limits(1), limits(2)]);
    ylim([limits(3), limits(4)]);
    axis_limits = axis;
    
    % legend config
    leg = findobj(gcf, 'Type', 'Legend');
    leg.ItemTokenSize(1) = 15;
    
    
    % plot DS
    streamlines_plt = plot_streamlines_for_f(f, axis_limits);
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
    
    
    % plot Barrier zero-level curve
    if options.enable_barrier
        Bfun_eval = reshape(B(XY), size(X));
        contour(X, Y, Bfun_eval, [0, 0], 'linewidth', 1, 'color', '#7D0675', 'displayname', '$B(\xi)=0$');
        hold on;
        
        % plot initial set
        for p = 1:length(initial_set)
            fgptilde = sdpvar2fun(initial_set{p}, xi);
            Fgptilde = reshape(fgptilde(XY), size(X));
            [cctilde, hhtilde] = contourf(X, Y, Fgptilde, [0, 0], 'linewidth', 1, 'color', 'cyan', 'facecolor', 'cyan', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{D}_{0;', int2str(p), '}$'));
            hold on;
        end
    
        % plot unsafe set
        for m = 1:length(unsafe_set)
            fgm = sdpvar2fun(unsafe_set{m}, xi);
            Fgm = reshape(fgm(XY), size(X));
            [cc, hh] = contourf(X, Y, Fgm, [0, 0], 'linewidth', 1, 'color', 'black', 'facecolor', 'black', 'facealpha', 0.35, 'displayname', strcat('$\mathcal{D}_{u;', int2str(m), '}$'));
            hold on;
        end
    end
    
    
    % sample trajectories from estimated DS
    initial_set_center_est = rd.xi0_mean;
    % initial_set_center_est = rd.Data(1:rd.M, rd.indivTrajStartIndices(1:end-1));

    rd_est = RefData;
    [Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est] = generateRefData(f, initial_set_center_est);
    rd_est.directInit(Data_est, Target_est, indivTrajStartIndices_est, Timestamps_est, false);
    plt_sampled_traj = rd_est.plotLines(plot_dim1, plot_dim2, {'markeredgecolor', 'none', 'markerfacecolor', 'none', 'handlevisibility', 'off'}, {'color', '#ff00ff'});
    
    
    cbar = colorbar('ticks', plt_V.LevelList, 'ticklabels', plt_V.LevelList);
    yt = get(cbar, 'ytick');
    set(cbar, 'yticklabel', sprintf('%.2g\n', yt));
    set(cbar, 'ticklabelinterpreter', 'latex');

    % save figure
    filename = strcat(output_root_path, timestamp, '.png'); % TODO: change to '.pdf'
    exportgraphics(fig, filename, 'contenttype', 'image', 'resolution', global_options.image_resolution); % TODO: change to 'vector'
    
    
    %% plot V
    filename = strcat(output_root_path, timestamp, '_V.png'); % TODO: change to '.pdf'
    figV = plot_V(plot_dim1, plot_dim2, X, Y, Vfun_eval, axis_limits, streamlines_plt, filename, global_options.image_resolution);
    close(figV);

    
    %% plot B
    if options.enable_barrier
        filename = strcat(output_root_path, timestamp, '_B.png'); % TODO: change to '.pdf'
        figB = plot_B(plot_dim1, plot_dim2, X, Y, Bfun_eval, axis_limits, streamlines_plt, filename, global_options.image_resolution);
        close(figB);
    end
    
    
    %%
    close(fig);
end


%%
diary off;

sources_filename = strcat(output_root_path, timestamp, '_sources.txt');
[sources_logfile_fid, msg] = fopen(sources_filename, 'at');
assert(sources_logfile_fid >= 3, msg)
append_to_logfile(sources_logfile_fid, {'main2.m', 'run_experiment.m', 'fvb.m'});
fclose(sources_logfile_fid);


end