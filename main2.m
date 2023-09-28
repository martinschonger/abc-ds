% Main script to configure and run experiments to compute polynomial dynamical systems f, polynomial
% Lyapunov functions V, and polynomial Barrier certificates B. Simply add elements to the
% `experiments` cell array, obtaining the initial options struct from `fvbsettings`.
%
% To get started, you can run this script without change (to run the `two_obstacle` robot scenario),
% or uncomment any of the other provided sample experiments.
%
%
% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


% cleanup
close('all');
clear all;

global_options = struct;
global_options.seed = 13;
global_options.generator = 'twister';
global_options.output_root_path = 'output/';
global_options.generate_plots = true;

rng(global_options.seed, global_options.generator);


local_seed = 1;
experiments = {};


%%% [BEGIN] Robot: two_obstacle
experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
    'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
    'dataset', 'robot', 'dataset_opts', struct('exp_list', "datasets/robot_recordings/two_obstacle/"), ...
    'enable_regularization', false, ...
    'deg_f', 5, 'deg_V', 2, 'deg_B', 4, 'deg_B_slack', 2, ...
    'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
    'sdpoptions_penbmi', struct('PBM_MAX_ITER', 100, 'PEN_UP', 0.0, 'UM_MAX_ITER', 250, 'PRECISION_2', 1e-9), ...
    'yz', false);
experiments{end}.pre = {};

experiments{end}.pre{end+1} = experiments{end}.options;
experiments{end}.pre{end}.unmatched.restrict_to_convex = 0;
experiments{end}.pre{end}.deg_f = 1;
experiments{end}.pre{end}.enable_barrier = false;
experiments{end}.pre{end}.unmatched.keep_fc = -1;
%%% [END] Robot: two_obstacle


%%% [BEGIN] Robot: box
% experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
%     'epsilon', 1e-4, 'init_Bc_var', false, 'constraint_version', 3, ...
%     'dataset', 'robot', 'dataset_opts', struct('exp_list', "datasets/robot_recordings/box/"), ...
%     'enable_regularization', false, ...
%     'deg_f', 5, 'deg_V', 2, 'deg_B', 4, 'deg_B_slack', 2, ...
%     'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%     'sdpoptions_penbmi', struct('PBM_MAX_ITER', 600, 'PEN_UP', 0.0, 'UM_MAX_ITER', 250, 'PRECISION_2', 1e-9), ...
%     'yz', true);
% experiments{end}.pre = {};
%%% [END] Robot: box


%%% [BEGIN] LASA multi: [1, 3, 5, 11, 14, 15, 19, 22, 24], no obstacles
%%% Angle: 1
%%% C: 3
%%% G: 5
%%% Leaf2: 11
%%% N: 14
%%% P: 15
%%% Sine: 19
%%% S: 22
%%% Worm: 24
% lasa_ids = [1, 3, 5, 11, 14, 15, 19, 22, 24];
% for i = 1:length(lasa_ids)
%     lid = lasa_ids(i);
% 
%     experiments{end+1}.options = fvbsettings('enable_barrier', false, ...
%         'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
%         'dataset', 'lasa', 'dataset_opts', struct('idx', lid, 'obstacle_repr', "ellipse_axis"), ...
%         'enable_regularization', false, ...
%         'deg_f', 6, 'deg_V', 4, 'deg_B', 2, 'deg_B_slack', 2, ...
%         'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%         'sdpoptions_penbmi', struct('PBM_MAX_ITER', 256, 'PRECISION', 1e-5));
%     experiments{end}.pre = {};
% 
%     experiments{end}.pre{end+1} = experiments{end}.options;
%     experiments{end}.pre{end}.unmatched.restrict_to_convex = 0;
%     experiments{end}.pre{end}.deg_f = 1;
%     experiments{end}.pre{end}.enable_barrier = false;
%     experiments{end}.pre{end}.unmatched.keep_fc = -1;
% end
%%% [END] LASA multi: [1, 3, 5, 11, 14, 15, 19, 22, 24]


%%% [BEGIN] LASA 11: Leaf2, elliptic obstacle
% experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
%     'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
%     'dataset', 'lasa', 'dataset_opts', struct('idx', 11, 'obstacle_repr', "ellipse_axis"), ...
%     'enable_regularization', false, ...
%     'deg_f', 5, 'deg_V', 4, 'deg_B', 4, 'deg_B_slack', 2, ...
%     'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%     'sdpoptions_penbmi', struct('PBM_MAX_ITER', 50, 'PEN_UP', 0.0, 'UM_MAX_ITER', 250));
% experiments{end}.pre = {};
%%% [END] LASA 11: Leaf2


%%% [BEGIN] LASA 15: P, elliptic obstacle
% experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
%     'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
%     'dataset', 'lasa', 'dataset_opts', struct('idx', 15, 'obstacle_repr', "ellipse_axis"), ...
%     'enable_regularization', false, ...
%     'deg_f', 4, 'deg_V', 2, 'deg_B', 3, 'deg_B_slack', 2, ...
%     'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%     'sdpoptions_penbmi', struct('PBM_MAX_ITER', 50, 'PEN_UP', 0.0, 'UM_MAX_ITER', 250));
% experiments{end}.pre = {};
% 
% experiments{end}.pre{end+1} = experiments{end}.options;
% experiments{end}.pre{end}.unmatched.restrict_to_convex = 0;
% experiments{end}.pre{end}.deg_f = 1;
% experiments{end}.pre{end}.enable_barrier = false;
% experiments{end}.pre{end}.unmatched.keep_fc = -1;
%%% [END] LASA 15: P


%%% [BEGIN] LASA 22: S, elliptic obstacle
% experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
%     'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
%     'dataset', 'lasa', 'dataset_opts', struct('idx', 22, 'obstacle_repr', "ellipse_axis"), ...
%     'enable_regularization', false, ...
%     'deg_f', 5, 'deg_V', 2, 'deg_B', 3, 'deg_B_slack', 2, ...
%     'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%     'sdpoptions_penbmi', struct('PBM_MAX_ITER', 100, 'PEN_UP', 0.0, 'UM_MAX_ITER', 250));
% experiments{end}.pre = {};
% 
% experiments{end}.pre{end+1} = experiments{end}.options;
% experiments{end}.pre{end}.unmatched.restrict_to_convex = 0;
% experiments{end}.pre{end}.deg_f = 1;
% experiments{end}.pre{end}.enable_barrier = false;
% experiments{end}.pre{end}.unmatched.keep_fc = -1;
%%% [END] LASA 22: S


%%% [BEGIN] LASA 24: Worm, elliptic obstacle
% experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
%     'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
%     'dataset', 'lasa', 'dataset_opts', struct('idx', 24, 'obstacle_repr', "ellipse_axis"), ...
%     'enable_regularization', false, ...
%     'deg_f', 4, 'deg_V', 2, 'deg_B', 3, 'deg_B_slack', 2, ...
%     'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%     'sdpoptions_penbmi', struct('PBM_MAX_ITER', 50, 'PEN_UP', 0.0));
% experiments{end}.pre = {};
% 
% experiments{end}.pre{end+1} = experiments{end}.options;
% experiments{end}.pre{end}.unmatched.restrict_to_convex = 0;
% experiments{end}.pre{end}.deg_f = 1;
% experiments{end}.pre{end}.enable_barrier = false;
% experiments{end}.pre{end}.unmatched.keep_fc = -1;
% 
% experiments{end}.pre{end+1} = experiments{end}.options;
% experiments{end}.pre{end}.enable_barrier = false;
% experiments{end}.pre{end}.sdpoptions_penbmi.PBM_MAX_ITER = 50;
% experiments{end}.pre{end}.unmatched.keep_fc = -1;
% experiments{end}.pre{end}.unmatched.keep_Vc = -1;
%%% [END] LASA 24: Worm


%%% [BEGIN] LASA 19: Sine, u-shaped obstacle
% experiments{end+1}.options = fvbsettings('enable_barrier', true, ...
%     'epsilon', 1e-3, 'init_Bc_var', false, 'constraint_version', 3, ...
%     'dataset', 'lasa', 'dataset_opts', struct('idx', 19, 'obstacle_repr', "poly"), ...
%     'enable_regularization', false, ...
%     'deg_f', 5, 'deg_V', 2, 'deg_B', 4, 'deg_B_slack', 2, ...
%     'enable_extra_constraint', true, 'regularization_factor', 0.01, 'seed', local_seed, ...
%     'sdpoptions_penbmi', struct('PBM_MAX_ITER', 100, 'PEN_UP', 0.0, 'UM_MAX_ITER', 250));
% experiments{end}.pre = {};
%%% [END] LASA 19: Sine


for curr_exp_idx = 1:length(experiments)
    rng_savepoint = rng;
    fprintf("Experiment %u of %u:\n", curr_exp_idx, length(experiments));
    try
        options = experiments{curr_exp_idx}.options;
        
        [result, options] = run_experiment(global_options, options, experiments{curr_exp_idx}.pre);

        % unscaled poly str to file
        M = length(result.f_fh_str_arr);
        poly_str_filename = strcat(global_options.output_root_path, result.timestamp, '_generated_DS_unscaled.txt');
        [poly_str_fid, msg] = fopen(poly_str_filename, 'at');
        assert(poly_str_fid >= 3, msg)
        for m = 1:M
            poly_str = xistr2xystr(result.f_fh_str_arr{m}, "xi1", "x");
            poly_str = xistr2xystr(poly_str, "xi2", "y");
            if m < M
                fprintf(poly_str_fid, '%s\n', poly_str);
            else
                fprintf(poly_str_fid, '%s', poly_str);
            end
        end
        fclose(poly_str_fid);
        
        % scaled poly str to file
        if strcmp(options.dataset, 'robot')
            M = length(result.f_fh_str_arr_scaled);
            poly_str_filename = strcat(global_options.output_root_path, result.timestamp, '_generated_DS.txt');
            [poly_str_fid, msg] = fopen(poly_str_filename, 'at');
            assert(poly_str_fid >= 3, msg)
            for m = 1:M
                poly_str = xistr2xystr(result.f_fh_str_arr_scaled{m}, "xi1", "x");
                poly_str = xistr2xystr(poly_str, "xi2", "y");
                if m < M
                    fprintf(poly_str_fid, '%s\n', poly_str);
                else
                    fprintf(poly_str_fid, '%s', poly_str);
                end
            end
            fclose(poly_str_fid);
        end

        % unscaled V str to file
        M = length(result.V_fh_str_arr);
        poly_str_filename = strcat(global_options.output_root_path, result.timestamp, '_V_unscaled.txt');
        [poly_str_fid, msg] = fopen(poly_str_filename, 'at');
        assert(poly_str_fid >= 3, msg)
        for m = 1:M
            poly_str = xistr2xystr(result.V_fh_str_arr{m}, "xi1", "x");
            poly_str = xistr2xystr(poly_str, "xi2", "y");
            if m < M
                fprintf(poly_str_fid, '%s\n', poly_str);
            else
                fprintf(poly_str_fid, '%s', poly_str);
            end
        end
        fclose(poly_str_fid);
        
        % scaled V str to file
        if strcmp(options.dataset, 'robot')
            M = length(result.V_fh_str_arr_scaled);
            poly_str_filename = strcat(global_options.output_root_path, result.timestamp, '_V.txt');
            [poly_str_fid, msg] = fopen(poly_str_filename, 'at');
            assert(poly_str_fid >= 3, msg)
            for m = 1:M
                poly_str = xistr2xystr(result.V_fh_str_arr_scaled{m}, "xi1", "x");
                poly_str = xistr2xystr(poly_str, "xi2", "y");
                if m < M
                    fprintf(poly_str_fid, '%s\n', poly_str);
                else
                    fprintf(poly_str_fid, '%s', poly_str);
                end
            end
            fclose(poly_str_fid);
        end

        % unscaled B str to file
        M = length(result.B_fh_str_arr);
        poly_str_filename = strcat(global_options.output_root_path, result.timestamp, '_B_unscaled.txt');
        [poly_str_fid, msg] = fopen(poly_str_filename, 'at');
        assert(poly_str_fid >= 3, msg)
        for m = 1:M
            poly_str = xistr2xystr(result.B_fh_str_arr{m}, "xi1", "x");
            poly_str = xistr2xystr(poly_str, "xi2", "y");
            if m < M
                fprintf(poly_str_fid, '%s\n', poly_str);
            else
                fprintf(poly_str_fid, '%s', poly_str);
            end
        end
        fclose(poly_str_fid);
        
        % scaled B str to file
        if strcmp(options.dataset, 'robot')
            M = length(result.B_fh_str_arr_scaled);
            poly_str_filename = strcat(global_options.output_root_path, result.timestamp, '_B.txt');
            [poly_str_fid, msg] = fopen(poly_str_filename, 'at');
            assert(poly_str_fid >= 3, msg)
            for m = 1:M
                poly_str = xistr2xystr(result.B_fh_str_arr_scaled{m}, "xi1", "x");
                poly_str = xistr2xystr(poly_str, "xi2", "y");
                if m < M
                    fprintf(poly_str_fid, '%s\n', poly_str);
                else
                    fprintf(poly_str_fid, '%s', poly_str);
                end
            end
            fclose(poly_str_fid);
        end
    catch ME
        fprintf(2, '[ERROR] %s (More info: %s)\n', ME.identifier, ME.message);
        fprintf(2, '%s\n', getReport(ME, 'extended'));
    
        diary off;
    
        pause(2);
    end
    rng(rng_savepoint);
end