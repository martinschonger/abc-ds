% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


% Note: You may have to restart MATLAB after running 'generate_plots_video.m' and before running
%       this script. Otherwise, the generated PDF's may be cropped.

% create plots
output_root_path = "output/figures/";

%% simulation LASA poly fit performance, vs LPV
generate_poly_perf_plot(output_root_path, "2023-09-28_112106", "perf_lasa_angle");
generate_poly_perf_plot(output_root_path, "2023-09-28_112205", "perf_lasa_c");
generate_poly_perf_plot(output_root_path, "2023-09-28_112258", "perf_lasa_g");
generate_poly_perf_plot(output_root_path, "2023-09-28_112340", "perf_lasa_leaf2");
generate_poly_perf_plot(output_root_path, "2023-09-28_112528", "perf_lasa_n");
generate_poly_perf_plot(output_root_path, "2023-09-28_112714", "perf_lasa_p");
generate_poly_perf_plot(output_root_path, "2023-09-28_112750", "perf_lasa_sine");
generate_poly_perf_plot(output_root_path, "2023-09-28_112848", "perf_lasa_s");
generate_poly_perf_plot(output_root_path, "2023-09-28_113008", "perf_lasa_worm");


%% sim plots
generate_sim_plot(output_root_path, "2023-09-28_060738", "sim_lasa_leaf");
generate_sim_plot(output_root_path, "2023-09-28_055510", "sim_lasa_p");
generate_sim_plot(output_root_path, "2023-09-28_055810", "sim_lasa_s");
generate_sim_plot(output_root_path, "2023-09-28_060616", "sim_lasa_worm", true);
generate_sim_plot_ushape_v2(output_root_path, "2023-09-28_061826", "abcds_ushape_v2");


%% 2D and 3D overlays
generate_sim_plot_lyapunov_v2(output_root_path, "2023-09-28_025714", "datasets/robot_recordings/two_obstacle/closed_loop/Recording_2023-09-11_19-31-25", "sim_robot_v3");