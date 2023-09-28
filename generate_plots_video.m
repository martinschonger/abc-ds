% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


% create plots
output_root_path = "output/video/";


%% 2-obstacle experiment
%%
P1 = [[-0.138; -0.07], [-0.23; -0.015], [-0.094; 0.21], [-0.005; 0.155]];
P2 = [[-0.273; 0.299], [-0.309; 0.344], [-0.103; 0.502], [-0.066; 0.46]];
% shift obstacles
shift_dir1 = [-0.23; -0.015] - [-0.094; 0.21];
shift_dir1 = shift_dir1 ./ norm(shift_dir1);
shift_dir1b = [-0.23; -0.015] - [-0.138; -0.07];
shift_dir1b = shift_dir1b ./ norm(shift_dir1b);
shift_dir2 = [-0.103; 0.502] - [-0.309; 0.344];
shift_dir2 = shift_dir2 ./ norm(shift_dir2);
shift_dir2b = [-0.309; 0.344] - [-0.273; 0.299];
shift_dir2b = shift_dir2b ./ norm(shift_dir2b);
P1 = P1 + (shift_dir1 * 0.04);
P1 = P1 + (shift_dir1b * 0.02);
P2 = P2 + (shift_dir2 * 0.04);
P2 = P2 + (shift_dir2b * 0.02);
obstacle_polygon_cellarr = {};
obstacle_polygon_cellarr{1} = P1;
obstacle_polygon_cellarr{2} = P2;

generate_PFV_recording_large(output_root_path, "2023-09-28_025714", "2_recording_large_1", "datasets/robot_recordings/two_obstacle/unsafe_set_poly.json", obstacle_polygon_cellarr, 6);

generate_PFV_recording_multi(output_root_path, "2023-09-28_025714", "2_recording_multi_1", "datasets/robot_recordings/two_obstacle/unsafe_set_poly.json", obstacle_polygon_cellarr);

%%
generate_PFV_trajectory_following(output_root_path, "2023-09-28_025714", "datasets/robot_recordings/two_obstacle/closed_loop_improved/Recording_2023-09-20_14-47-15", "3_clc_1", "datasets/robot_recordings/two_obstacle/unsafe_set_poly.json", obstacle_polygon_cellarr);

%%
% generate_PFV_interaction_helperplot(output_root_path, "2023-09-28_025714", "datasets/robot_recordings/two_obstacle/closed_loop_with_perturbation/Recording_2023-09-20_15-42-02", "4_interaction_helper_1");

interaction_timepoints = [55, 80; 159, 198; 276, 326; 390, 423; 457, 544];  % for 1000 frames

generate_PFV_interaction(output_root_path, "2023-09-28_025714", "datasets/robot_recordings/two_obstacle/closed_loop_with_perturbation/Recording_2023-09-20_15-42-02", "4_interaction_1", "datasets/robot_recordings/two_obstacle/unsafe_set_poly.json", obstacle_polygon_cellarr, interaction_timepoints, 675);



%% box experiment
%%
P2 = [[14.75; 0], [14.75; 27.9], [13.4; 27.9], [13.4; 0]];
P3 = [[14.75; 0], [14.75; 2.3], [-14.75; 2.3], [-14.75; 0]];
P4 = [[-13.4; 0], [-13.4; 27.9], [-14.75; 27.9], [-14.75; 0]];

offset_vect_x3 = [5; 0];
P2 = P2 + offset_vect_x3;
offset_vect_x4 = [1, 1, -1, -1] .* offset_vect_x3;
P3 = P3 + offset_vect_x4;
P4 = P4 - offset_vect_x3;


translation_vect = [38; -5.2];  % large box
P2 = P2 + translation_vect;
P3 = P3 + translation_vect;
P4 = P4 + translation_vect;

% convert to meters
P2 = P2 ./ 100;
P3 = P3 ./ 100;
P4 = P4 ./ 100;

obstacle_expansion_fact = 0.0; % meters
offset_vect_x = [1, 1, -1, -1] * obstacle_expansion_fact;
offset_vect_y = [-1, 0, 0, -1] * 0.0;
offset_vect = [offset_vect_x; offset_vect_y];

P2 = P2 + offset_vect;
P3 = P3 + offset_vect;
P4 = P4 + offset_vect;

P5 = [P2(:, 1:3), [P2(1, 4); P2(2, 4) + 0.023], [P4(1, 1); P4(2, 1) + 0.023], P4(:, 2:4)];

obstacle_polygon_cellarr = {};
obstacle_polygon_cellarr{1} = P5;


generate_PFV_recording_large(output_root_path, "2023-09-28_033613", "2_recording_large_2", "datasets/robot_recordings/box/unsafe_set_poly.json", obstacle_polygon_cellarr, 5, true);

generate_PFV_recording_multi(output_root_path, "2023-09-28_033613", "2_recording_multi_2", "datasets/robot_recordings/box/unsafe_set_poly.json", obstacle_polygon_cellarr, true);

%%
generate_PFV_trajectory_following(output_root_path, "2023-09-28_033613", "datasets/robot_recordings/box/closed_loop/Recording_2023-09-22_18-32-12", "3_clc_2", "datasets/robot_recordings/box/unsafe_set_poly.json", obstacle_polygon_cellarr, true);

%%
% generate_PFV_interaction_helperplot(output_root_path, "2023-09-28_033613", "datasets/robot_recordings/box/closed_loop_with_perturbation/Recording_2023-09-22_19-35-42", "4_interaction_helper_2", true);

interaction_timepoints = [85, 120; 208, 294; 444, 477; 689, 770];  % for 1000 frames

generate_PFV_interaction(output_root_path, "2023-09-28_033613", "datasets/robot_recordings/box/closed_loop_with_perturbation/Recording_2023-09-22_19-35-42", "4_interaction_2", "datasets/robot_recordings/box/unsafe_set_poly.json", obstacle_polygon_cellarr, interaction_timepoints, 0, true);