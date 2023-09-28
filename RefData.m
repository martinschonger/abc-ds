% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


classdef RefData < handle
    %DEMONSTRATIONS Uniform interface to reference data
    %   N               number of trajectories
    %   T_n             number of samples of nth trajectory
    %   T=T_1+...+T_N   total number of samples

    properties
        Target  % equilibrium/attractor
        Data  % stacked x and xdot, concatenated trajectories; shape (2*M,T)
        M  % dimension of state space D
        N  % number of trajectories
        indivTrajStartIndices % indices in second dim of Data where the different ref. trajs. start
        Timestamps

        % DEBUG:
        state_maxnorm
        vel_maxnorm
        shift
    end

    methods
        % load from file:
        % loadedObj = load(filename).origObjName;

        function T = T(obj)
            T = size(obj.Data, 2);
        end

        function directInit(obj, Data, Target, indivTrajStartIndices, timestamps, scale_me, scale_fact_state, scale_fact_vels)
            arguments
                obj;
                Data;
                Target;
                indivTrajStartIndices;
                timestamps;
                scale_me;
                scale_fact_state = 0.0;
                scale_fact_vels = 0.0;
            end

            obj.Data = Data;
            obj.M = size(obj.Data, 1) / 2;
            obj.Target = Target;
            obj.indivTrajStartIndices = indivTrajStartIndices;
            obj.N = length(obj.indivTrajStartIndices) - 1;
            obj.Timestamps = timestamps;

            if scale_me
                % [Optional] Scale workspace and velocities to range/extent [0, 1]:
                if scale_fact_state == 0.0
                    obj.state_maxnorm = max(sqrt(sum(obj.Data(1:obj.M, :).^2, 1)));
                    obj.vel_maxnorm = max(sqrt(sum(obj.Data(obj.M+1:end, :).^2, 1)));
                else
                    obj.state_maxnorm = scale_fact_state;
                    obj.vel_maxnorm = scale_fact_vels;
                end
                tmp_fact = [repmat(obj.state_maxnorm, obj.M, 1); repmat(obj.vel_maxnorm, obj.M, 1)];
                obj.Data = obj.Data ./ tmp_fact;
            end
        end

        function loadCustom(obj, filepath)
            S = load(filepath);
            obj.Data = S.Data_sh;
            obj.M = size(obj.Data, 1) / 2;
            obj.Target = zeros(obj.M, 1);
            obj.N = length(S.data);
            indivTrajStartIndices_tmp = [1];
            for n = 1:obj.N
                indivTrajStartIndices_tmp(end+1) = length(S.data{n});
            end
            obj.indivTrajStartIndices = cumsum(indivTrajStartIndices_tmp);
        end

        function loadLasa(obj, idx)
            sub_sample = 10;  % take 100 samples per trajectory
            nb_trajectories = 7;  % take all 7 trajectories
            [Data, Data_sh, att, x0_all, data, dt] = load_LASA_dataset_shape_DS(idx, sub_sample, nb_trajectories);
            obj.Data = Data_sh;
            obj.M = size(obj.Data, 1) / 2;

            % [Optional] Scale workspace and velocities to range/extent [0, 1]:
            obj.state_maxnorm = max(sqrt(sum(obj.Data(1:obj.M, :).^2, 1)));
            obj.vel_maxnorm = max(sqrt(sum(obj.Data(obj.M+1:end, :).^2, 1)));
            tmp_fact = [repmat(obj.state_maxnorm, obj.M, 1); repmat(obj.vel_maxnorm, obj.M, 1)];
            obj.Data = obj.Data ./ tmp_fact;

            obj.Target = zeros(obj.M, 1);
            obj.N = length(data);
            indivTrajStartIndices_tmp = [1];
            for n = 1:obj.N
                indivTrajStartIndices_tmp(end+1) = length(data{n});
            end
            obj.indivTrajStartIndices = cumsum(indivTrajStartIndices_tmp);

            warning('The Data and data fields of LASA datasets may not be aligned.');
        end

        function res = xi0_mean(obj)
            res = mean(obj.Data(1:obj.M, obj.indivTrajStartIndices(1:end-1)), 2);
        end

        function plt_objs = plot(obj, dim1, dim2)
            arguments
                obj RefData;
                dim1(1, 1) {mustBeInteger} = 1;
                dim2(1, 1) {mustBeInteger} = 2;
            end
            plt_objs.trajectories = plot(obj.Data(dim1, :), obj.Data(dim2, :), '.', 'color', '#3392ff', 'markersize', 10, 'displayname', '$\xi^{\mathrm{ref}}$');
            hold on;
            plt_objs.equilibrium = scatter(obj.Target(dim1), obj.Target(dim2), 100, 'black', 'x', 'linewidth', 1, 'displayname', '$\xi^*$');
            hold on;
        end

        function plt_objs = plotLines(obj, dim1, dim2, pltopts_equi, pltopts_traj, plot_velocities, plot_lines)
            arguments
                obj RefData;
                dim1(1, 1) {mustBeInteger} = 1;
                dim2(1, 1) {mustBeInteger} = 2;
                pltopts_equi = {};
                pltopts_traj = {};
                plot_velocities = false;
                plot_lines = true;
            end
            plt_objs.equilibrium = scatter(obj.Target(dim1), obj.Target(dim2), 200, 'o', 'markeredgecolor', 'black', 'markerfacecolor', 'black', 'linewidth', 1, 'displayname', '$\xi^*$', pltopts_equi{:});
            hold on;
            if plot_lines
                tmp_xiref = obj.getCellArrs();
                for n = 1:obj.N
                    plt_objs.trajectories{n} = plot(tmp_xiref{n}(dim1, :), tmp_xiref{n}(dim2, :), 'color', '#3392ff', 'linewidth', 1, 'displayname', '$\xi^{\mathrm{ref}}$', pltopts_traj{:});
                    if n > 1
                        set(plt_objs.trajectories{n}, 'handlevisibility', 'off');
                    end
                    hold on;
                end
            end
            if plot_velocities
                plt_objs.velocities = quiver(obj.Data(dim1, :), obj.Data(dim2, :), obj.Data(obj.M+dim1, :), obj.Data(obj.M+dim2, :), 'color', 'black', 'linewidth', 0.5, 'displayname', '$\dot{\xi}^{\mathrm{ref}}$');
            end
        end

        function [plt_objs, abs_vels] = plot2DVelocity(obj, dim1, dim2, pltopts)
            arguments
                obj RefData;
                dim1(1, 1) {mustBeInteger} = 1;
                dim2(1, 1) {mustBeInteger} = 2;
                pltopts = {};
            end

            abs_vels = {};
            [~, vels] = obj.getCellArrs;
            for n = 1:obj.N
                abs_vels{end+1} = sqrt(sum(vels{n}.^2, 1));
                plot(obj.Timestamps{n}, abs_vels{n}, pltopts{:});
                hold on;
            end
        end

        function [xiref, xiref_dot] = getCellArrs(obj)
            xiref = {};
            xiref_dot = {};
            for n = 1:obj.N
                xiref{end+1} = obj.Data(1:obj.M, obj.indivTrajStartIndices(n):obj.indivTrajStartIndices(n+1)-1);
                xiref_dot{end+1} = obj.Data(obj.M+1:end, obj.indivTrajStartIndices(n):obj.indivTrajStartIndices(n+1)-1);
            end
        end

        function res = mean_distance(obj, other, dist_fun)
            assert(obj.N == other.N, 'Expected same number of trajectories.');
            res = 0;
            for i = 1:obj.N
                [xiref_obj, ~] = obj.getCellArrs();
                [xiref_other, ~] = other.getCellArrs();
                res = res + dist_fun(xiref_obj, xiref_other);
            end
            res = res / obj.N;
        end

        % function s = saveobj(obj)
        %     s.Target = obj.Target;
        %     s.Data = obj.Data;
        %     s.M = obj.M;
        %     s.N = obj.N;
        %     s.indivTrajStartIndices = obj.indivTrajStartIndices;
        % end
    end

    % methods(Static)
    %     function obj = loadobj(s)
    %         if isstruct(s)
    %             newObj = RefData;
    %             newObj.Target = s.Target;
    %             newObj.Data = s.Data;
    %             newObj.M = s.M;
    %             newObj.N = s.N;
    %             newObj.indivTrajStartIndices = s.indivTrajStartIndices;
    %             obj = newObj;
    %         else
    %             obj = s;
    %         end
    %     end
    % end
end