% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [f_fh, V_fh, dVdx_fh, B_fh, dBdx_fh, debug_output, fc, Vc, Bc] = fvb(rd, restrict_to_convex, xi, initial_set, unsafe_set, options)
% Run the main optimization problem to compute a polynomial dynamical system f, a polynomial
%   Lyapunov function V, and a polynomial Barrier certificate B.
%
% :param rd: RefData object encapsulating the reference trajectories.
% :param restrict_to_convex: Whether to solve the full optimization problem or a convex variant.
%   0: Solve easier convex problem with fixed Lyapunov function.
%   1: Solve full problem but initialize with solution of easier problem.
% :param xi: sdpvar corresponding to the state space.
% :param initial_set: Cell array of sdpvar expressions of polynomials (which are taken as >= 0) that
%   define a semi-algebraic set representation of the initial set.
% :param unsafe_set: Cell array of sdpvar expressions of polynomials (which are taken as >= 0) that
%   define a semi-algebraic set representation of the unsafe set.
% :param options: Options struct to configure this function's behavior, see :func:`fvbsettings` for
%   a subset of supported options.
% :returns: [f_fh, V_fh, dVdx_fh, B_fh, dBdx_fh, debug_output, fc, Vc, Bc], where the individual
%   parts elements mean the following:
%   f_fh: Function handle to the dynamical system f.
%   V_fh: Function handle to the Lyapunov function V.
%   dVdx_fh: Function handle to the Jacobian of V.
%   B_fh: Function handle to the Barrier certificate B.
%   dBdx_fh: Function handle to the Jacobian of B.
%   debug_output: Struct containing function-internal data as well as result data, such as solver output and MSE.
%   fc: Polynomial coefficients of f.
%   Vc: Polynomial coefficients of V.
%   Bc: Polynomial coefficients of B.


%
% assumption: attractor is at the origin

epsilon = options.epsilon;

debug_output = struct;


% Reference trajectories
M = rd.M; % # states
T = rd.T; % time/sample index (all demonstrations concatenated)


%% Define variables
deg_f = options.deg_f;
f = [];
fc_var = [];
for m = 1:M
    [f_tmp, fc_var_tmp, f_monomials] = polynomial(xi, deg_f, 1);
    f = [f; f_tmp];
    fc_var = [fc_var; fc_var_tmp];
end
debug_output.f_monomials = f_monomials;

% Lyapunov function
if restrict_to_convex == 0
    V = xi' * xi;
    [Vc_var, V_monomials] = coefficients(V, xi);
elseif restrict_to_convex == 1
    deg_V = options.deg_V;
    [V, Vc_var, V_monomials] = polynomial(xi, deg_V, 1);
end
debug_output.V_monomials = V_monomials;
dVdx = jacobian(V, xi)';

% Barrier certificate
deg_B = options.deg_B;
[B, Bc_var, B_monomials] = polynomial(xi, deg_B);
debug_output.B_monomials = B_monomials;
dBdx = jacobian(B, xi)';


% init sdpvars from options struct
if isfield(options, 'unmatched')
    if isfield(options.unmatched, 'fc_init')
        assert(isfield(options.unmatched, 'fc_init_monomials'), 'options.unmatched.fc_init_monomials is required when options.unmatched.fc_init is set');
        assert(isfield(options.unmatched, 'fc_init_sdpvar'), 'options.unmatched.fc_init_sdpvar is required when options.unmatched.fc_init is set');

        fc_init_tmp = zeros(length(f_monomials), M);
        % Assign the computed coefficient values to the matching coefficients of the higher-deg polynomial
        fc_init_reshaped = reshape(options.unmatched.fc_init, length(options.unmatched.fc_init)/M, M);
        fc_init_monomials_replaced = replace(options.unmatched.fc_init_monomials, options.unmatched.fc_init_sdpvar, xi);
        [idx_from, idx_to] = match_monomials(sdisplay(fc_init_monomials_replaced), sdisplay(f_monomials));
        for m = 1:M
            fc_init_tmp(idx_to, m) = fc_init_reshaped(idx_from, m);
        end
        assign(fc_var, fc_init_tmp(:));
    end
    if isfield(options.unmatched, 'Vc_init')
        assert(isfield(options.unmatched, 'Vc_init_monomials'), 'options.unmatched.Vc_init_monomials is required when options.unmatched.Vc_init is set');
        assert(isfield(options.unmatched, 'Vc_init_sdpvar'), 'options.unmatched.Vc_init_sdpvar is required when options.unmatched.Vc_init is set');

        Vc_init_tmp = zeros(size(Vc_var));
        % Assign the computed coefficient values to the matching coefficients of the higher-deg polynomial
        Vc_init_monomials_replaced = replace(options.unmatched.Vc_init_monomials, options.unmatched.Vc_init_sdpvar, xi);
        [idx_from, idx_to] = match_monomials(sdisplay(Vc_init_monomials_replaced), sdisplay(V_monomials));
        Vc_init_tmp(idx_to) = options.unmatched.Vc_init(idx_from);
        assign(Vc_var, Vc_init_tmp);
    end
    if isfield(options.unmatched, 'Bc_init')
        assert(isfield(options.unmatched, 'Bc_init_monomials'), 'options.unmatched.Bc_init_monomials is required when options.unmatched.Bc_init is set');
        assert(isfield(options.unmatched, 'Bc_init_sdpvar'), 'options.unmatched.Bc_init_sdpvar is required when options.unmatched.Bc_init is set');

        Bc_init_tmp = zeros(size(Bc_var));
        % Assign the computed coefficient values to the matching coefficients of the higher-deg polynomial
        Bc_init_monomials_replaced = replace(options.unmatched.Bc_init_monomials, options.unmatched.Bc_init_sdpvar, xi);
        [idx_from, idx_to] = match_monomials(sdisplay(Bc_init_monomials_replaced), sdisplay(B_monomials));
        Bc_init_tmp(idx_to) = options.unmatched.Bc_init(idx_from);
        assign(Bc_var, Bc_init_tmp);
    end
end


%% Objective
xi_dot = sdpvar(M, T, 'full');
for t = 1:T
    xi_dot(:, t) = replace(f, xi, rd.Data(1:M, t));
end
xi_dot_error = xi_dot - rd.Data(M+1:end, :);

mse = sum(sum(xi_dot_error.^2)) / (2 * T);
Objective = mse;


%% Constraints
Constraints = [];
epxi = epsilon * sum(xi.^2, 1);

% Stability
if restrict_to_convex == 1
    Constraints = [Constraints, (sos(V-epxi)):'V > 0'];
end
Vdot = sum(dVdx.*f, 1);
% Lie derivative
Constraints = [Constraints, (sos(-Vdot-epxi)):'-Vdot < 0'];

% Barrier
if options.enable_barrier
    deg_B_slack = options.deg_B_slack;

    % initial set
    tau_arr = {};
    tauc_arr = {};
    tau_monomials_arr = {};
    sos_safe = -B;
    for p = 1:length(initial_set)
        [tau_arr{p}, tauc_arr{p}, tau_monomials_arr{p}] = polynomial(xi, deg_B_slack);
        Constraints = [Constraints, sos(tau_arr{p})];
        sos_safe = sos_safe - tau_arr{p} * initial_set{p};
    end
    Constraints = [Constraints, sos(sos_safe)];
    debug_output.tau_arr = tau_arr;
    debug_output.tauc_arr = tauc_arr;
    debug_output.tau_monomials_arr = tau_monomials_arr;

    % unsafe set
    sig_arr = {};
    sigc_arr = {};
    sig_monomials_arr = {};
    sos_unsafe = B - epsilon;
    for m = 1:length(unsafe_set)
        [sig_arr{m}, sigc_arr{m}, sig_monomials_arr{m}] = polynomial(xi, deg_B_slack);
        Constraints = [Constraints, sos(sig_arr{m})];
        sos_unsafe = sos_unsafe - sig_arr{m} * unsafe_set{m};
    end
    Constraints = [Constraints, sos(sos_unsafe)];
    debug_output.sig_arr = sig_arr;
    debug_output.sigc_arr = sigc_arr;
    debug_output.sig_monomials_arr = sig_monomials_arr;

    % Lie derivative
    Bdot = sum(dBdx.*f, 1);
    switch options.constraint_version
        case 1
            Constraints = [Constraints, sos(-Bdot-epxi)];
        case 2
            [slackvar, slackvarc, slackvar_monomials] = polynomial(xi, deg_B_slack);
            Constraints = [Constraints, sos(slackvar)];
            Constraints = [Constraints, sos(-Bdot-epxi-slackvar*(epsilon - B))];
            debug_output.slackvar = slackvar;
            debug_output.slackvarc = slackvarc;
            debug_output.slackvar_monomials = slackvar_monomials;
        case 3
            [slackvar, slackvarc, slackvar_monomials] = polynomial(xi, deg_B_slack);
            Constraints = [Constraints, sos(slackvar)];
            Constraints = [Constraints, sos(-Bdot-epxi-slackvar*(epsilon - B.^2))]; % restrict to B^{-1}(0)
            debug_output.slackvar = slackvar;
            debug_output.slackvarc = slackvarc;
            debug_output.slackvar_monomials = slackvar_monomials;
    end

    if options.enable_extra_constraint
        B_tmp3 = 0;
        for t = 1:T
            B_tmp3 = B_tmp3 + max([0, replace(B, xi, rd.Data(1:M, t))]);
        end
        Constraints = [Constraints, B_tmp3 <= 0]; % TODO: maybe use -epsilon on right side
    end
end

%% Optimize
params = [];
params = [params; fc_var(:)];
if restrict_to_convex == 1
    params = [params; Vc_var(:)];
end
if options.enable_barrier
    params = [params; Bc_var];
    for p = 1:length(initial_set)
        params = [params; tauc_arr{p}];
    end
    for m = 1:length(unsafe_set)
        params = [params; sigc_arr{m}];
    end
    if options.constraint_version == 2 || options.constraint_version == 3
        params = [params; slackvarc];
    end
end


if options.enable_regularization
    Objective = Objective + options.regularization_factor * norm(fc_var(:), 'fro')^2;
end


tmp_sdpoptions = options.sdpoptions;
tmp_sdpoptions = namedargs2cell(tmp_sdpoptions);
tmp_sdpoptions_penbmi = options.sdpoptions_penbmi;
tmp_sdpoptions_penbmi = namedargs2cell(tmp_sdpoptions_penbmi);
for idx = 1:length(tmp_sdpoptions_penbmi)
    if mod(idx, 2) == 1
        tmp_sdpoptions_penbmi(idx) = append('penbmi.', tmp_sdpoptions_penbmi(idx));
    end
end

final_sdpoptions = [tmp_sdpoptions, tmp_sdpoptions_penbmi];


% Solver options
if restrict_to_convex == 0
    sdp_options = sdpsettings('verbose', 1);
elseif restrict_to_convex == 1
    sdp_options = sdpsettings('solver', 'penbmi', 'verbose', 1, 'debug', 1, 'showprogress', 1, 'usex0', 1, 'warmstart', 1);
    %                                                                           Default     More Info
    % sdp_options = sdpsettings(sdp_options, 'penbmi.PBM_MAX_ITER', 10000);     % 50        maximum number of iterations of the overall algorithm
    sdp_options = sdpsettings(sdp_options, 'penbmi.OUTPUT', 3);                 % 1
    % sdp_options = sdpsettings(sdp_options, 'penbmi.ALPHA', 1e-6);             % 0.01      stopping criterium for unconstrained minimization; YALMIP lower bounds it by 1e-6
    
    % sdp_options = sdpsettings(sdp_options, 'penbmi.UM_MAX_ITER', 1000);       % 100       maximum number of iterations for the unconstrained minimization
    % sdp_options = sdpsettings(sdp_options, 'penbmi.LS', 1);                   % 0
    % sdp_options = sdpsettings(sdp_options, 'penbmi.NWT_SYS_MODE', 1);         % 0
    % sdp_options = sdpsettings(sdp_options, 'penbmi.PREC_TYPE', 1);            % 0
    % sdp_options = sdpsettings(sdp_options, 'penbmi.TR_MODE', 1);              % 0
    % sdp_options = sdpsettings(sdp_options, 'penbmi.PRECISION', 1e-6);         % 1e-6      stopping criterium for the overall algorithm (likely corresponds to PBM_EPS in PENBMI doc)
    % sdp_options = sdpsettings(sdp_options, 'penbmi.P_EPS', 1e-6);             % 1e-4
    % sdp_options = sdpsettings(sdp_options, 'penbmi.P0', 0.01);                % 0.1
    % sdp_options = sdpsettings(sdp_options, 'penbmi.PEN_UP', 0.0);             % 0.5
    % if options.enable_barrier
    %     warning('Check if it makes sense to lower the initial penalty when including a Barrier function -> is the initial guess feasible?');
    % end
    % sdp_options = sdpsettings(sdp_options, 'penbmi.ALPHA_UP', 1.0);           % 1.0
    % sdp_options = sdpsettings(sdp_options, 'penbmi.PRECISION_2', 1e-2);       % 1e-6      precision of the KKT conditions

    sdp_options = sdpsettings(sdp_options, final_sdpoptions{:});
end


% Call solver
if sum(logical(is(Constraints, 'sos'))) > 0
    [sol, v, Q, res] = solvesos(Constraints, Objective, sdp_options, params);
else
    sol = optimize(Constraints, Objective, sdp_options);
    v = [];
    Q = [];
    res = [];
end
debug_output.v = v;
debug_output.Q = Q;
debug_output.res = res;

if sol.problem ~= 0
    yalmiperror(sol.problem);
end


% Optimization result
sol.info
check(Constraints)
fprintf('Total error: %2.2f\nComputation Time: %2.2f\n', value(Objective), sol.solvertime);
debug_output.sol = sol;


% Output variables
f_fh = replace(f, fc_var, value(fc_var));
f_fhm = {};
f_fh_str = '@(xi) [';
for m = 1:M
    f_fhm{m} = sdpvar2fun(f_fh(m), xi);
    f_fh_str = strcat(f_fh_str, 'f_fhm{', num2str(m), '}(xi);');
end
f_fh_str = strcat(f_fh_str, ']');
f_fh = eval(f_fh_str);

if restrict_to_convex == 1
    V_fh = replace(V, Vc_var, value(Vc_var));
    V_fh = sdpvar2fun(V_fh, xi);

    dVdx_fh = replace(dVdx, Vc_var, value(Vc_var));
    dVdx_fhm = {};
    dVdx_fh_str = '@(xi) [';
    for m = 1:M
        dVdx_fhm{m} = sdpvar2fun(dVdx_fh(m), xi);
        dVdx_fh_str = strcat(dVdx_fh_str, 'dVdx_fhm{', num2str(m), '}(xi);');
    end
    dVdx_fh_str = strcat(dVdx_fh_str, ']');
    dVdx_fh = eval(dVdx_fh_str);
else % TODO: unify with above case
    V_fh = sdpvar2fun(V, xi);

    dVdx_fhm = {};
    dVdx_fh_str = '@(xi) [';
    for m = 1:M
        dVdx_fhm{m} = sdpvar2fun(dVdx(m), xi);
        dVdx_fh_str = strcat(dVdx_fh_str, 'dVdx_fhm{', num2str(m), '}(xi);');
    end
    dVdx_fh_str = strcat(dVdx_fh_str, ']');
    dVdx_fh = eval(dVdx_fh_str);
end

B_fh = replace(B, Bc_var, value(Bc_var));
B_fh = sdpvar2fun(B_fh, xi);

dBdx_fh = replace(dBdx, Bc_var, value(Bc_var));
dBdx_fhm = {};
dBdx_fh_str = '@(xi) [';
for m = 1:M
    dBdx_fhm{m} = sdpvar2fun(dBdx_fh(m), xi);
    dBdx_fh_str = strcat(dBdx_fh_str, 'dBdx_fhm{', num2str(m), '}(xi);');
end
dBdx_fh_str = strcat(dBdx_fh_str, ']');
dBdx_fh = eval(dBdx_fh_str);

fc = value(fc_var);
Vc = value(Vc_var);
Bc = value(Bc_var);


[f_fh_str, ~] = mvfun2str(f_fh);
fprintf(f_fh_str);
fprintf(mvfun2str(V_fh));
fprintf(mvfun2str(dVdx_fh));
fprintf(mvfun2str(B_fh));
fprintf(mvfun2str(dBdx_fh));

debug_output.mse = value(mse);


% For SOS constraints: eigenvalues of Gramians -> compare with size of residual
fprintf('[DEBUG] Validation of SOS constraints:\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('| Idx| Min eigval of Q| Primal residual|  Relative fact.|\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
min_q_eigvals = [];
primal_res = check(Constraints(logical(is(Constraints, 'sos'))));
relative_factors = [];
for i = 1:length(Q)
    min_eigval = min(eig(Q{i}));
    min_q_eigvals(end+1) = min_eigval;
    rel_fact = primal_res(i) / min_eigval;
    relative_factors(end+1) = rel_fact;
    fprintf('|%4.0f|%16g|%16g|%16g|  All eigvals of Q: %s\n', i, min_eigval, primal_res(i), rel_fact, mat2str(eig(Q{i})));
end
debug_output.min_q_eigvals = min_q_eigvals;
debug_output.primal_residuals = primal_res';
debug_output.relative_factors = relative_factors;
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('[DEBUG_END]\n');

end