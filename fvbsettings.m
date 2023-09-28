% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function options = fvbsettings(varargin)

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'seed', 13);
addParameter(p, 'enable_barrier', true);
addParameter(p, 'epsilon', 1e-4);  % required for numerical stability of the optimization
addParameter(p, 'init_Bc_var', true);  % currently unused
addParameter(p, 'constraint_version', 1);
addParameter(p, 'enable_regularization', true);
addParameter(p, 'deg_f', 2);
addParameter(p, 'deg_V', 2);
addParameter(p, 'deg_B', 2);
addParameter(p, 'deg_B_slack', 2);
addParameter(p, 'enable_extra_constraint', true);  % explicitly force B to be <= 0 at the reference trajectories
addParameter(p, 'regularization_factor', 0.01);
addParameter(p, 'dataset', 'lasa');  % can also be set to 'robot'
addParameter(p, 'dataset_opts', struct('idx', 14));  % see `main2.m` for other possible options
addParameter(p, 'sdpoptions', struct);  % passed to YALMIP
addParameter(p, 'sdpoptions_penbmi', struct);  % passed to PENBMI
parse(p, varargin{:});

options = p.Results;
if ~isempty(fieldnames(p.Unmatched))
    options.unmatched = p.Unmatched;
end

% check with schema
str = fileread("fvbsettings_schema.txt");
schema_fieldnames = regexp(str, '\r\n|\r|\n', 'split')';

actual_fieldnames = fieldnames(options);
diff_to_schema1 = setdiff(actual_fieldnames, schema_fieldnames);
diff_to_schema2 = setdiff(schema_fieldnames, actual_fieldnames);

if ~(isempty(diff_to_schema1) && isempty(diff_to_schema2))
    warning('fvbsettings_schema violated');
end

end