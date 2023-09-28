% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [options, result, f, V, B] = read_exp(timestamp_str)

options_json_filename = "output/" + timestamp_str + "_options.json";
options_json_raw = fileread(options_json_filename);
options = jsondecode(options_json_raw);

result_json_filename = "output/" + timestamp_str + "_result.json";
result_json_raw = fileread(result_json_filename);
result = jsondecode(result_json_raw);

f = vectorize_mvfun_str(result.f_fh_str_arr);
V = vectorize_mvfun_str(result.V_fh_str_arr);
B = vectorize_mvfun_str(result.B_fh_str_arr);

end