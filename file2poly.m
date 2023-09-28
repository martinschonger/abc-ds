% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [coefs, monomials] = file2poly(filename, xi)

tmp_json_in_filename = filename;
tmp_json_in_raw = fileread(tmp_json_in_filename);
tmp_json_in = jsondecode(tmp_json_in_raw);

num_monos = length(tmp_json_in.monomials);
monomial_list = sdpvar(num_monos, 1, 'full');
% xi needs to exist as sdpvar
for i = 1:num_monos
    monomial_list(i) = eval(tmp_json_in.monomials{i});
end
sdisplay(monomial_list);

coefs = tmp_json_in.coefs;
monomials = monomial_list;

end