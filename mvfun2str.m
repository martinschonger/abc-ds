% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [str_repr, str_arr] = mvfun2str(mvfun, precision)
% multivariate symbolic -> string

arguments
    mvfun;
    precision = 16;
end

str_repr = append(inputname(1), "(xi) = ");
indent_str = blanks(strlength(str_repr));
xi = sym('xi', [2, 1]);
str_arr = string(vpa(sym(mvfun(xi)), precision));
for m = 1:length(str_arr)
    if m == 1
        str_repr = append(str_repr, "[", str_arr(m), "]\n");
    else
        str_repr = append(str_repr, indent_str, "[", str_arr(m), "]\n");
    end
end

end