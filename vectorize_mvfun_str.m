% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [func_handle] = vectorize_mvfun_str(mvfun_str_cellarr, num_dims)

arguments
    mvfun_str_cellarr;
    num_dims = 2;
end

if ~isa(mvfun_str_cellarr, "cell")
    mvfun_str_cellarr = {mvfun_str_cellarr};
end

mvfun_str_cellarr = replace(mvfun_str_cellarr, '*', '.*');
mvfun_str_cellarr = replace(mvfun_str_cellarr, '^', '.^');
for i = 1:num_dims
    mvfun_str_cellarr = replace(mvfun_str_cellarr, strcat("xi", int2str(i)), strcat("xi(", int2str(i), ",:)"));
end

num_out_dims = length(mvfun_str_cellarr);
mvfun_str = "";
for i = 1:num_out_dims
    mvfun_str = append(mvfun_str, mvfun_str_cellarr{i});
    if i < num_out_dims
        mvfun_str = append(mvfun_str, "; ");
    end
end

func_handle = eval("@(xi) [" + mvfun_str + "]");

end