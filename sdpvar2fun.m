% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [func_handle] = sdpvar2fun(sdpexpr, xi)

func_handle = sdisplay(sdpexpr);
func_handle = func_handle{1};
func_handle = replace(func_handle, '*', '.*');
func_handle = replace(func_handle, '^', '.^');

for i = 1:length(xi)
    func_handle = replace(func_handle, strcat("xi(", int2str(i), ")"), strcat("xi(", int2str(i), ",:)"));
end

func_handle = eval(['@(xi)', func_handle]);

end