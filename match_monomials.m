% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [ind_from, ind_to] = match_monomials(monolist_from, monolist_to)
% returns indices in monolist_to which correspond to in 

ind_from = [];
ind_to = [];

for i = 1:length(monolist_from)
    tmp = find(strcmp(monolist_to, monolist_from{i}));
    if ~isempty(tmp)
        ind_from(end+1) = i;
        ind_to(end+1) = tmp;
        assert(length(tmp) == 1, 'Found multiple correspondences.');
    else
        warning('Did not find correspondence for: %s', monolist_from{i});
    end
end

end