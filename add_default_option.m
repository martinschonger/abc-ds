% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function options = add_default_option(options,fieldname,value)

if ~isfield(options,fieldname)
    options.(fieldname) = value;
end

end
