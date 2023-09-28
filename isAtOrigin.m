% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [value, isterminal, direction] = isAtOrigin(T, Y)
    epsilon = 1e-4;
    value      = all(Y > -epsilon) && all(Y < epsilon);
    isterminal = 1;  % Stop the integration
    direction  = 0;
end