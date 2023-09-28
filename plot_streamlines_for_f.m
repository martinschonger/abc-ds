% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function [plt_handle, plt_data] = plot_streamlines_for_f(f, limits, res, density)

arguments
    f;
    limits;
    res = 200;
    density = 5;
end

[grid_x, grid_y] = meshgrid(linspace(limits(1), limits(2), res), linspace(limits(3), limits(4), res));

plt_data = feval(f, [grid_x(:), grid_y(:)]');

plt_handle = streamslice(grid_x, grid_y, reshape(plt_data(1, :), res, res), reshape(plt_data(2, :), res, res), density, 'method', 'cubic');
set(plt_handle, 'linewidth', 0.5);
set(plt_handle, 'color', 'black');

end