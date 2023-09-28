% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function figV = plot_V(plot_dim1, plot_dim2, X, Y, Vfun_eval, axis_limits, streamlines_plt, filename, image_resolution)

figV = figure;

surf(X, Y, Vfun_eval, 'edgecolor', 'none', 'displayname', '$V(\xi)$');

cmap = crameri('lajolla');
cmap = cmap(1:160, :);
colormap(cmap);

hold on;

[~, isolines_V] = contour3(X, Y, Vfun_eval, 'edgecolor', 'black');
set(isolines_V, 'handlevisibility', 'off');

view([-227 67]);

leg = legend();
set(leg, 'Interpreter', 'latex');

xlabel(strcat('\xi_', int2str(plot_dim1)));
ylabel(strcat('\xi_', int2str(plot_dim2)));

% axis limits (depending on ref data)
xlim([axis_limits(1), axis_limits(2)]);
ylim([axis_limits(3), axis_limits(4)]);

% legend config
leg = findobj(gcf, 'Type', 'Legend');
leg.ItemTokenSize(1) = 15;

axV = gca;
streamlines_cpyV_plt = copyobj(streamlines_plt, axV);

for i=1:length(streamlines_cpyV_plt)
    zi = interp2(X, Y, Vfun_eval, streamlines_cpyV_plt(i).XData, streamlines_cpyV_plt(i).YData);
    streamlines_cpyV_plt(i).ZData = zi;
end


exportgraphics(figV, filename, 'contenttype', 'image', 'resolution', image_resolution); % TODO: change to 'vector'

end