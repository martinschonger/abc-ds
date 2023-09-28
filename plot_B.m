% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function figB = plot_B(plot_dim1, plot_dim2, X, Y, Bfun_eval, axis_limits, streamlines_plt, filename, image_resolution)

figB = figure;
        
surf(X, Y, Bfun_eval, 'edgecolor', 'none', 'displayname', '$B(\xi)$');

cmap = crameri('-devon');
cmap = cmap(1:end-8, :);
colormap(cmap);

hold on;

[~, isolines_B] = contour3(X, Y, Bfun_eval, 'edgecolor', 'white');
set(isolines_B, 'handlevisibility', 'off');

contour3(X, Y, Bfun_eval, [0, 0], 'color', '#B40101', 'linewidth', 1, 'displayname', '$B^{-1}(0)$');
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

axB = gca;
streamlines_cpyB_plt = copyobj(streamlines_plt, axB);

for i=1:length(streamlines_cpyB_plt)
    zi = interp2(X, Y, Bfun_eval, streamlines_cpyB_plt(i).XData, streamlines_cpyB_plt(i).YData);
    streamlines_cpyB_plt(i).ZData = zi;
end


exportgraphics(figB, filename, 'contenttype', 'image', 'resolution', image_resolution); % TODO: change to 'vector'

end