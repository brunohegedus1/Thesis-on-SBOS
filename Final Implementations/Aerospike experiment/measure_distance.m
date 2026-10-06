function d = measure_distance(ax, x, y)
%MEASURE_DISTANCE  Distance between two points clicked on a plot.
%
%   MEASURE_DISTANCE waits for two clicks on the current axes, draws the
%   line between them, labels it with the distance and prints the result.
%   It then waits for the next pair of clicks, so several distances can be
%   measured in a row. Press Enter or Escape to stop.
%
%   MEASURE_DISTANCE(AX) measures on the axes AX instead of the current
%   one. AX may also be a figure, in which case its current axes is used.
%
%   MEASURE_DISTANCE(AX, X, Y) skips the clicking and measures between the
%   two points given by the two-element vectors X and Y. Useful for
%   repeating a measurement from a script.
%
%   D = MEASURE_DISTANCE(...) returns the last distance measured, in the
%   units of the axes (mm for the aerospike plots).
%
%   Example, after running plot_aerospike_data or aerospike_density_gradient:
%       measure_distance          % click two points on the figure
%
%   The printout gives the two points, the x and y separation, the distance
%   and the angle of the line, measured counter-clockwise from the x-axis.

    if nargin < 1 || isempty(ax)
        ax = gca;
    end
    if isgraphics(ax, 'figure')
        ax = get(ax, 'CurrentAxes');
    end
    if ~isgraphics(ax, 'axes')
        error('measure_distance:badAxes', ...
            'First argument must be an axes or a figure holding one.');
    end

    d = NaN;

    % Given points: measure once and return.
    if nargin >= 3
        if numel(x) ~= 2 || numel(y) ~= 2
            error('measure_distance:badPoints', ...
                'X and Y must each hold two values.');
        end
        d = drawMeasurement(ax, x, y);
        return
    end

    % Interactive: one measurement per pair of clicks.
    fig = ancestor(ax, 'figure');
    figure(fig);
    fprintf('Click two points to measure. Press Enter or Escape to stop.\n');
    while isgraphics(ax)
        [x, y] = ginput(2);
        if numel(x) < 2
            break
        end
        d = drawMeasurement(ax, x, y);
    end
    fprintf('Done measuring.\n');
end

function d = drawMeasurement(ax, x, y)
% Draws the line, labels it and prints the numbers.
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    d  = hypot(dx, dy);
    ang = atan2d(dy, dx);

    held = ishold(ax);
    hold(ax, 'on');
    plot(ax, x, y, 'o-', 'Color', [0 0 0], 'MarkerSize', 5, ...
        'MarkerFaceColor', [0 0 0], 'LineWidth', 1.2);
    text(ax, mean(x), mean(y), sprintf('  %.3f', d), ...
        'Color', [0 0 0], 'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 ], 'Margin', 1, ...
        'VerticalAlignment', 'bottom');
    if ~held
        hold(ax, 'off');
    end

    fprintf(['  (%.3f, %.3f) to (%.3f, %.3f):  dx %+.3f  dy %+.3f  ' ...
             'distance %.3f  angle %+.2f deg\n'], ...
        x(1), y(1), x(2), y(2), dx, dy, d, ang);
end
