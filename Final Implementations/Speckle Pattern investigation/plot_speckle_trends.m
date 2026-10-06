%% PLOT_SPECKLE_TRENDS  Measured speckle size against the swept variable.
%
% Reads the measured sizes from speckle_size_batch.csv and the optical
% configuration of each point from test_matrix.csv, then draws one curve per
% sweep: the points of a sweep share every setting except the variable on
% the horizontal axis.
%
% Three figures are produced.
%   1. Ground glass distance sweep, points 5 and 7 to 10 at a fixed f-number
%      of 11.31 and a fixed magnification of 0.274.
%   2. f-number sweeps, the four groups in which the ground glass distance
%      was held fixed. Predicted sizes are drawn as dashed lines.
%   3. Magnification sweep, points 14 to 28 at a fixed f-number of 22.63.
%
% To add a sweep, append an entry to the SERIES array below. Points in an
% entry are sorted by the chosen x variable before plotting.
%
% Output: three figures, saved as PNG next to the data.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc; close all;

folder = fileparts(mfilename('fullpath'));
B = readtable(fullfile(folder, 'speckle_size_batch.csv'));
T = readtable(fullfile(folder, 'test_matrix.csv'));

opt.savePng   = false;
opt.showTheory = false;      % overlay the Rayleigh limit prediction
opt.showInterp = true;       % interpolated direct speckles on figure 2
opt.showTitle  = false;      % the report carries a caption instead
opt.figSize    = [5.5 2.4];  % [in]; ~61 pt per inch on the page, so this is close to 0.75 of the text width
opt.fontSize   = 9;          % [pt]; the page enlarges the export by ~1.17
opt.lineWidth  = 1.5;
opt.markerSize = 6;
opt.yPad       = 0.12;       % vertical margin, as a fraction of the data span
opt.legendSpace  = 0.7;     % [in] extra height for a legend below the axes
opt.legendColumns = 2;      % entries per row in that legend

% Definition of the sweeps. XVAR selects the column of test_matrix.csv used
% on the horizontal axis, and the remaining settings are constant within an
% entry by construction of the test matrix.
series(1).figure = 1;
series(1).points = [5 7 8 9 10];
series(1).xvar   = 'd_GG_mm';
series(1).xlabel = 'd_{GG}  [mm]';
series(1).name   = 'f-number = 11.31,  M_{ag} = 0.296';
series(1).title  = 'Ground glass distance sweep';

series(2).figure = 2;
series(2).points = [1 2 3];
series(2).xvar   = 'f_number';
series(2).xlabel = 'f-number  [--]';
series(2).name   = 'd_{GG} = 740 mm,  M_{ag} = 0.964';
series(2).title  = 'f-number sweeps at fixed ground glass distance';

series(3).figure = 2;
series(3).points = [4 5 6];
series(3).xvar   = 'f_number';
series(3).xlabel = 'f-number  [--]';
series(3).name   = 'd_{GG} = 740 mm,  M_{ag} = 0.296';
series(3).title  = 'f-number sweeps at fixed ground glass distance';

series(4).figure = 2;
series(4).points = [11 12 13];
series(4).xvar   = 'f_number';
series(4).xlabel = 'f-number  [--]';
series(4).name   = 'd_{GG} = 300 mm,  M_{ag} = 0.877';
series(4).title  = 'f-number sweeps at fixed ground glass distance';

series(5).figure = 2;
series(5).points = [29 30 31];
series(5).xvar   = 'f_number';
series(5).xlabel = 'f-number  [--]';
series(5).name   = 'd_{GG} = 300 mm,  M_{ag} = 0.457 (projected)';
series(5).title  = 'f-number sweeps at fixed ground glass distance';

series(6).figure = 3;
series(6).points = 14:28;
series(6).xvar   = 'M_ag';
series(6).xlabel = 'M_{ag}  [--]';
series(6).name   = 'f-number = 22.63,  d_{GG} = 300 mm';
series(6).title  = 'Magnification sweep';

col = lines(7);

for i = 1:numel(series)
    s = series(i);
    [x, y, yth] = local_gather(B, T, s.points, s.xvar);

    f = figure(s.figure);
    set(f, 'Color', 'w', 'Name', s.title);
    hold on; grid on; box on;

    c = col(mod(i-1, size(col,1)) + 1, :);
    plot(x, y, '-o', 'Color', c, 'MarkerFaceColor', c, ...
         'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', s.name);

    if opt.showTheory && numel(unique(yth)) > 1
        plot(x, yth, '--', 'Color', c, 'LineWidth', 1.0, ...
             'HandleVisibility', 'off');
    end

    xlabel(s.xlabel);
    ylabel('\Delta s_{exp}  [px]');
    title(s.title);
    ylim([0 inf]);
    xl = xlim; xlim([0 xl(2)]);
    if numel(series([series.figure] == s.figure)) > 1
        legend('Location', 'best');
    end
end

% Direct speckles at the configuration of the projected points 29 to 31
% (d_GG = 300 mm, M_ag = 0.457), interpolated from measured direct points:
%   f/4     from points 1, 4 and 11  (magnification correction from point 11)
%   f/11.31 from points 2, 5 and 8   (plane through the three points)
%   f/22.63 from points 13 and 20    (linear in magnification)
if opt.showInterp
    figure(2);
    ci = col(6,:);
    plot([4.00 11.31 22.63], [2.34 2.47 3.67], '-.d', 'Color', ci, ...
         'MarkerFaceColor', ci, 'LineWidth', 1.4, 'MarkerSize', 7, ...
         'DisplayName', 'd_{GG} = 300 mm,  M_{ag} = 0.457 (direct, interpolated)');
    legend('Location', 'best');
end

if opt.showTheory
    figure(2);
    plot(nan, nan, 'k--', 'LineWidth', 1.0, ...
         'DisplayName', '\Delta s predicted, Equation 4.x');
    legend('Location', 'best');
end

% Styling for print. The figures go into the report at \textwidth, about
% 5.4 in wide, so a 7.5 in export is reduced by roughly 0.72 and a 15 pt
% label lands near 11 pt on the page. Lines and markers are thickened by the
% same reasoning, and the axes title is dropped because the report carries a
% caption instead.
for k = unique([series.figure])
    f = figure(k);
    % Figures with a legend get extra height, because the legend is placed
    % under the axes rather than inside them, where it would cover the data.
    sz = opt.figSize;
    if numel(series([series.figure] == k)) > 1
        sz(2) = sz(2) + opt.legendSpace;
    end
    set(f, 'Units', 'inches', 'Position', [1 1 sz]);
    ax = findall(f, 'Type', 'axes');
    if ~opt.showTitle, title(ax, ''); end
    set(ax, 'FontSize', opt.fontSize, 'LineWidth', 0.9, 'GridAlpha', 0.18);
    set(findall(f, 'Type', 'line'),  'LineWidth', opt.lineWidth, ...
                                     'MarkerSize', opt.markerSize);
    set(findall(f, 'Type', 'text'),  'FontSize', opt.fontSize);
    lg = findall(f, 'Type', 'legend');
    for h = reshape(lg, 1, [])
        h.FontSize  = opt.fontSize - 2;
        h.Location  = 'southoutside';
        h.NumColumns = opt.legendColumns;
        h.Box       = 'off';
    end
    set([get(ax, 'XLabel'); get(ax, 'YLabel')], 'FontSize', opt.fontSize);

    % Vertical range. Anchoring the axis at zero leaves most of the panel
    % empty when the data covers a narrow band, so the single sweeps are
    % framed around their own data with equal margins above and below. The
    % f-number figure keeps the zero baseline, since its legend sits in the
    % lower corner, and only gains headroom above the highest point.
    yd = [];
    for h = reshape(findall(ax, 'Type', 'line'), 1, [])
        v = h.YData(:);  yd = [yd; v(isfinite(v))]; %#ok<AGROW>
    end
    if ~isempty(yd)
        lo = min(yd);  hi = max(yd);  span = max(hi - lo, eps);
        if k == 2
            ylim(ax, [0, hi + opt.yPad*hi]);
        else
            ylim(ax, [lo - opt.yPad*span, hi + opt.yPad*span]);
        end
    end
end

if opt.savePng
    for k = unique([series.figure])
        f = figure(k);
        name = fullfile(folder, sprintf('speckle_trend_fig%d.png', k));
        exportgraphics(f, name, 'Resolution', 300);
        fprintf('Saved %s\n', name);
    end
end

%% ---------------------------------------------------------------- helper

function [x, y, yth] = local_gather(B, T, points, xvar)
% Measured diameter, prediction and abscissa for a list of points, sorted by
% the abscissa so the connecting line runs in order.
    n = numel(points);
    x = nan(n,1); y = nan(n,1); yth = nan(n,1);
    for i = 1:n
        p  = points(i);
        jb = find(B.point == p, 1);
        jt = find(T.point == p, 1);
        if isempty(jb) || isempty(jt)
            warning('Point %d is missing from one of the input files.', p);
            continue
        end
        x(i)   = T.(xvar)(jt);
        y(i)   = B.twoa_1e_px(jb);
        yth(i) = T.delta_s_th_px(jt);
    end
    [x, ord] = sort(x);
    y = y(ord);  yth = yth(ord);
end
