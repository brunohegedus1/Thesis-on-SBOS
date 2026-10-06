%% PLOT_SPECKLE_METHODS  Measured speckle size per point, both estimators.
%
% Plots the average speckle size of every point of the speckle investigation
% campaign as returned by the one-dimensional autocorrelation method of
% Hamarova et al. (2016) and by the DaVis estimator. Values are those of
% Tables 6.7 (points 1 to 28, directed speckles) and 6.8 (points 29 to 31,
% projected speckles) of the thesis. DaVis returned no value for point 10,
% which is left as a gap.
%
% Dashed vertical lines separate the three blocks of the test matrix:
% aperture and ground glass sweeps (1 to 13), magnification sweep (14 to 28)
% and projected speckles (29 to 31).
%
% Output: one figure, saved as PNG next to this script.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc; close all;

folder = fileparts(mfilename('fullpath'));

opt.savePng    = true;
opt.bars       = 'none';     % 'dev'  : bar from each value to Delta s_th
                            % 'utot' : u_tot on the Hamarova curve
                            % 'none' : no bars
opt.xOffset    = 0.15;      % horizontal shift so the two sets of bars separate
opt.outFile    = fullfile(folder, 'speckle_methods_per_point.png');

% Styling for print, as in plot_speckle_trends.m. The report scales the
% export to \textwidth, which enlarges it by about 1.17, so a 9 pt label
% lands near the 11 pt body text and the 8 pt legend near the caption.
opt.figSize    = [6.0 3.1];  % [in]; close to 0.85 of the text width, with room for the legend below the axes
opt.fontSize   = 9;          % [pt]
opt.lineWidth  = 1.3;
opt.markerSize = 5;

pt = (1:31)';

% Delta s_exp [px], Tables 6.7 and 6.8
hamarova = [1.74 2.73 3.83 1.62 2.71 3.81 2.29 2.45 2.91 4.59 ...
            2.42 2.98 4.15 3.05 3.14 3.26 2.66 2.80 2.39 3.60 ...
            2.16 2.38 2.53 2.88 3.06 2.72 2.76 3.09 1.59 2.21 3.61]';
davis    = [1.20 1.48 2.00 1.19 1.56 2.04 1.35 1.40 1.59  NaN ...
            1.45 1.67 2.41 1.81 1.76 1.86 1.60 1.65 1.47 1.68 ...
            1.42 1.48 1.59 1.74 1.91 1.65 1.71 1.87 1.22 1.44 2.09]';

% Combined standard uncertainty u_tot [%] of the Hamarova value
uTot     = [6.35 1.50 11.45 2.62 7.11 4.87 6.60 2.87 2.79 3.12 ...
            1.50 1.50 1.55 1.52 1.50 1.50 6.43 4.78 2.80 2.39 ...
            4.13 3.39 4.66 3.33 3.42 3.05 1.82 5.47 1.54 1.80 2.08]';

% Deviation from the Rayleigh prediction, Delta s_exp - Delta s_th [px]
devH     = [ 0.82  0.14 -1.36  1.02  1.00  0.39  0.58  0.74  1.20  2.88 ...
             1.54  0.49 -0.83 -2.33 -1.99 -1.58 -1.81 -1.42 -1.64 -0.16 ...
            -1.42 -1.07 -0.76 -0.25  0.06 -0.17 -0.07  0.31  0.92  0.31 -0.20]';
devD     = [ 0.28 -1.11 -3.19  0.59 -0.15 -1.39 -0.36 -0.31 -0.12   NaN ...
             0.57 -0.82 -2.56 -3.57 -3.38 -2.98 -2.88 -2.57 -2.56 -2.08 ...
            -2.16 -1.97 -1.70 -1.39 -1.09 -1.23 -1.12 -0.91  0.55 -0.46 -1.72]';

col = lines(2);

f = figure('Color', 'w', 'Name', 'Speckle size per point');
hold on; grid on; box on;

xH = pt - opt.xOffset;
xD = pt + opt.xOffset;
mk = {'LineWidth', 1.4, 'MarkerSize', 6, 'CapSize', 4};

switch opt.bars
    case 'dev'
        % One-sided bars that end on the prediction: a positive deviation
        % hangs below the measured value, a negative one rises above it.
        errorbar(xH, hamarova, max(devH,0), max(-devH,0), '-o', ...
                 'Color', col(1,:), 'MarkerFaceColor', col(1,:), mk{:}, ...
                 'DisplayName', 'Hamarová et al.');
        errorbar(xD, davis, max(devD,0), max(-devD,0), '-s', ...
                 'Color', col(2,:), 'MarkerFaceColor', col(2,:), mk{:}, ...
                 'DisplayName', 'DaVis');
        yTop = max([hamarova + max(-devH,0); davis + max(-devD,0)]);
    case 'utot'
        errorbar(xH, hamarova, hamarova .* uTot/100, '-o', ...
                 'Color', col(1,:), 'MarkerFaceColor', col(1,:), mk{:}, ...
                 'DisplayName', 'Hamarová et al.');
        plot(xD, davis, '-s', 'Color', col(2,:), 'MarkerFaceColor', col(2,:), ...
             'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', 'DaVis');
        yTop = max(hamarova .* (1 + uTot/100));
    otherwise
        plot(xH, hamarova, '-o', 'Color', col(1,:), 'MarkerFaceColor', col(1,:), ...
             'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', 'Hamarová et al.');
        plot(xD, davis, '-s', 'Color', col(2,:), 'MarkerFaceColor', col(2,:), ...
             'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', 'DaVis');
        yTop = max(hamarova);
end

xlim([0.5 31.5]);
ylim([0 ceil(yTop*1.15)]);
yl = ylim;
for xb = [13.5 28.5]
    plot([xb xb], yl, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
text(7,  yl(2)*0.93, 'f-number and d_{GG} sweeps', 'HorizontalAlignment', 'center');
text(21, yl(2)*0.93, 'M_{ag} sweep',             'HorizontalAlignment', 'center');
text(31.3, yl(2)*0.93, 'projected',              'HorizontalAlignment', 'right');

set(gca, 'XTick', 1:31, 'TickLabelInterpreter', 'tex');
xlabel('Point  [--]');
ylabel('\Delta s_{exp}  [px]');
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

set(f, 'Units', 'inches', 'Position', [1 1 opt.figSize]);
ax = findall(f, 'Type', 'axes');
set(ax, 'FontSize', opt.fontSize, 'LineWidth', 0.9, 'GridAlpha', 0.18);
set(findall(f, 'Type', 'line'), 'LineWidth', opt.lineWidth, ...
                                'MarkerSize', opt.markerSize);
set(findall(f, 'Type', 'text'), 'FontSize', opt.fontSize);
set([get(ax, 'XLabel'); get(ax, 'YLabel')], 'FontSize', opt.fontSize);
for h = reshape(findall(f, 'Type', 'legend'), 1, [])
    h.FontSize = opt.fontSize - 1;
end

if opt.savePng
    exportgraphics(f, opt.outFile, 'Resolution', 300);
    fprintf('Saved %s\n', opt.outFile);
end
