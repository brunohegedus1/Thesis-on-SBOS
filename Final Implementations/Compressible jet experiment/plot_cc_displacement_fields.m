% PLOT_CC_DISPLACEMENT_FIELDS
%
%   Reads the six DaVis cross-correlation exports in CC_DIR and plots each
%   displacement field in its own figure. Every figure has three panels:
%   x-displacement, y-displacement and displacement magnitude. Set
%   MAGNITUDE_ONLY to true to plot only the magnitude panel.
%
%   All six figures share one colour scale per quantity, so you can compare
%   them directly. The limits come from the 99.5th percentile over all files
%   of the valid vectors, rounded up to a round number. The few spurious
%   vectors near 0.05 mm (about 3 px) saturate instead of washing out the
%   rest of the field, which sits well below one pixel.
%
%   x and y displacement use a symmetric diverging map centred on zero.
%   The magnitude uses a sequential map starting at zero.
%
%   Masked vectors (all columns zero in the export) are shown in grey.
%
%   Each figure carries a legend with the nozzle pressure ratio of its
%   case, NPR = (p_gauge + P_AMB)/P_AMB. The files sort as 0001, 0001_1,
%   ..., 0001_5, which is the 3 to 8 bar order in GAUGE_BAR.

clear; clc

%% Settings ---------------------------------------------------------------
CC_DIR    = fullfile(fileparts(mfilename('fullpath')), 'CC results');   % next to this script
PATTERN   = 'BOS_3cm_12x12*.csv';
UNITS     = 'px';        % 'px' or 'mm'
pxPerMm   = 64.154;      % 3 cm runs, same value as PLOT_DISPLACEMENT_HISTOGRAM
PCTL      = 99.5;        % percentile that sets the colour limits
SAVE_FIGS = true;        % true writes a PNG per figure next to this script
MAGNITUDE_ONLY = true;   % true plots only the displacement magnitude
GAUGE_BAR = 3:8;         % tap gauge pressure of each file, in sorted file order [bar]
P_AMB     = 1.013;       % ambient pressure [bar]
NPR_LEGEND = false;      % true writes the case NPR in a legend on each plot
FONT_SIZE = 22;          % axes, colorbar and legend; sized for a 1/3 text width subfigure
FIG_SIZE  = [560 340];   % figure width and height per panel [px]

%% Read -------------------------------------------------------------------
files = dir(fullfile(CC_DIR, PATTERN));
if isempty(files)
    error('No files matching %s in:\n    %s', PATTERN, CC_DIR);
end
[~, order] = sort({files.name});
files = files(order);
nF    = numel(files);
if numel(GAUGE_BAR) ~= nF
    error('GAUGE_BAR has %d entries but %d files were found.', numel(GAUGE_BAR), nF);
end
NPR = (GAUGE_BAR + P_AMB) / P_AMB;

switch UNITS
    case 'px', scale = pxPerMm;
    case 'mm', scale = 1;
    otherwise, error('UNITS must be ''px'' or ''mm''.');
end

F = struct('name', {}, 'xv', {}, 'yv', {}, 'dx', {}, 'dy', {}, 'dm', {});
for k = 1:nF
    fprintf('reading %d/%d  %s\n', k, nF, files(k).name);
    D = readmatrix(fullfile(CC_DIR, files(k).name), ...
                   'Delimiter', ';', 'NumHeaderLines', 1);

    % Rebuild the regular grid from the vector list
    [xv, ~, ix] = unique(D(:,1));
    [yv, ~, iy] = unique(D(:,2));
    idx  = sub2ind([numel(yv) numel(xv)], iy, ix);
    mask = D(:,5) == 0;                         % masked vectors

    F(k).name = files(k).name;
    F(k).xv   = xv;
    F(k).yv   = yv;
    F(k).dx   = toGrid(D(:,3) * scale, idx, mask, numel(yv), numel(xv));
    F(k).dy   = toGrid(D(:,4) * scale, idx, mask, numel(yv), numel(xv));
    F(k).dm   = toGrid(D(:,5) * scale, idx, mask, numel(yv), numel(xv));
end

% The pressure order should show as a growing displacement
fprintf('\n%-28s %8s %6s %14s\n', 'file', 'p [barg]', 'NPR', 'P99 |d| [px]');
for k = 1:nF
    fprintf('%-28s %8d %6.2f %14.3f\n', F(k).name, GAUGE_BAR(k), NPR(k), ...
            percentile(F(k).dm(~isnan(F(k).dm)) * pxPerMm / scale, 99));
end

%% Shared colour limits ---------------------------------------------------
allDx = cell2mat(arrayfun(@(s) s.dx(~isnan(s.dx)), F(:), 'UniformOutput', false));
allDy = cell2mat(arrayfun(@(s) s.dy(~isnan(s.dy)), F(:), 'UniformOutput', false));
allDm = cell2mat(arrayfun(@(s) s.dm(~isnan(s.dm)), F(:), 'UniformOutput', false));

limX = niceCeil(percentile(abs(allDx), PCTL));
limY = niceCeil(percentile(abs(allDy), PCTL));
limM = niceCeil(percentile(allDm,      PCTL));

fprintf('\ncolour limits (%g th percentile, %s)\n', PCTL, UNITS);
fprintf('  x-displacement  : +/- %g\n', limX);
fprintf('  y-displacement  : +/- %g\n', limY);
fprintf('  magnitude       : 0 to %g\n\n', limM);

%% Plot -------------------------------------------------------------------
cmapDiv = divergingMap(256);
cmapSeq = parula(256);
nanGrey = [0.85 0.85 0.85];

for k = 1:nF
    panels = { F(k).dx, [-limX limX], cmapDiv, 'x-displacement'
               F(k).dy, [-limY limY], cmapDiv, 'y-displacement'
               F(k).dm, [0 limM],     cmapSeq, 'Displacement' };
    if MAGNITUDE_ONLY
        panels = panels(3,:);
    end
    nP = size(panels, 1);

    fig = figure('Color', 'w', 'Position', [60+25*k 60+25*k FIG_SIZE(1)*nP FIG_SIZE(2)], ...
                 'Name', F(k).name, 'NumberTitle', 'off');
    t = tiledlayout(fig, 1, nP, 'TileSpacing', 'compact', 'Padding', 'tight');

    % Crop to the region that holds valid vectors
    valid = ~isnan(F(k).dm);
    xr = F(k).xv([find(any(valid,1),1) find(any(valid,1),1,'last')]);
    yr = F(k).yv([find(any(valid,2),1) find(any(valid,2),1,'last')]);

    for p = 1:nP
        ax = nexttile(t);
        G  = panels{p,1};
        imagesc(ax, F(k).xv, F(k).yv, G, 'AlphaData', ~isnan(G));
        set(ax, 'YDir', 'normal', 'Color', nanGrey, 'FontSize', FONT_SIZE);
        axis(ax, 'image');
        xlim(ax, xr); ylim(ax, yr);
        clim(ax, panels{p,2});
        colormap(ax, panels{p,3});
        cb = colorbar(ax);
        cb.FontSize = FONT_SIZE;
        cb.Label.String = sprintf('%s [%s]', panels{p,4}, UNITS);
        cb.Label.FontSize = FONT_SIZE;
        xlabel(ax, 'x [mm]');
        ylabel(ax, 'y [mm]');
        ax.Toolbar.Visible = 'off';

        % Legend with the case NPR; the invisible line only carries the label
        if NPR_LEGEND
            hold(ax, 'on');
            h = plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none');
            hold(ax, 'off');
            lgd = legend(ax, h, sprintf('NPR = %.2f', NPR(k)), ...
                         'Location', 'southwest', 'FontSize', FONT_SIZE);
            lgd.ItemTokenSize = [1 1];
            lgd.Color = 'w';
        end
    end

    if SAVE_FIGS
        here = fileparts(mfilename('fullpath'));
        exportgraphics(fig, fullfile(here, strrep(F(k).name, '.csv', '_displacement.png')), ...
                       'Resolution', 300);
    end
end

%% Local functions --------------------------------------------------------
function G = toGrid(v, idx, mask, ny, nx)
    v(mask) = NaN;
    G = NaN(ny, nx);
    G(idx) = v;
end

function p = percentile(v, pct)
    % Nearest-rank percentile, so the script needs no toolbox
    v = sort(v(:));
    p = v(max(1, min(numel(v), ceil(pct/100 * numel(v)))));
end

function y = niceCeil(x)
    % Round up to 1, 1.5, 2, 2.5 or 5 times a power of ten
    if x <= 0, y = 1; return; end
    e = 10^floor(log10(x));
    steps = [1 1.5 2 2.5 5 10];
    y = steps(find(steps*e >= x, 1)) * e;
end

function c = divergingMap(n)
    % Blue - white - red, symmetric about the centre
    anchors = [0.230 0.299 0.754
               1.000 1.000 1.000
               0.706 0.016 0.150];
    c = interp1([0 0.5 1], anchors, linspace(0, 1, n));
end
