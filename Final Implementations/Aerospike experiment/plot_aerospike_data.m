%% plot_aerospike_data.m
% Plots the BOS displacement fields exported from DaVis for the aerospike
% runs (SBOS1.csv ... SBOSn.csv).
%
% Each CSV holds a regular grid in long format, separated by ';':
%   x [mm]; y [mm]; x-displacement [mm]; y-displacement [mm];
%   Displacement [mm]; Uncertainty V [mm]; Uncertainty Vx [mm]; Uncertainty Vy [mm]
% Grid points outside the processing mask are exported as all zeros. The
% script sets them to NaN so they render as blank.
%
% In the settings you choose which runs and which variables to plot, and
% whether the variables share one figure per run or each get their own
% figure. Figures are saved as PNG and FIG in the "Figures" folder next
% to this script.

clear; clc;

%% Settings
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'Aerospike Data');
figDir    = fullfile(scriptDir, 'Figures');

% Runs to plot, by number: [1 4 10] plots SBOS1, SBOS4 and SBOS10.
% Leave empty ([]) to plot every CSV file in dataDir.
runsToPlot = [8 9 10 11];

% Variables to plot, in the order they appear. Available names:
%   'dx'   x-displacement            'unc'   Uncertainty V
%   'dy'   y-displacement            'uncX'  Uncertainty Vx
%   'mag'  Displacement magnitude    'uncY'  Uncertainty Vy
variablesToPlot = {'mag'};

% 'combined': one figure per run with a panel for each variable.
% 'separate': one figure per run and variable.
figureLayout = 'combined';

saveFigures    = false;  % write PNG and FIG files to figDir
clipQuantile   = 0.995;  % colour limits clip the top 0.5 % of |values|
quiverStep     = 0;      % vector every quiverStep grid points on 'mag' (0 = off)
quiverScale    = 0;      % arrow length multiplier for the quiver plot
alignToTopEdge = true;   % rotate image and vectors so the top edge is horizontal

%% Variable definitions
% name, colour bar label, colormap, symmetric colour limits around zero
varDefs = struct( ...
    'dx',   struct('label', 'x-displacement [mm]',         'cmap', @blueWhiteRed, 'symmetric', true), ...
    'dy',   struct('label', 'y-displacement [mm]',         'cmap', @blueWhiteRed, 'symmetric', true), ...
    'mag',  struct('label', 'Displacement magnitude [mm]', 'cmap', @parula,       'symmetric', false), ...
    'unc',  struct('label', 'Uncertainty V [mm]',          'cmap', @hot,          'symmetric', false), ...
    'uncX', struct('label', 'Uncertainty Vx [mm]',         'cmap', @hot,          'symmetric', false), ...
    'uncY', struct('label', 'Uncertainty Vy [mm]',         'cmap', @hot,          'symmetric', false));

unknown = setdiff(variablesToPlot, fieldnames(varDefs));
if ~isempty(unknown)
    error('Unknown variable(s): %s. Choose from: %s', ...
        strjoin(unknown, ', '), strjoin(fieldnames(varDefs).', ', '));
end
if ~ismember(figureLayout, {'combined', 'separate'})
    error('figureLayout must be ''combined'' or ''separate''.');
end

if saveFigures && ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Find the files, sorted by run number (SBOS2 before SBOS10)
files = dir(fullfile(dataDir, '*.csv'));
if isempty(files)
    error('No CSV files found in %s', dataDir);
end
runNum = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), {files.name});
[runNum, order] = sort(runNum);
files = files(order);

if ~isempty(runsToPlot)
    missing = setdiff(runsToPlot, runNum);
    if ~isempty(missing)
        warning('No CSV file found for run(s): %s', num2str(missing));
    end
    files = files(ismember(runNum, runsToPlot));
    if isempty(files)
        error('None of the selected runs exist in %s', dataDir);
    end
end

%% Loop over files
for k = 1:numel(files)
    fname = files(k).name;
    [~, runName] = fileparts(fname);
    fprintf('Reading %s (%d/%d)\n', fname, k, numel(files));

    D = readmatrix(fullfile(dataDir, fname), ...
        'FileType', 'text', 'Delimiter', ';', 'NumHeaderLines', 1);

    F = csvToGrid(D);

    if alignToTopEdge
        theta = topEdgeAngle(F);
        F = rotateFields(F, theta);
        fprintf('  Top edge at %.3f deg, rotated by %.3f deg\n', ...
            rad2deg(theta), -rad2deg(theta));
    end

    % Shift the coordinates so (0,0) sits at the centre of the unmasked
    % data, which is also the centre of the plotted image
    valid = ~isnan(F.mag);
    xr = F.x(any(valid, 1));
    yr = F.y(any(valid, 2));
    x0 = (xr(1) + xr(end)) / 2;
    y0 = (yr(1) + yr(end)) / 2;
    F.x = F.x - x0;
    F.y = F.y - y0;

    % Axis limits: extent of the unmasked data
    axLim = [xr(1) xr(end) yr(1) yr(end)] - [x0 x0 y0 y0];

    nVar = numel(variablesToPlot);
    switch figureLayout
        case 'combined'
            fig = figure('Name', runName, 'Color', 'w', 'Position', [100 100 1400 900]);
            [nRows, nCols] = tileGrid(nVar);
            tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
            axs = gobjects(1, nVar);
            for v = 1:nVar
                axs(v) = nexttile(tl);
                drawVariable(axs(v), F, variablesToPlot{v}, varDefs, ...
                    clipQuantile, quiverStep, quiverScale);
            end
            linkaxes(axs);
            axis(axs(1), axLim);
            if saveFigures
                saveFigure(fig, figDir, runName);
            end

        case 'separate'
            for v = 1:nVar
                name = variablesToPlot{v};
                fig = figure('Name', [runName ' ' name], 'Color', 'w', ...
                    'Position', [100 100 900 650]);
                ax = axes(fig);
                drawVariable(ax, F, name, varDefs, ...
                    clipQuantile, quiverStep, quiverScale);
                axis(ax, axLim);
                if saveFigures
                    saveFigure(fig, figDir, [runName '_' name]);
                end
            end
    end
end

fprintf('Done. Plotted %d files.\n', numel(files));

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function F = csvToGrid(D)
% Converts the long-format DaVis export into 2-D arrays on the x-y grid.
% Rows where every value column is zero lie outside the mask and become NaN.
    [xu, ~, ix] = unique(D(:, 1));
    [yu, ~, iy] = unique(D(:, 2));
    idx = sub2ind([numel(yu) numel(xu)], iy, ix);

    masked = all(D(:, 3:end) == 0, 2);

    F.x = xu.';
    F.y = yu.';
    F.dx   = toGrid(D(:, 3));
    F.dy   = toGrid(D(:, 4));
    F.mag  = toGrid(D(:, 5));
    F.unc  = toGrid(D(:, 6));
    F.uncX = toGrid(D(:, 7));
    F.uncY = toGrid(D(:, 8));

    function G = toGrid(v)
        v(masked) = NaN;
        G = nan(numel(yu), numel(xu));
        G(idx) = v;
    end
end

function theta = topEdgeAngle(F)
% Returns the angle [rad] of the top edge of the unmasked region,
% measured counter-clockwise from the x-axis. For every column the top
% edge is the highest valid grid point. A line is fitted through these
% points, and points that do not lie on the straight part of the edge
% (for example the nozzle ramp) are rejected iteratively.
    valid = ~isnan(F.mag);
    cols  = find(any(valid, 1));
    if numel(cols) < 2
        theta = 0;
        return
    end
    xt = F.x(cols);
    yt = zeros(size(xt));
    for i = 1:numel(cols)
        yt(i) = F.y(find(valid(:, cols(i)), 1, 'last'));
    end

    h    = median(diff(F.y));          % grid spacing
    keep = true(size(xt));
    for iter = 1:10
        p   = polyfit(xt(keep), yt(keep), 1);
        res = yt - polyval(p, xt);
        s   = 1.4826 * median(abs(res(keep) - median(res(keep))));
        newKeep = abs(res) < max(3 * s, 2 * h);
        if isequal(newKeep, keep) || nnz(newKeep) < 2
            break
        end
        keep = newKeep;
    end
    theta = atan(p(1));
end

function R = rotateFields(F, theta)
% Rotates the fields by -theta about the origin, so an edge at angle theta
% becomes horizontal. The rotated fields are interpolated onto a new
% regular grid with the original spacing. The displacement vectors are
% rotated by the same angle so they stay consistent with the new axes.
% The uncertainty components are rotated assuming uncorrelated x and y
% errors: u'^2 = (c*ux)^2 + (s*uy)^2.
    c = cos(theta);
    s = sin(theta);

    % Extent of the valid data in the rotated frame
    [X, Y] = meshgrid(F.x, F.y);
    valid  = ~isnan(F.mag);
    xr = c * X(valid) + s * Y(valid);
    yr = -s * X(valid) + c * Y(valid);

    h = median(diff(F.x));
    R.x = min(xr):h:max(xr);
    R.y = min(yr):h:max(yr);

    % Position of each new grid point in the original frame
    [Xn, Yn] = meshgrid(R.x, R.y);
    Xs = c * Xn - s * Yn;
    Ys = s * Xn + c * Yn;

    interp = @(Z) interp2(F.x, F.y, Z, Xs, Ys, 'linear', NaN);
    dx = interp(F.dx);
    dy = interp(F.dy);
    ux = interp(F.uncX);
    uy = interp(F.uncY);

    R.dx   =  c * dx + s * dy;
    R.dy   = -s * dx + c * dy;
    R.mag  = interp(F.mag);
    R.unc  = interp(F.unc);
    R.uncX = sqrt((c * ux).^2 + (s * uy).^2);
    R.uncY = sqrt((s * ux).^2 + (c * uy).^2);
end

function drawVariable(ax, F, name, varDefs, clipQuantile, quiverStep, quiverScale)
% Draws one variable with its colormap and colour limits. The magnitude
% plot gets displacement vectors on top when quiverStep > 0.
    def = varDefs.(name);
    Z   = F.(name);
    plotField(ax, F.x, F.y, Z, def.label);
    lim = symLimit(Z, clipQuantile);
    if def.symmetric
        clim(ax, [-lim lim]);
    else
        clim(ax, [0 lim]);
    end
    colormap(ax, def.cmap(256));

    if strcmp(name, 'mag') && quiverStep > 0
        hold(ax, 'on');
        iy = 1:quiverStep:numel(F.y);
        ix = 1:quiverStep:numel(F.x);
        [Xq, Yq] = meshgrid(F.x(ix), F.y(iy));
        quiver(ax, Xq, Yq, F.dx(iy, ix), F.dy(iy, ix), quiverScale, ...
            'k', 'LineWidth', 0.6);
        hold(ax, 'off');
    end
end

function plotField(ax, x, y, Z, label)
% Draws a field with the masked area transparent and equal axis scaling.
    imagesc(ax, x, y, Z, 'AlphaData', ~isnan(Z));
    set(ax, 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
    axis(ax, 'image');
    xlabel(ax, 'x [mm]');
    ylabel(ax, 'y [mm]');
    title(ax, label);
    cb = colorbar(ax);
    cb.Label.String = label;
end

function [nRows, nCols] = tileGrid(n)
% Picks a tile layout for n panels: 1x1, 1x2, 1x3, 2x2, 2x3.
    nCols = min(n, 3);
    if n == 4
        nCols = 2;
    end
    nRows = ceil(n / nCols);
end

function saveFigure(fig, figDir, baseName)
% Hides the axes toolbars so they do not end up in the PNG.
    axs = findall(fig, 'Type', 'axes');
    for i = 1:numel(axs)
        axs(i).Toolbar.Visible = 'off';
    end
    exportgraphics(fig, fullfile(figDir, [baseName '.png']), 'Resolution', 200);
    savefig(fig, fullfile(figDir, [baseName '.fig']));
end

function lim = symLimit(Z, q)
% Returns the q-quantile of |Z|, ignoring NaN. Avoids the Statistics
% Toolbox dependency of prctile.
    v = sort(abs(Z(~isnan(Z))));
    if isempty(v)
        lim = 1;
        return
    end
    lim = v(max(1, round(q * numel(v))));
    if lim == 0
        lim = max(v);
    end
    if lim == 0
        lim = 1;
    end
end

function cmap = blueWhiteRed(n)
% Diverging colormap: blue for negative, white at zero, red for positive.
    h = floor(n / 2);
    down = [linspace(0.23, 1, h)' linspace(0.30, 1, h)' linspace(0.75, 1, h)'];
    up   = [linspace(1, 0.70, n - h)' linspace(1, 0.02, n - h)' linspace(1, 0.15, n - h)'];
    cmap = [down; up];
end
