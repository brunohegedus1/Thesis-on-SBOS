%% compute_aerospike_density_gradient.m
% Computes the density gradient of the aerospike flow from the BOS
% displacements exported by DaVis (SBOS1.csv ... SBOSn.csv), assuming a
% planar flow. Based on run_bos_gradients_unfolded.m.
%
% METHOD
% The measured displacement D is converted to a deflection angle with the
% effective lever arm of the setup:
%
%     l_eff = l + G + W/2                         [mm]
%     eps   = D / l_eff                           [rad]
%
% where l is the defocus distance, G the gap between the test section and
% the background and W the width of the planar object along the line of
% sight. The planar relation of the thesis (Eq. 2.22)
%
%     eps_x = (K*W/n0) * drho/dx
%
% is then inverted for each component:
%
%     drho/dx = n0 * eps_x / (K*W)
%     drho/dy = n0 * eps_y / (K*W)
%
% K is the Gladstone-Dale constant (the thesis writes it G; here G is the
% gap distance). With W in mm the gradient comes out in kg m^-3 per mm.
%
% ASSUMPTION
% The planar relation assumes the gradient does not vary along the line of
% sight over the width W. Where the real flow is not two-dimensional, for
% example near the end walls, the result is the equivalent planar
% gradient: the uniform gradient over W that bends a ray by the measured
% amount.
%
% The script reads each selected run, rotates it so the top edge of the
% measurement region is horizontal, centres (0,0) in the image, computes
% the gradients and plots d(rho)/dx, d(rho)/dy and, optionally, |grad(rho)|.

clear; clc;

%% Settings
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'Aerospike Data');
figDir    = fullfile(scriptDir, 'Figures', 'Density gradient');

% Runs to process, by number: [4 5 6 7] processes SBOS4 to SBOS7.
% Leave empty ([]) to process every CSV file in dataDir.
runsToProcess = [5];

% --- optical geometry [mm] ---------------------------------------------
l = 55;     % defocus distance
G = 30;     % gap between test section and background
W = 150;    % width of the planar object along the line of sight
l_eff = l + G + W/2;

%--- INSIDE TS ---
%l_eff = 126.292079; % For dataset 2
%l_eff = 123.485602; % For points 8-11 


% --- gas properties (air) ----------------------------------------------
K       = 2.26e-4;          % Gladstone-Dale constant [m^3/kg]
rho_amb = 1.2;              % ambient density [kg/m^3]
n0      = 1 + K * rho_amb;  % ambient refractive index

% --- processing -----------------------------------------------------------
alignToTopEdge = true;  % rotate image and vectors so the top edge is horizontal
signFlip       = false; % flip the sign of the displacements
smoothWindow   = 0;     % moving-average window on the displacements, in
                        % grid points (0 = off)

% --- printed text size ------------------------------------------------------
% The 'separate' figures are included in the report as subfigures of
% 0.49\textwidth. Drawing them at exactly that width and setting every label
% to fontSize means LaTeX does not rescale the image, so its letters print at
% fontSize points. textWidth is the text width of the report class
% (a4paper, hscale = 0.75, so 0.75*210 = 157.5 mm).
fontSize  = 8;              % pt, all text in the 'separate' figures
textWidth = 15.75;          % cm, \textwidth of the report class
widthFrac = 0.49;           % subfigure width as a fraction of \textwidth
figW      = widthFrac * textWidth;
figH      = 0.72 * figW;    % keeps the aspect of the former figures
printDpi  = 600;            % raster resolution of the exported jpg

% --- output ----------------------------------------------------------------
clipQuantile = 0.995;   % colour limits clip the top 0.5 % of |values|
nPeakAvg     = 100;     % the printed peak is also averaged over the
                        % nPeakAvg largest |values|, which is less
                        % sensitive to a single outlying point
plotMagnitude = false;   % false plots only d(rho)/dx and d(rho)/dy
figureLayout = 'combined';  % 'combined': one figure per run, panels side by side
                            % 'separate': one figure per gradient
saveFigures  = false; % write PNG and FIG files to figDir
saveResults  = false;   % write a MAT file per run to figDir

if (saveFigures || saveResults) && ~exist(figDir, 'dir')
    mkdir(figDir);
end

uLab = 'kg m^{-3} mm^{-1}';

%% Find the files, sorted by run number (SBOS2 before SBOS10)
files = dir(fullfile(dataDir, '*.csv'));
if isempty(files)
    error('No CSV files found in %s', dataDir);
end
runNum = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), {files.name});
[runNum, order] = sort(runNum);
files = files(order);

if ~isempty(runsToProcess)
    missing = setdiff(runsToProcess, runNum);
    if ~isempty(missing)
        warning('No CSV file found for run(s): %s', num2str(missing));
    end
    files = files(ismember(runNum, runsToProcess));
    if isempty(files)
        error('None of the selected runs exist in %s', dataDir);
    end
end

fprintf('l_eff = l + G + W/2 = %g + %g + %g/2 = %g mm\n', l, G, W, l_eff);
fprintf('K = %.3e m^3/kg, n0 = %.6f\n\n', K, n0);

%% Loop over runs
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
        fprintf('  top edge at %.3f deg, rotated by %.3f deg\n', ...
            rad2deg(theta), -rad2deg(theta));
    end

    % Centre (0,0) in the image
    valid = ~isnan(F.mag);
    xr = F.x(any(valid, 1));
    yr = F.y(any(valid, 2));
    x0 = (xr(1) + xr(end)) / 2;
    y0 = (yr(1) + yr(end)) / 2;
    F.x = F.x - x0;
    F.y = F.y - y0;
    axLim = [xr(1) xr(end) yr(1) yr(end)] - [x0 x0 y0 y0];

    %% Displacement -> deflection -> density gradient
    Dx = F.dx;
    Dy = F.dy;
    if signFlip
        Dx = -Dx;
        Dy = -Dy;
    end
    if smoothWindow > 1
        Dx = smoothNaN(Dx, smoothWindow);
        Dy = smoothNaN(Dy, smoothWindow);
    end

    eps_x = Dx / l_eff;                          % [rad]
    eps_y = Dy / l_eff;

    R.x        = F.x;
    R.y        = F.y;
    R.drho_dx  = n0 * eps_x / (K * W);           % [kg m^-3 mm^-1]
    R.drho_dy  = n0 * eps_y / (K * W);
    R.grad_mag = hypot(R.drho_dx, R.drho_dy);
    % Standard uncertainty of each component, from the DaVis uncertainty
    % of the displacement. It covers the correlation only, not l_eff or W.
    R.u_drho_dx = n0 * (F.uncX / l_eff) / (K * W);
    R.u_drho_dy = n0 * (F.uncY / l_eff) / (K * W);
    R.settings  = struct('l', l, 'G', G, 'W', W, 'l_eff', l_eff, 'K', K, ...
        'n0', n0, 'signFlip', signFlip, 'smoothWindow', smoothWindow);

    fprintf('  peak |deflection| : %.3e rad\n', max(hypot(eps_x(:), eps_y(:))));
    % The mean is over the signed values, so it shows any offset of the
    % field; the RMS is about zero, not about the mean.
    fprintf('  peak |drho/dx|    : %8.4f   top%d %8.4f   mean %+8.4f   RMS %8.4f   median u %8.4f\n', ...
        max(abs(R.drho_dx(:))), nPeakAvg, topMean(R.drho_dx, nPeakAvg), ...
        mean(R.drho_dx(:), 'omitnan'), ...
        rmsNaN(R.drho_dx), median(R.u_drho_dx(:), 'omitnan'));
    fprintf('  peak |drho/dy|    : %8.4f   top%d %8.4f   mean %+8.4f   RMS %8.4f   median u %8.4f\n', ...
        max(abs(R.drho_dy(:))), nPeakAvg, topMean(R.drho_dy, nPeakAvg), ...
        mean(R.drho_dy(:), 'omitnan'), ...
        rmsNaN(R.drho_dy), median(R.u_drho_dy(:), 'omitnan'));
    fprintf('  peak |grad rho|   : %8.4f   top%d %8.4f   mean %+8.4f   RMS %8.4f   [kg m^-3 mm^-1]\n', ...
        max(R.grad_mag(:)), nPeakAvg, topMean(R.grad_mag, nPeakAvg), ...
        mean(R.grad_mag(:), 'omitnan'), rmsNaN(R.grad_mag));

    %% Figure
    % One entry per panel: field, title, file name suffix, colormap and
    % whether the colour limits are symmetric about zero.
    panels = { ...
        R.drho_dx,  '\partial\rho/\partial x', 'dx',  @blueWhiteRed, true; ...
        R.drho_dy,  '\partial\rho/\partial y', 'dy',  @blueWhiteRed, true};
    if plotMagnitude
        panels(end+1, :) = ...
            {R.grad_mag, '|\nabla\rho|',       'mag', @parula,       false};
    end
    nPanels = size(panels, 1);

    switch figureLayout
        case 'combined'
            fig = figure('Name', runName, 'Color', 'w', ...
                'Position', [60 60 500 * nPanels 480]);
            tl  = tiledlayout(fig, 1, nPanels, 'TileSpacing', 'compact', ...
                'Padding', 'compact');
            axs = gobjects(1, nPanels);
            for p = 1:nPanels
                axs(p) = nexttile(tl);
                drawPanel(axs(p), R.x, R.y, panels(p, :), uLab, clipQuantile);
            end
            linkaxes(axs);
            axis(axs(1), axLim);
            if saveFigures
                saveFigure(fig, figDir, [runName '_density_gradient']);
            end

        case 'separate'
            % Drawn at the size it occupies in the report, so the exported
            % file needs no rescaling and its text keeps fontSize points.
            for p = 1:nPanels
                fig = figure('Name', [runName ' ' panels{p, 3}], 'Color', 'w', ...
                    'Units', 'centimeters', 'Position', [4 4 figW figH], ...
                    'PaperPositionMode', 'auto');
                ax = axes(fig);
                drawPanel(ax, R.x, R.y, panels(p, :), uLab, clipQuantile);
                axis(ax, axLim);
                setFontSize(fig, fontSize);
                if saveFigures
                    printFigure(fig, figDir, ...
                        [runName '_density_gradient_' panels{p, 3}], printDpi);
                end
            end

        otherwise
            error('figureLayout must be ''combined'' or ''separate''.');
    end

    %% Save
    if saveResults
        save(fullfile(figDir, [runName '_density_gradient.mat']), '-struct', 'R');
    end
    fprintf('\n');
end

fprintf('Done. Processed %d runs.\n', numel(files));

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

    function Gr = toGrid(v)
        v(masked) = NaN;
        Gr = nan(numel(yu), numel(xu));
        Gr(idx) = v;
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

    [X, Y] = meshgrid(F.x, F.y);
    valid  = ~isnan(F.mag);
    xr = c * X(valid) + s * Y(valid);
    yr = -s * X(valid) + c * Y(valid);

    h = median(diff(F.x));
    R.x = min(xr):h:max(xr);
    R.y = min(yr):h:max(yr);

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

function Z = smoothNaN(Z, w)
% Moving average over a w-by-w window that ignores NaN and keeps the
% masked region masked.
    mask = isnan(Z);
    Z = smoothdata2(Z, 'movmean', w, 'omitmissing');
    Z(mask) = NaN;
end

function drawPanel(ax, x, y, panel, unitLabel, clipQuantile)
% Draws one gradient panel with its colormap and colour limits.
% PANEL is one row of the panels list: {field, title, suffix, cmap, symmetric}.
    [Z, name, ~, cmap, symmetric] = panel{:};
    plotField(ax, x, y, Z, name, unitLabel);
    lim = symLimit(Z, clipQuantile);
    if symmetric
        clim(ax, [-lim lim]);
    else
        clim(ax, [0 lim]);
    end
    colormap(ax, cmap(256));
end

function plotField(ax, x, y, Z, name, unitLabel)
% Draws a field with the masked area transparent and equal axis scaling.
    imagesc(ax, x, y, Z, 'AlphaData', ~isnan(Z));
    set(ax, 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
    axis(ax, 'image');
    xlabel(ax, 'x [mm]');
    ylabel(ax, 'y [mm]');
    title(ax, name);
    cb = colorbar(ax);
    cb.Label.String = sprintf('%s [%s]', name, unitLabel);
end

function setFontSize(fig, pt)
% One size for every letter: tick labels, axis labels, title and the
% colorbar label. MATLAB draws axis labels and titles larger than the axes
% font size by default, so both multipliers are reset first.
    axs = findall(fig, 'Type', 'axes');
    set(axs, 'FontSize', pt, 'LabelFontSizeMultiplier', 1, ...
             'TitleFontSizeMultiplier', 1);
    set(findall(fig, '-property', 'FontSize'), 'FontSize', pt);
end

function printFigure(fig, figDir, baseName, dpi)
% PRINT, not EXPORTGRAPHICS: the latter crops the margins, which makes the
% file narrower than the panel. LaTeX would then stretch it back to
% 0.49\textwidth and the letters would print larger than fontSize.
    axs = findall(fig, 'Type', 'axes');
    for i = 1:numel(axs)
        axs(i).Toolbar.Visible = 'off';
    end
    % InvertHardcopy off keeps the grey axes background that shows through
    % the masked region; PRINT whitens backgrounds by default. '-image'
    % forces the raster renderer, which honours the AlphaData of the mask.
    set(fig, 'InvertHardcopy', 'off');
    print(fig, fullfile(figDir, [baseName '.jpg']), '-djpeg95', '-image', ...
          sprintf('-r%d', dpi));
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

function v = topMean(A, n)
% Mean of the n largest |values| of A, ignoring NaN. A peak read from a
% single grid point is set by noise; averaging the n highest points gives
% a peak that a single outlier cannot move.
    a = sort(abs(A(isfinite(A))), 'descend');
    if isempty(a)
        v = NaN;
        return
    end
    v = mean(a(1:min(n, numel(a))));
end

function v = rmsNaN(A)
    a = A(isfinite(A));
    if isempty(a), v = NaN; else, v = sqrt(mean(a.^2)); end
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
