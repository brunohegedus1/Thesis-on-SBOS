% PLOT_DISPLACEMENT_HISTOGRAM
%
%   Reads the DaVis vector exports (.csv) in CC_DIR and plots the
%   displacement histograms, one figure per file.
%
%   The export stores displacements in millimetres in the object plane, so
%   they are converted to pixels here. Peak locking is a pixel-domain effect
%   and only shows up once the data is in pixels.
%
%   Two separate corrections are applied before the histograms are built.
%
%   1. BACKGROUND OFFSET (SUBTRACT_BG). The near-uniform displacement vector
%   that an image shift adds everywhere is subtracted from every vector. It
%   is the median of U and V inside BG_WINDOW, a rectangle of quiescent air,
%   computed by BOS_BACKGROUND_OFFSET, the same function BOS_PIPELINE uses
%   in its step 1. Leave BG_WINDOW empty to use the same default window.
%   Masked vectors are identified from the RAW export (both components
%   exactly zero) before the offset is subtracted.
%
%   2. JET REGION (REMOVE_BG). Only jet vectors are kept, so the histograms
%   show the jet only. The jet region is found on a smoothed |d| field of
%   the offset-corrected vectors: vectors
%   count as jet where the smoothed magnitude exceeds the background level
%   (its median) by N_SIGMA robust standard deviations. Only the connected
%   region holding the strongest signal is kept, so isolated noise patches
%   drop out. The region is then grown by GROW_MM so the weak vectors at its
%   edge stay in. Detection works on the smoothed field, not on each vector,
%   so small displacements inside the jet are kept. A band of EDGE_MM along
%   every mask edge is dropped first, because DaVis returns outliers there.
%   Set SHOW_MASK to check what was removed.
%
%   Top row     displacement distribution, after both corrections.
%   Bottom row  fractional part of the displacement. This is the standard
%               peak-locking check: it should be flat. Spikes at 0 and 1 with
%               a dip at 0.5 mean displacements are being pulled onto whole
%               pixel values.
%
%   The bottom row also uses the offset-corrected displacements. Peak
%   locking acts on the RAW correlation result, so with SUBTRACT_BG on, any
%   locking spikes sit shifted from 0 and 1 by the removed offset (printed
%   in px for each file). Set SUBTRACT_BG = false to see them at whole
%   pixels.
%
%   Read the warning the script prints. The bottom row is only meaningful
%   when a decent share of vectors exceed one pixel. If nearly everything is
%   sub-pixel, the fractional part piles up at 0 and 1 simply because the
%   distribution is centred on zero and narrower than a pixel, which looks
%   like locking but is not.

clear; clc

%% Settings ---------------------------------------------------------------
CC_DIR  = fullfile(fileparts(mfilename('fullpath')), 'CC results');   % next to this script
PATTERN = 'BOS_3cm_12x12*.csv';   % use a full file name to plot a single file
%pxPerMm = 72.43;      % use 64.154 for the 3 cm runs
pxPerMm = 64.154;

SUBTRACT_BG = true;    % subtract the background offset vector from every vector
BG_WINDOW   = [];      % [x1 y1; x2 y2] in mm; [] = BOS_BACKGROUND_OFFSET default,
                       % the same window BOS_PIPELINE uses

REMOVE_BG = true;      % keep only the jet region; false keeps every unmasked vector
SMOOTH_MM = 1.0;       % averaging window for the jet detection [mm]
N_SIGMA   = 5;         % jet threshold above the background [robust std]
GROW_MM   = 0.5;       % grow the detected jet region by this distance [mm]
EDGE_MM   = 0.3;       % band dropped along every mask edge [mm]
SHOW_MASK = true;      % extra figure with the jet region outlined

files = dir(fullfile(CC_DIR, PATTERN));
if isempty(files)
    error('No files matching %s in:\n    %s', PATTERN, CC_DIR);
end
[~, order] = sort({files.name});
files = files(order);

for k = 1:numel(files)
    file = files(k).name;

    %% Read ---------------------------------------------------------------
    D  = readmatrix(fullfile(CC_DIR, file), 'Delimiter', ';', 'NumHeaderLines', 1);
    dx = D(:,3);
    dy = D(:,4);

    keep = dx ~= 0 | dy ~= 0;          % drop masked vectors (raw zeros)

    %% Subtract background offset ------------------------------------------
    % Same window and estimator as step 1 of BOS_PIPELINE, in mm.
    if SUBTRACT_BG
        [u0, v0, nBg, bgWin] = bos_background_offset(D(keep,1), D(keep,2), ...
                                   dx(keep), dy(keep), BG_WINDOW, 'median');
        dx = dx - u0;
        dy = dy - v0;
    end

    dx   = dx(keep) * pxPerMm;         % mm -> px
    dy   = dy(keep) * pxPerMm;
    d    = hypot(dx, dy);
    nAll = numel(d);

    %% Remove background --------------------------------------------------
    % Rebuild the regular grid so the field can be smoothed
    [xv, ~, ix] = unique(D(:,1));
    [yv, ~, iy] = unique(D(:,2));
    idx = sub2ind([numel(yv) numel(xv)], iy(keep), ix(keep));
    h   = mean(diff(xv));              % vector spacing [mm]

    M = NaN(numel(yv), numel(xv));
    M(idx) = d;

    if REMOVE_BG
        % Drop a band along every mask edge, where DaVis returns outliers
        valid = conv2(double(~isnan(M)), disk(EDGE_MM / h), 'same') ...
                >= nnz(disk(EDGE_MM / h)) - 0.5;
        Mv = M;
        Mv(~valid) = NaN;

        n = max(1, round(SMOOTH_MM / h));
        S = movmean(movmean(Mv, n, 1, 'omitnan'), n, 2, 'omitnan');
        S(~valid) = NaN;

        bgLevel = median(S(:), 'omitnan');
        bgStd   = 1.4826 * median(abs(S(:) - bgLevel), 'omitnan');
        above   = S > bgLevel + N_SIGMA * bgStd;

        % Keep only the connected region that holds the strongest signal,
        % which is the jet; isolated noise patches are dropped
        Sa = S; Sa(~above) = -Inf;
        [~, seed] = max(Sa(:));
        jetGrid = false(size(above));
        jetGrid(seed) = true;
        while true
            grown = conv2(double(jetGrid), ones(9), 'same') > 0 & above;
            if isequal(grown, jetGrid), break; end
            jetGrid = grown;
        end

        jetGrid = conv2(double(jetGrid), disk(GROW_MM / h), 'same') > 0 & valid;

        jet = jetGrid(idx);
        dx  = dx(jet);
        dy  = dy(jet);
        d   = d(jet);
    end

    fprintf('file        : %s\n', file);
    if SUBTRACT_BG
        xw = sort(bgWin(:,1));  yw = sort(bgWin(:,2));
        fprintf('bg offset   : %+.3f px (x)   %+.3f px (y), subtracted; median of %d vectors in x [%.2f, %.2f], y [%.2f, %.2f] mm\n', ...
                u0*pxPerMm, v0*pxPerMm, nBg, xw(1), xw(2), yw(1), yw(2));
    end
    if REMOVE_BG
        fprintf('jet thresh. : background |d| %.3f px, robust std %.3f px (smoothed |d|)\n', bgLevel, bgStd);
        fprintf('vectors     : %d of %d kept as jet (%.1f %%)\n', numel(d), nAll, 100*numel(d)/nAll);
    else
        fprintf('vectors     : %d\n', numel(d));
    end
    fprintf('mean offset : %+.3f px (x)   %+.3f px (y)\n', mean(dx), mean(dy));
    fprintf('spread, std : %.3f px (x)   %.3f px (y)\n', std(dx), std(dy));
    fprintf('median |d|  : %.3f px\n', median(d));
    fprintf('max         : %.2f px\n\n', max(d));

    % What decides whether the sub-pixel panels mean anything is the SPREAD,
    % not the magnitude. A field sitting at a uniform offset of several
    % pixels still occupies a single pixel bin if its scatter is small, and
    % folding it by mod 1 then reproduces its own peak rather than revealing
    % locking.
    spread = max(std(dx), std(dy));
    if spread < 0.5
        fprintf(['WARNING: the displacement spread is only %.2f px.\n' ...
                 'The distribution is narrower than one pixel, so the sub-pixel\n' ...
                 'panels below reproduce its own shape and cannot show locking.\n' ...
                 'Treat them as uninformative for this file.\n\n'], spread);
    end

    % A large background offset is a shift between the reference and the
    % measurement image rather than flow. With SUBTRACT_BG on it has been
    % removed, and the sub-pixel panels are shifted by it; without, the
    % mean displacement is the best available estimate of it.
    if SUBTRACT_BG
        bulk = [u0 v0] * pxPerMm;
    else
        bulk = [mean(dx) mean(dy)];
    end
    if max(abs(bulk)) > 0.3
        fprintf(['NOTE: the background offset is %+.2f, %+.2f px. A bulk offset this\n' ...
                 'size usually indicates reference-to-measurement misalignment.\n\n'], ...
                 bulk(1), bulk(2));
    end

    %% Plot ---------------------------------------------------------------
    if REMOVE_BG && SHOW_MASK
        figure('Color', 'w', 'Position', [60+25*k 60+25*k 800 450], ...
               'Name', [file ' - jet region'], 'NumberTitle', 'off');
        imagesc(xv, yv, M, 'AlphaData', ~isnan(M));
        set(gca, 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
        axis image; clim([0 1.5]); colorbar
        hold on
        contour(xv, yv, double(jetGrid), [0.5 0.5], 'r', 'LineWidth', 1.2);
        hold off
        xlabel('x [mm]'); ylabel('y [mm]')
        title('|d| [px], red = jet region kept for the histograms')
    end

    figure('Color', 'w', 'Position', [100+25*k 100+25*k 950 620], ...
           'Name', file, 'NumberTitle', 'off');
    sgtitle(strrep(file, '.csv', ''), 'Interpreter', 'none');

    subplot(2,2,1)
    histogram(dx, 200)
    xlabel('x-displacement [px]'); ylabel('count')
    title('x-displacement'); grid on

    subplot(2,2,2)
    histogram(dy, 200)
    xlabel('y-displacement [px]'); ylabel('count')
    title('y-displacement'); grid on

    subplot(2,2,3)
    histogram(mod(dx,1), 0:0.02:1, 'Normalization', 'pdf', 'EdgeColor', 'none')
    hold on; yline(1, 'r--', 'LineWidth', 1.5)
    xlabel('fractional part of x [px]'); ylabel('density')
    title('Sub-pixel histogram, x   (flat = no peak locking)'); grid on

    subplot(2,2,4)
    histogram(mod(dy,1), 0:0.02:1, 'Normalization', 'pdf', 'EdgeColor', 'none')
    hold on; yline(1, 'r--', 'LineWidth', 1.5)
    xlabel('fractional part of y [px]'); ylabel('density')
    title('Sub-pixel histogram, y   (flat = no peak locking)'); grid on
end

%% Local functions --------------------------------------------------------
function k = disk(r)
    % Circular structuring element of radius r grid cells
    r = max(0, round(r));
    [u, v] = meshgrid(-r:r);
    k = double(u.^2 + v.^2 <= r^2);
end
