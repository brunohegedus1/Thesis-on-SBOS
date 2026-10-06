%% SPECKLE_SIZE  Mean speckle size from a recorded intensity map.
%
% Method (Hamarova, Horvath, Smid & Hrabovsky, "A new approach for
% determination of a mean speckle size in simulated speckle pattern",
% Measurement 88 (2016) 271-277):
%
%   The mean speckle size is the width of the normalized autocorrelation
%   function of intensity,
%
%       r_I(dx,dy) = ( <I(x,y) I(x-dx,y-dy)> - <I>^2 ) / ( <I^2> - <I>^2 )
%
%   read as the lag at which r_I drops to 1/e (eq. 1-3 of the paper).
%
% Three estimators are computed and compared:
%   (A) 2D autocorrelation, profiles cut through the origin along x and y.
%       This is the conventional estimator (refs. [2,3] of the paper).
%   (B) Average 1D correlation method (eq. 4-5): one 1D ACF per row and per
%       column, threshold crossing found for each, then averaged.
%   (C) 1D correlation method of a long vector (Fig. 2b): all rows (columns)
%       concatenated into a single vector, one ACF per direction. This is
%       the estimator the paper recommends; it stays within ~1.5% of theory,
%       whereas (B) underestimates the size when few speckles are sampled.
%
% Optional decimation by a factor DELTA (eq. 6-7) is supported for (B), (C).
%
% Input file: a semicolon separated DaVis export with one header line. Both
% coordinate formats are accepted, see read_intensity_map.m:
%   x [pixel]; y [pixel]; Average [counts]; Number of data points [samples]
%   x [mm];    y [mm];    Average [];       Number of data points [samples]
% When the file carries millimetres, its sample spacing is adopted as the
% pixel pitch unless opt.pixelPitch_um has been set by hand.
%
% Output: printed summary, three figures, and speckle_size_results.mat
%
% Note on the illumination envelope. The recorded beam is not uniform, and
% that slow brightness variation is itself correlated over hundreds of
% pixels. Left in the data it lifts the whole ACF onto a positive pedestal
% (r_I stays near +0.08 out to a lag of 40 px) and widens the measured size.
% The script therefore divides the image by a Gaussian-smoothed copy of
% itself. The measured size depends on the smoothing width, for
% intensity_6.txt:
%
%   sigma [px]      0      25      50      75     100     120     150
%   a(1/e) [px]     10.99   9.14    9.86   10.06   10.17   10.24   10.34
%   ACF tail        +0.083   -       -     -0.005  +0.003  +0.009  +0.017
%
% (ACF tail = mean r_I over lags 51-60 px, far outside the correlation peak.)
%
% sigma must be large enough not to attenuate the speckle itself, and small
% enough to still flatten the beam. Too small and the contrast collapses
% (0.71 at sigma = 25 px against 0.90 raw) and the size comes out too low.
% sigma = 100 px, about ten speckle radii, is the default: it drives the
% tail to +0.003. The residual spread over sigma = 75 to 150 px is about
% +-1.5% on a. The diagnostic line in the output reports that tail; retune
% sigma if it is not close to zero.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc; close all;

est = speckle_estimators();   % shared routines, see speckle_estimators.m

%% ------------------------------------------------------------------ input
opt.dataFile      = 'I3_2.txt';  % intensity map, semicolon separated
opt.cacheFile     = 'I3_2.mat';  % parsed image cache (built on 1st run)

% Both sit next to this script, so they are resolved against its folder and
% the script runs from any working directory.
opt.folder    = fileparts(mfilename('fullpath'));
opt.dataFile  = fullfile(opt.folder, opt.dataFile);
opt.cacheFile = fullfile(opt.folder, opt.cacheFile);

opt.pixelPitch_um = NaN;    % camera pixel pitch [um]. NaN -> results in pixels
opt.roi           = [];     % [x0 y0 width height] in pixels, [] -> whole frame
opt.detrend       = true;   % remove the low-frequency illumination envelope
opt.detrendSigma  = 20;    % [px] Gaussian sigma of the envelope estimate.
                            % Keep it >> speckle size and << beam size. The
                            % result is sensitive to it: too small a sigma
                            % eats into the speckle itself and shrinks the
                            % measured size. Check the pedestal diagnostic
                            % printed below and see the note in the header.
opt.maxLag        = 60;     % [px] lag window kept for searching and plotting
opt.thresholds    = [1/exp(1), 0.5];   % ACF levels at which the size is read
opt.decimation    = 1;      % DELTA of eq. (6)-(7); 1 = every row and column
opt.savePlots     = true;

%% ------------------------------------------------------- load / build image
if isfile(opt.cacheFile)
    fprintf('Loading cached image from %s ...\n', opt.cacheFile);
    S = load(opt.cacheFile, 'I', 'pitch_um');
    I = S.I;
    if isfield(S, 'pitch_um'), filePitch_um = S.pitch_um; else, filePitch_um = NaN; end
else
    fprintf('Reading %s (5.5e6 rows, this takes a while) ...\n', opt.dataFile);
    [I, filePitch_um] = read_intensity_map(opt.dataFile);
    pitch_um = filePitch_um;                    
    save(opt.cacheFile, 'I', 'pitch_um', '-v7.3');
    fprintf('Cached parsed image to %s\n', opt.cacheFile);
end

% Exports that carry millimetre coordinates define their own scale. It is
% adopted only when no pitch has been set by hand, so an explicit value in
% the input block always wins.
if isnan(opt.pixelPitch_um) && ~isnan(filePitch_um)
    opt.pixelPitch_um = filePitch_um;
    fprintf('Sample pitch %.3f um taken from the coordinates in the file.\n', ...
            filePitch_um);
end

[ny, nx] = size(I);
fprintf('Image: %d x %d px, mean = %.2f counts, std = %.2f counts\n', ...
        nx, ny, mean(I(:)), std(I(:)));

if ~isempty(opt.roi)
    r = opt.roi;
    I = I(r(2)+(1:r(4)), r(1)+(1:r(3)));
    [ny, nx] = size(I);
    fprintf('ROI applied: %d x %d px\n', nx, ny);
end

%% -------------------------------------------- flatten illumination envelope
% A slowly varying beam profile is a deterministic, strongly correlated
% component. Left in, it adds a wide pedestal to r_I and biases the speckle
% size upwards. Dividing by a smoothed copy leaves the speckle fluctuation.
Iraw = I;
if opt.detrend
    env = est.gaussblur(I, opt.detrendSigma);
    env(env <= 0) = eps('single');
    I = I ./ env * mean(env(:));
    fprintf('Illumination flattened (Gaussian envelope, sigma = %g px).\n', ...
            opt.detrendSigma);
end

% Speckle contrast, a sanity check: ~1 for fully developed speckle.
contrast = std(I(:)) / mean(I(:));
fprintf('Speckle contrast std/mean = %.3f\n', contrast);

%% ---------------------------------------------------- (A) 2D autocorrelation
fprintf('\nComputing 2D autocorrelation ...\n');
[acf2, lagx, lagy] = est.acf2(I, opt.maxLag);

profX = acf2(lagy == 0, :);    % cut along x' (horizontal)
profY = acf2(:, lagx == 0).';  % cut along y' (vertical)

res.acf2D.x = est.widths(lagx, profX, opt.thresholds);
res.acf2D.y = est.widths(lagy, profY, opt.thresholds);

%% ------------------------------ (B), (C) 1D correlation methods of the paper
D = max(1, round(opt.decimation));
Idec_rows = I(1:D:end, :);     % keep every D-th row    (eq. 6)
Idec_cols = I(:, 1:D:end);     % keep every D-th column (eq. 7)

fprintf('Computing 1D correlation methods (decimation DELTA = %d) ...\n', D);

% (B) average 1D correlation method, eq. (4)-(5)
[res.avg1D.x, res.avg1D.xStd, nUsedX] = ...
    est.avg1d(Idec_rows,   opt.maxLag, opt.thresholds);   % rows -> x'
[res.avg1D.y, res.avg1D.yStd, nUsedY] = ...
    est.avg1d(Idec_cols.', opt.maxLag, opt.thresholds);   % cols -> y'

% (C) 1D correlation method of a long vector, Fig. 2b
vecX = reshape(Idec_rows.', [], 1);   % rows concatenated left to right
vecY = reshape(Idec_cols,   [], 1);   % columns concatenated top to bottom

% The segment lengths tell the estimator where the joints between rows and
% columns lie, so that the uncorrelated pairs straddling them are corrected.
[acfX, lag1] = est.acf1(vecX, opt.maxLag, size(Idec_rows, 2));
[acfY, ~   ] = est.acf1(vecY, opt.maxLag, size(Idec_cols, 1));

res.long1D.x = est.widths(lag1, acfX, opt.thresholds);
res.long1D.y = est.widths(lag1, acfY, opt.thresholds);

%% ----------------------------------------------------------------- report
px = opt.pixelPitch_um;
useUm = ~isnan(px);
if useUm, unit = 'um'; k = px; else, unit = 'px'; k = 1; end

fprintf('\n===================== MEAN SPECKLE SIZE =====================\n');
fprintf('File: %s   (%d x %d px', opt.dataFile, nx, ny);
if useUm, fprintf(', pixel pitch %.3f um', px); end
fprintf(')\n');
fprintf('Speckle contrast: %.3f\n\n', contrast);

fprintf('Definition: lag where the normalized intensity ACF falls to the level.\n');
fprintf('  a  = half width, the alpha of Hamarova et al. (2016)\n');
fprintf('  2a = full width of the correlation peak, i.e. speckle diameter\n\n');

lvlName = {'1/e', '1/2 (HWHM)'};
methods = {'2D ACF profile      ', 'acf2D'; ...
           'average 1D (eq. 4-5)', 'avg1D'; ...
           'long vector 1D (2b) ', 'long1D'};

for m = 1:size(methods, 1)
    f = methods{m,2};
    fprintf('%s\n', methods{m,1});
    for t = 1:numel(opt.thresholds)
        ax = res.(f).x(t) * k;  ay = res.(f).y(t) * k;
        fprintf(['   r_I = %-10s   a_x = %6.3f %s   a_y = %6.3f %s' ...
                 '   |   2a_x = %6.3f %s   2a_y = %6.3f %s\n'], ...
                lvlName{t}, ax, unit, ay, unit, 2*ax, unit, 2*ay, unit);
    end
    if strcmp(f, 'avg1D')
        fprintf(['   (spread over profiles at 1/e: %.3f %s in x, %.3f %s in y;' ...
                 ' %d/%d rows, %d/%d cols used)\n'], ...
                res.avg1D.xStd(1)*k, unit, res.avg1D.yStd(1)*k, unit, ...
                nUsedX(1), size(Idec_rows,1), nUsedY(1), size(Idec_cols,2));
    end
end

% Pedestal diagnostic: far from the peak the ACF of a flattened speckle
% field must sit at zero. A clearly positive tail means the illumination
% envelope is still in the data, so raise opt.detrendSigma; a clearly
% negative tail means the smoothing has eaten into the speckle, so lower it.
tail = mean([acfX(end-9:end), acfY(end-9:end)]);
res.pedestal = tail;
fprintf('\nPedestal check: mean r_I over lags %d-%d = %+.4f (want |r_I| < 0.01)\n', ...
        lag1(end-9), lag1(end), tail);

aMean = mean([res.long1D.x(1), res.long1D.y(1)]);
fprintf('\nRecommended value (long-vector 1D method, 1/e level, x-y average):\n');
fprintf('   mean speckle size  a  = %.3f %s\n', aMean*k, unit);
fprintf('   speckle diameter  2a  = %.3f %s\n', 2*aMean*k, unit);
fprintf('=============================================================\n');

%% ------------------------------------------------------------------ plots
if opt.savePlots
    figure('Name', 'Speckle pattern', 'Color', 'w');
    subplot(1,2,1);
    imagesc(Iraw); axis image; colormap gray; colorbar;
    title('Raw intensity [counts]'); xlabel('x [px]'); ylabel('y [px]');
    subplot(1,2,2);
    n = min([1256, ny, nx]);
    imagesc(I(1000:n, 1000:n)); axis image; colormap gray; colorbar;
    title(sprintf('Detail %d x %d px (flattened)', n, n));
    xlabel('x [px]'); ylabel('y [px]');

    figure('Name', '2D autocorrelation', 'Color', 'w');
    imagesc(lagx, lagy, acf2); axis image; colorbar;
    xlabel('\Deltax'' [px]'); ylabel('\Deltay'' [px]');
    title('Normalized autocorrelation r_I(\Deltax'', \Deltay'')');

    figure('Name', 'ACF profiles', 'Color', 'w'); hold on; grid on;
    plot(lagx, profX, 'b-',  'LineWidth', 1.2, 'DisplayName', '2D ACF, x''');
    plot(lagy, profY, 'r-',  'LineWidth', 1.2, 'DisplayName', '2D ACF, y''');
    plot(lag1, acfX,  'b--', 'LineWidth', 1.2, 'DisplayName', 'long vector, x''');
    plot(lag1, acfY,  'r--', 'LineWidth', 1.2, 'DisplayName', 'long vector, y''');
    yline(1/exp(1), 'k:', '1/e');
    yline(0.5,      'k:', '1/2');
    xlabel('lag [px]'); ylabel('r_I'); xlim([0 opt.maxLag]);
    title('Normalized 1D autocorrelation of intensity');
    legend('Location', 'northeast');
end

save('speckle_size_results.mat', 'res', 'opt', 'acf2', 'lagx', 'lagy', ...
     'acfX', 'acfY', 'lag1', 'contrast');
fprintf('Results written to speckle_size_results.mat\n');

