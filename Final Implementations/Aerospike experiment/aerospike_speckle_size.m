%% AEROSPIKE_SPECKLE_SIZE  Mean speckle size of the aerospike reference images.
%
% Applies the method of speckle_size.m to every DaVis intensity export in the
% folder "Aerospike Speckle Data" and collects the results in one table.
%
% Method (Hamarova, Horvath, Smid & Hrabovsky, Measurement 88 (2016) 271-277):
% the mean speckle size is the lag at which the normalized autocorrelation of
% intensity falls to 1/e. Three estimators are evaluated for each file:
%   (A) 2D autocorrelation, profiles through the origin along x and y,
%   (B) average 1D correlation, one ACF per row and column (eq. 4-5),
%   (C) 1D correlation of a long vector, all rows (columns) concatenated
%       (Fig. 2b). This is the estimator the paper recommends and the one
%       reported as the result.
% Sizes are given as the full width 2a, the quantity that compares with the
% Rayleigh limit 1.22 lambda f# (1+M).
%
% Illumination flattening. Each image is divided by a Gaussian-smoothed copy
% of itself before the correlation is computed, to remove the beam envelope.
% With opt.detrendSigma = 'auto' the smoothing width is chosen per file as the
% one that drives the tail of the autocorrelation (mean over lags 51 to 60 px)
% to zero, as in speckle_size_batch.m: a bisection towards the sign change of
% the tail when one exists, otherwise a golden section search on its
% magnitude. The CROSSING column records which case applied. A number fixes
% the same width for every file instead.
%
% Scale. The exports carry millimetre coordinates, so their sample spacing is
% the pixel pitch in the calibrated plane and sizes are also reported in
% micrometres at that pitch. Set opt.pixelPitch_um to override it, for example
% with the 6.5 um sensor pixel to obtain the size on the sensor.
%
% Requires speckle_estimators.m and read_intensity_map.m from the speckle size
% investigation folder, which is added to the path below.
%
% Output: aerospike_speckle_size.csv and aerospike_speckle_acf.png next to
% this script, plus a printed table.
%
% Bruno Hegedus, MSc thesis - aerospike post processing.

clear; clc;

%% ------------------------------------------------------------------ input
opt.folder        = fileparts(mfilename('fullpath'));
opt.dataFolder    = fullfile(opt.folder, 'Aerospike Speckle Data');
opt.pattern       = '*.csv';
% speckle_estimators.m and read_intensity_map.m are kept in the sibling
% folder of this collection, so both experiments share one implementation.
opt.libFolder     = fullfile(fileparts(opt.folder), 'Speckle Pattern investigation');

opt.useCache      = true;     % keep parsed images as .mat in a cache subfolder
opt.pixelPitch_um = NaN;      % NaN -> spacing of the file's mm coordinates
opt.roi           = [];       % [x0 y0 width height] in px, [] -> whole frame
opt.detrendSigma  = 'auto';   % 'auto' or a fixed width in px
opt.sigmaGrid     = [6 10 15 22 32 45 65 90 130 180 250];  % search grid [px]
opt.refineIter    = 4;        % refinement steps after the grid
opt.maxLag        = 60;       % [px] lag window
opt.tailLags      = 51:60;    % [px] lags of the tail diagnostic
opt.thresholds    = [1/exp(1), 0.5];   % 1/e -> size, 1/2 -> FWHM
opt.savePlot      = false;
opt.outFile       = fullfile(opt.folder, 'aerospike_speckle_size.csv');
opt.plotFile      = fullfile(opt.folder, 'aerospike_speckle_acf.png');

assert(isfile(fullfile(opt.libFolder, 'speckle_estimators.m')) && ...
       isfile(fullfile(opt.libFolder, 'read_intensity_map.m')), ...
       'speckle_estimators.m and read_intensity_map.m not found in %s', ...
       opt.libFolder);
addpath(opt.libFolder);
est = speckle_estimators();

files = dir(fullfile(opt.dataFolder, opt.pattern));
assert(~isempty(files), 'No files matching %s in %s', opt.pattern, opt.dataFolder);

cacheFolder = fullfile(opt.dataFolder, 'cache');
if opt.useCache && ~isfolder(cacheFolder), mkdir(cacheFolder); end

%% ------------------------------------------------------------ process files
n = numel(files);
R = table('Size', [n 16], ...
    'VariableTypes', [{'string'}, repmat({'double'}, 1, 15)], ...
    'VariableNames', {'file', 'nx', 'ny', 'pitch_um', 'sigma_px', 'crossing', ...
                      'contrast', 'tail', 'twoa_px', 'twoa_um', 'fwhm_px', ...
                      'twoa_x_px', 'twoa_y_px', 'twoa_2D_px', 'twoa_avg1D_px', ...
                      'samples'});
acfStore = cell(n, 1);

tAll = tic;
for k = 1:n
    tf = tic;
    name = files(k).name;
    [~, stem] = fileparts(name);
    fprintf('%-18s ', name);

    % --- load, from cache when available ---------------------------------
    cacheFile = fullfile(cacheFolder, [stem '.mat']);
    if opt.useCache && isfile(cacheFile)
        S = load(cacheFile, 'I', 'pitch_um');
        I = S.I;  filePitch = S.pitch_um;
    else
        [I, filePitch] = read_intensity_map(fullfile(opt.dataFolder, name));
        if opt.useCache
            pitch_um = filePitch; %#ok<NASGU>
            save(cacheFile, 'I', 'pitch_um', '-v7.3');
        end
    end
    pitch = opt.pixelPitch_um;
    if isnan(pitch), pitch = filePitch; end

    if ~isempty(opt.roi)
        r = opt.roi;
        I = I(r(2)+(1:r(4)), r(1)+(1:r(3)));
    end
    [ny, nx] = size(I);

    % --- flattening width ---------------------------------------------------
    if ischar(opt.detrendSigma) || isstring(opt.detrendSigma)
        [sigma, crossing] = local_select_sigma(est, I, opt);
    else
        sigma = opt.detrendSigma;  crossing = NaN;
    end
    J = local_flatten(est, I, sigma);

    % --- (C) long vector, the reported estimator -------------------------
    [pX, lag] = est.acf1(reshape(J.', [], 1), opt.maxLag, nx);
    [pY, ~  ] = est.acf1(J(:),                opt.maxLag, ny);
    wX = est.widths(lag, pX, opt.thresholds);
    wY = est.widths(lag, pY, opt.thresholds);

    % --- (A) 2D and (B) average 1D, for comparison -------------------------
    [acf2, lx, ly] = est.acf2(J, opt.maxLag);
    w2x = est.widths(lx, acf2(ly == 0, :),   opt.thresholds);
    w2y = est.widths(ly, acf2(:, lx == 0).', opt.thresholds);
    wBx = est.avg1d(J,   opt.maxLag, opt.thresholds);
    wBy = est.avg1d(J.', opt.maxLag, opt.thresholds);

    twoa = 2*mean([wX(1) wY(1)]);
    tail = mean([pX(opt.tailLags+1) pY(opt.tailLags+1)]);
    K    = std(J(:))/mean(J(:));

    R(k,:) = {string(name), nx, ny, pitch, sigma, crossing, K, tail, ...
              twoa, twoa*pitch, 2*mean([wX(2) wY(2)]), 2*wX(1), 2*wY(1), ...
              2*mean([w2x(1) w2y(1)]), 2*mean([wBx(1) wBy(1)]), ...
              local_samples(fullfile(opt.dataFolder, name))};
    acfStore{k} = struct('lag', lag, 'x', pX, 'y', pY);

    fprintf('sigma %3d  2a = %6.3f px = %7.2f um  K %.3f  tail %+.4f  (%.0f s)\n', ...
            sigma, twoa, twoa*pitch, K, tail, toc(tf));
end

%% ----------------------------------------------------------------- report
writetable(R, opt.outFile);

fprintf('\n======================= AEROSPIKE SPECKLE SIZE =======================\n');
fprintf('Long vector method, 1/e level, full width 2a, x and y averaged.\n');
fprintf('Pitch from the file coordinates unless set by hand.\n\n');
fprintf('%-18s %6s %5s %8s %9s %8s %8s %8s %8s %8s\n', 'file', 'sigma', ...
        'cross', 'K', '2a [px]', '2a [um]', 'FWHM px', '2a_x', '2a_y', 'pitch');
for k = 1:n
    fprintf('%-18s %6d %5s %8.3f %9.3f %8.2f %8.3f %8.3f %8.3f %8.2f\n', ...
            R.file(k), R.sigma_px(k), local_yesno(R.crossing(k)), R.contrast(k), ...
            R.twoa_px(k), R.twoa_um(k), R.fwhm_px(k), R.twoa_x_px(k), ...
            R.twoa_y_px(k), R.pitch_um(k));
end
fprintf('\nCross-check, 2a [px] by estimator:\n');
fprintf('%-18s %10s %10s %10s\n', 'file', 'long 1D', '2D', 'avg 1D');
for k = 1:n
    fprintf('%-18s %10.3f %10.3f %10.3f\n', R.file(k), R.twoa_px(k), ...
            R.twoa_2D_px(k), R.twoa_avg1D_px(k));
end
fprintf('\nTotal time %.1f min. Table written to %s\n', toc(tAll)/60, opt.outFile);

%% ------------------------------------------------------------------ plot
if opt.savePlot
    f = figure('Color', 'w', 'Name', 'Aerospike speckle ACF');
    hold on; grid on; box on;
    col = lines(n);
    for k = 1:n
        a = acfStore{k};
        plot(a.lag, (a.x + a.y)/2, '-', 'Color', col(k,:), 'LineWidth', 1.4, ...
             'DisplayName', sprintf('%s,  2a = %.2f px', R.file(k), R.twoa_px(k)));
    end
    yline(1/exp(1), 'k:', '1/e', 'HandleVisibility', 'off');
    xlim([0 20]);
    xlabel('lag  [px]');
    ylabel('r_I  [--]');
    legend('Location', 'northeast', 'Interpreter', 'none');
    exportgraphics(f, opt.plotFile, 'Resolution', 200);
    fprintf('Saved %s\n', opt.plotFile);
end

%% --------------------------------------------------------------- helpers

function J = local_flatten(est, I, sigma)
% Divide out the illumination envelope, estimated by Gaussian smoothing.
    e = est.gaussblur(I, sigma);  e(e <= 0) = eps('single');
    J = double(I) ./ double(e) * mean(e(:));
end

function t = local_tail(est, I, sigma, opt)
% Mean normalized autocorrelation over the tail lags, both directions.
    [ny, nx] = size(I);
    J  = local_flatten(est, I, sigma);
    p1 = est.acf1(reshape(J.', [], 1), opt.maxLag, nx);
    p2 = est.acf1(J(:),                opt.maxLag, ny);
    t  = mean([p1(opt.tailLags+1) p2(opt.tailLags+1)]);
end

function [sigma, crossing] = local_select_sigma(est, I, opt)
% Width that brings the autocorrelation tail to zero: bisection on its sign
% change when there is one, golden section on its magnitude otherwise.
    tailOf = @(s) local_tail(est, I, s, opt);
    g = opt.sigmaGrid;
    t = arrayfun(tailOf, g);
    idx = find(t(1:end-1) .* t(2:end) < 0, 1, 'last');
    crossing = ~isempty(idx);
    if crossing
        lo = g(idx);  hi = g(idx+1);  tlo = t(idx);
        for it = 1:opt.refineIter
            mid = round((lo + hi)/2);
            if mid <= lo || mid >= hi, break, end
            tm = tailOf(mid);
            if sign(tm) == sign(tlo), lo = mid; tlo = tm; else, hi = mid; end
        end
        sigma = round((lo + hi)/2);
    else
        [~, j] = min(abs(t));
        lo = g(max(j-1, 1));  hi = g(min(j+1, numel(g)));
        sigma = local_golden(tailOf, lo, hi, opt.refineIter + 2);
    end
end

function s = local_golden(fun, lo, hi, nIter)
% Golden section search for the integer width minimizing |fun|.
    phi = (sqrt(5) - 1)/2;
    a = lo;  b = hi;
    c = round(b - phi*(b - a));   d = round(a + phi*(b - a));
    fc = abs(fun(c));             fd = abs(fun(d));
    for it = 1:nIter
        if b - a <= 2, break, end
        if fc < fd
            b = d;  d = c;  fd = fc;
            c = round(b - phi*(b - a));  if c <= a, c = a + 1; end
            fc = abs(fun(c));
        else
            a = c;  c = d;  fc = fd;
            d = round(a + phi*(b - a));  if d >= b, d = b - 1; end
            fd = abs(fun(d));
        end
    end
    if fc < fd, s = c; else, s = d; end
end

function n = local_samples(f)
% Number of averaged frames, read from the fourth column of the first row.
    fid = fopen(f, 'r');  fgetl(fid);  l = fgetl(fid);  fclose(fid);
    v = sscanf(strrep(l, ';', ' '), '%f');
    if numel(v) >= 4, n = v(4); else, n = NaN; end
end

function s = local_yesno(c)
    if isnan(c), s = '--'; elseif c, s = 'yes'; else, s = 'no'; end
end
