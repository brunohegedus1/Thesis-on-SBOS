%% SPECKLE_SIZE_BATCH  Mean speckle size for every intensity export in a folder.
%
% Processes I1.csv ... I31.csv, which hold one intensity map each in the
% format exported by DaVis:
%
%   x [mm];y [mm];Average [];Number of data points [samples]
%
% Method and definitions follow speckle_size.m: the mean speckle size is the
% shift at which the normalized autocorrelation of intensity falls to 1/e,
% evaluated with the long-vector 1D method of Hamarova et al. (2016).
%
% Illumination flattening. Each image is divided by a Gaussian-smoothed copy
% of itself before the correlation is computed. The smoothing width sigma is
% chosen per file as the value that drives the tail of the autocorrelation
% (its mean over shifts of 51 to 60 px, far outside the correlation peak)
% closest to zero. A positive tail means beam envelope is still present, a
% negative one means the smoothing has attenuated the speckle itself, so the
% zero crossing separates the two errors. The search evaluates a coarse grid
% and then either bisects towards the crossing, when one exists, or runs a
% golden section search on the magnitude of the tail, when the tail keeps
% one sign over the whole range. The CROSSING column records which of the
% two applied: a file without a crossing does not satisfy the criterion, it
% only approaches it as closely as the range of widths allows, and its size
% should be treated as the less reliable of the two cases.
%
% Output: speckle_size_batch.csv in this folder, plus a printed table.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc;
est = speckle_estimators();

opt.folder     = fileparts(mfilename('fullpath'));

% The intensity maps, I4.csv to I31.csv plus the three corrected exports of
% points 1 to 3, sit next to this script.
opt.dataFolder = opt.folder;
opt.pattern    = 'I%d.csv';
opt.points     = 1:31;      % which points of the test matrix to process
opt.thresholds = [1/exp(1), 0.5];
opt.maxLag     = 60;
opt.tailLags   = 51:60;        % shifts used for the tail diagnostic [px]
opt.sigmaGrid  = [6 10 15 22 32 45 65 90 130 180 250];  % coarse search [px]
opt.refineIter = 4;            % bisection steps after the coarse grid
opt.outFile    = fullfile(opt.folder, 'speckle_size_batch.csv');

% Points whose data does not follow the default naming. Points 1 to 3 use the
% corrected export. The first three acquisitions were originally numbered one
% place out, which put the f/22.63 measurement below the f/11.31 one, an
% ordering the aperture cannot produce. Under the corrected numbering the
% sweep is monotonic and matches points 4 to 6, which share the same three
% f-numbers. The three files hold the same acquisitions as the originals,
% reassigned, not repeated measurements.
opt.override = { 1, 'I1_2.txt'
                 2, 'I2_2.txt'
                 3, 'I3_2.txt' };

% Rows for points that are not processed in this run are carried over from
% the existing file, so a single point can be redone without repeating the
% whole set.
hdr = ['point,file,sigma_px,a_1e_px,twoa_1e_px,fwhm_px,twoa_1e_mm,' ...
       'pitch_um,contrast,tail,crossing,twoa_x_px,twoa_y_px'];
kept = local_carry_over(opt.outFile, hdr, opt.points);

fid = fopen(opt.outFile, 'w');
fprintf(fid, '%s\n', hdr);
for i = 1:numel(kept), fprintf(fid, '%s\n', kept{i}); end
fclose(fid);

tAll = tic;
for k = opt.points
    name = sprintf(opt.pattern, k);
    for i = 1:size(opt.override, 1)
        if opt.override{i,1} == k, name = opt.override{i,2}; end
    end
    f = fullfile(opt.dataFolder, name);
    if ~isfile(f)
        fprintf('%-10s missing, skipped\n', name);
        continue
    end
    tf = tic;

    [I, pitch_um] = read_intensity_map(f);
    pitch_mm = pitch_um / 1000;
    [ny, nx] = size(I);

    % --- sigma that puts the autocorrelation tail closest to zero ---------
    tailOf = @(s) local_tail(est, I, s, opt);

    g = opt.sigmaGrid;
    t = arrayfun(tailOf, g);

    % Two cases. If the tail changes sign somewhere on the grid, the two
    % errors balance at the crossing and a bisection finds it exactly.
    % Otherwise the tail keeps the same sign over the whole range, the two
    % errors never balance, and the best available choice is the width that
    % minimizes the magnitude of the tail. A bisection is meaningless there,
    % so a golden section search on |tail| is used instead. The distinction
    % is recorded, because a case without a crossing does not satisfy the
    % criterion, it only comes as close to it as the range allows.
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

    % --- final evaluation at the selected sigma ---------------------------
    e = est.gaussblur(I, sigma);  e(e <= 0) = eps('single');
    J = double(I) ./ double(e) * mean(e(:));

    [px_, lag] = est.acf1(reshape(J.', [], 1), opt.maxLag, nx);
    [py_, ~  ] = est.acf1(J(:),                opt.maxLag, ny);
    ax = est.widths(lag, px_, opt.thresholds);
    ay = est.widths(lag, py_, opt.thresholds);

    a    = mean([ax(1) ay(1)]);
    fwhm = 2*mean([ax(2) ay(2)]);
    K    = std(J(:))/mean(J(:));
    tail = mean([px_(opt.tailLags+1) py_(opt.tailLags+1)]);

    fid = fopen(opt.outFile, 'a');
    fprintf(fid, '%d,%s,%d,%.4f,%.4f,%.4f,%.5f,%.4f,%.4f,%.5f,%d,%.4f,%.4f\n', ...
            k, name, sigma, a, 2*a, fwhm, 2*a*pitch_mm, pitch_mm*1000, ...
            K, tail, crossing, 2*ax(1), 2*ay(1));
    fclose(fid);

    fprintf(['pt%-3d %-10s %4d x %4d px  pitch %.3f um  sigma %3d  ' ...
             '2a = %6.3f px = %.4f mm  K %.3f  tail %+.4f  cross %d  (%.0f s)\n'], ...
            k, name, nx, ny, pitch_mm*1000, sigma, 2*a, 2*a*pitch_mm, ...
            K, tail, crossing, toc(tf));
end

sort_results_file(opt.outFile);
fprintf('\nTotal time %.1f min. Table written to %s\n', ...
        toc(tAll)/60, opt.outFile);

%% ---------------------------------------------------------------- helpers

function kept = local_carry_over(outFile, hdr, points)
% Rows of an existing results file whose point is not being reprocessed.
% Returns them in point order, so redoing one point leaves the rest intact.
    kept = {};
    if ~isfile(outFile), return, end
    txt = strsplit(fileread(outFile), newline);
    txt = txt(~cellfun(@isempty, strtrim(txt)));
    if isempty(txt) || ~strcmp(strtrim(txt{1}), hdr)
        warning(['%s has a different layout and is not carried over. ' ...
                 'Process every point to rebuild it.'], outFile);
        return
    end
    num = [];
    for i = 2:numel(txt)
        p = sscanf(txt{i}, '%d', 1);
        if ~isempty(p) && ~ismember(p, points)
            kept{end+1} = strtrim(txt{i}); %#ok<AGROW>
            num(end+1)  = p;               %#ok<AGROW>
        end
    end
    [~, ord] = sort(num);
    kept = kept(ord);
end

function t = local_tail(est, I, sigma, opt)
% Mean normalized autocorrelation over the tail shifts, both directions,
% for a given flattening width. This is the quantity driven to zero.
    [ny, nx] = size(I);
    e = est.gaussblur(I, sigma);  e(e <= 0) = eps('single');
    J = double(I) ./ double(e) * mean(e(:));
    p1 = est.acf1(reshape(J.', [], 1), opt.maxLag, nx);
    p2 = est.acf1(J(:),                opt.maxLag, ny);
    t  = mean([p1(opt.tailLags+1) p2(opt.tailLags+1)]);
end

function s = local_golden(fun, lo, hi, nIter)
% Golden section search for the integer width that minimizes |fun|, used
% when the tail keeps one sign over the whole range and therefore has no
% crossing to bisect towards.
    phi = (sqrt(5) - 1)/2;
    a = lo;  b = hi;
    c = round(b - phi*(b - a));   d = round(a + phi*(b - a));
    fc = abs(fun(c));             fd = abs(fun(d));
    for it = 1:nIter
        if b - a <= 2, break, end
        if fc < fd
            b = d;  d = c;  fd = fc;
            c = round(b - phi*(b - a));
            if c <= a, c = a + 1; end
            fc = abs(fun(c));
        else
            a = c;  c = d;  fc = fd;
            d = round(a + phi*(b - a));
            if d >= b, d = b - 1; end
            fd = abs(fun(d));
        end
    end
    if fc < fd, s = c; else, s = d; end
end

% The reader lives in read_intensity_map.m, shared with speckle_size.m.
