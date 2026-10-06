%% SPECKLE_SIZE_SENSITIVITY  Response of the measured speckle size to the
%% flattening width.
%
% For every intensity export, the measured speckle diameter is evaluated at
% the operating width sigma taken from speckle_size_batch.csv and at one
% pixel either side of it. The central difference
%
%     d(2a)/d(sigma)  =  [ 2a(sigma+1) - 2a(sigma-1) ] / 2
%
% is the sensitivity coefficient of the flattening term in the uncertainty
% budget: multiplied by the uncertainty on sigma itself, it gives the
% contribution of the flattening width to the uncertainty on the size.
%
% A large coefficient means the reported size is governed by an analysis
% choice rather than by the recorded pattern, so this quantity also serves
% as a quality indicator independent of the correlation tail.
%
% Requires speckle_size_batch.csv, which supplies the operating sigma.
% Output: speckle_size_sensitivity.csv plus a printed table.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc;
est = speckle_estimators();

opt.folder     = fileparts(mfilename('fullpath'));
opt.batchFile  = fullfile(opt.folder, 'speckle_size_batch.csv');
opt.outFile    = fullfile(opt.folder, 'speckle_size_sensitivity.csv');
opt.thresholds = 1/exp(1);
opt.maxLag     = 60;
opt.step       = 1;            % [px] the "one unit" of the question
opt.points     = 1:31;         % which points to process, the rest carried over

% The file for each point, and its operating width, come from the batch
% results, so an override recorded there (point 3 on its re-acquisition) is
% picked up here without repeating it.
B = readtable(opt.batchFile);

hdr = ['point,file,sigma_px,twoa_lo,twoa_mid,twoa_hi,' ...
       'd2a_dsigma_px,pct_per_px'];
kept = local_carry_over(opt.outFile, hdr, opt.points);

fid = fopen(opt.outFile, 'w');
fprintf(fid, '%s\n', hdr);
for i = 1:numel(kept), fprintf(fid, '%s\n', kept{i}); end
fclose(fid);

fprintf('%-5s %-10s %6s %9s %9s %9s %12s %11s\n', 'point', 'file', 'sigma', ...
        '2a(s-1)', '2a(s)', '2a(s+1)', 'd2a/dsigma', '[%/px]');

tAll = tic;
for k = 1:height(B)
    p = B.point(k);
    if ~ismember(p, opt.points), continue, end
    name = B.file{k};
    f = fullfile(opt.folder, name);
    if ~isfile(f), fprintf('%-9s missing, skipped\n', name); continue, end

    [I, ~] = read_intensity_map(f);
    s = B.sigma_px(k);

    a_lo  = local_size(est, I, s - opt.step, opt);
    a_mid = local_size(est, I, s,            opt);
    a_hi  = local_size(est, I, s + opt.step, opt);

    d   = (a_hi - a_lo) / (2*opt.step);       % [px of size per px of sigma]
    pct = d / a_mid * 100;

    fprintf('%-5d %-10s %6d %9.4f %9.4f %9.4f %12.5f %11.3f\n', ...
            p, name, s, a_lo, a_mid, a_hi, d, pct);

    fid = fopen(opt.outFile, 'a');
    fprintf(fid, '%d,%s,%d,%.5f,%.5f,%.5f,%.6f,%.4f\n', ...
            p, name, s, a_lo, a_mid, a_hi, d, pct);
    fclose(fid);
end
sort_results_file(opt.outFile);
fprintf('\nTotal time %.1f min. Written to %s\n', toc(tAll)/60, opt.outFile);

%% --------------------------------------------------------------- helpers

function kept = local_carry_over(outFile, hdr, points)
% Rows of an existing results file whose point is not being reprocessed.
    kept = {};
    if ~isfile(outFile), return, end
    txt = strsplit(fileread(outFile), newline);
    txt = txt(~cellfun(@isempty, strtrim(txt)));
    if isempty(txt) || ~strcmp(strtrim(txt{1}), hdr)
        warning(['%s has a different layout and is not carried over. ' ...
                 'Process every point to rebuild it.'], outFile);
        return
    end
    for i = 2:numel(txt)
        p = sscanf(txt{i}, '%d', 1);
        if ~isempty(p) && ~ismember(p, points)
            kept{end+1} = strtrim(txt{i}); %#ok<AGROW>
        end
    end
end

function a2 = local_size(est, I, sigma, opt)
% Speckle diameter at the 1/e level for one flattening width, x-y averaged.
    [ny, nx] = size(I);
    e = est.gaussblur(I, sigma);  e(e <= 0) = eps('single');
    J = double(I) ./ double(e) * mean(e(:));
    [p1, lag] = est.acf1(reshape(J.', [], 1), opt.maxLag, nx);
    [p2, ~  ] = est.acf1(J(:),                opt.maxLag, ny);
    ax = est.widths(lag, p1, opt.thresholds);
    ay = est.widths(lag, p2, opt.thresholds);
    a2 = 2 * mean([ax ay]);
end
