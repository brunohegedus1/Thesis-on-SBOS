%% SPECKLE_SIZE_SIGMA_BAND  Uncertainty of the flattening width, per point.
%
% The operating width is chosen by driving the correlation tail to zero. The
% tail is itself a measured quantity with a noise floor, so every width whose
% tail is statistically indistinguishable from zero is an admissible choice.
% The width of that admissible band is the uncertainty on sigma.
%
% The noise floor follows from the number of independent speckles in the
% frame. At a shift beyond the correlation length the autocorrelation has a
% standard error of about
%
%     eps = 1/sqrt(N_s),      N_s = (nx/2a) * (ny/2a)
%
% The band is the contiguous run of widths around the operating point for
% which |tail| <= eps, with the edges interpolated between grid values. When
% no width meets that criterion, the band is taken around the best width
% instead, as |tail| <= min|tail| + eps, and the CRITERION column records it.
%
% With no reason to prefer any width inside the band, a rectangular
% distribution gives the standard uncertainty
%
%     u(sigma) = half width / sqrt(3)
%
% which is what speckle_size_budget.m consumes.
%
% Output: speckle_size_sigma_band.csv plus a printed table.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc;
est = speckle_estimators();

folder = fileparts(mfilename('fullpath'));
B = readtable(fullfile(folder, 'speckle_size_batch.csv'));

opt.points   = 1:31;
opt.maxLag   = 60;
opt.tailLags = 51:60;
opt.grid     = [2 3 4 6 8 11 15 20 27 36 48 64 85 110 145 190 250];
opt.outFile  = fullfile(folder, 'speckle_size_sigma_band.csv');

% Rows for points that are not processed in this run are carried over from
% the existing file, so a subset can be redone without discarding the rest.
hdr = ['point,file,sigma_op,eps,sigma_lo,sigma_hi,halfwidth_px,' ...
       'u_sigma_px,open_edge,criterion'];
kept = local_carry_over(opt.outFile, hdr, opt.points);

fid = fopen(opt.outFile, 'w');
fprintf(fid, '%s\n', hdr);
for i = 1:numel(kept), fprintf(fid, '%s\n', kept{i}); end
fclose(fid);

fprintf('%5s %-10s %8s %8s %16s %10s %10s %s\n', 'point', 'file', ...
        'sigma', 'eps', 'band [px]', 'half', 'u(sigma)', 'note');

tAll = tic;
for i = 1:height(B)
    p = B.point(i);
    if ~ismember(p, opt.points), continue, end
    f = fullfile(folder, B.file{i});
    if ~isfile(f), fprintf('%5d %-10s missing\n', p, B.file{i}); continue, end

    I = read_intensity_map(f);
    [ny, nx] = size(I);
    twoa = B.twoa_1e_px(i);
    Ns   = (nx/twoa) * (ny/twoa);
    epsN = 1/sqrt(Ns);

    g = opt.grid;
    t = arrayfun(@(s) local_tail(est, I, s, opt), g);

    % Admissible widths, as a contiguous run containing the operating point.
    lim  = epsN;
    crit = 'noise floor';
    ok   = abs(t) <= lim;
    if ~any(ok)
        lim  = min(abs(t)) + epsN;
        crit = 'best + floor';
        ok   = abs(t) <= lim;
    end

    [~, iop] = min(abs(g - B.sigma_px(i)));      % nearest grid node
    if ~ok(iop), [~, iop] = min(abs(t)); end     % else start from the best
    lo = iop; while lo > 1 && ok(lo-1), lo = lo - 1; end
    hi = iop; while hi < numel(g) && ok(hi+1), hi = hi + 1; end

    % Interpolate the edges onto the |tail| = lim contour
    openEdge = 0;
    if lo > 1
        sLo = local_edge(g(lo-1), abs(t(lo-1)), g(lo), abs(t(lo)), lim);
    else
        sLo = g(1); openEdge = openEdge + 1;
    end
    if hi < numel(g)
        sHi = local_edge(g(hi), abs(t(hi)), g(hi+1), abs(t(hi+1)), lim);
    else
        sHi = g(end); openEdge = openEdge + 2;
    end

    half = (sHi - sLo)/2;
    uSig = half/sqrt(3);

    note = crit;
    if openEdge > 0, note = [note ', open edge']; end
    fprintf('%5d %-10s %8d %8.4f %7.1f to %-6.1f %10.2f %10.2f  %s\n', ...
            p, B.file{i}, B.sigma_px(i), epsN, sLo, sHi, half, uSig, note);

    fid = fopen(opt.outFile, 'a');
    fprintf(fid, '%d,%s,%d,%.6f,%.3f,%.3f,%.3f,%.3f,%d,%s\n', ...
            p, B.file{i}, B.sigma_px(i), epsN, sLo, sHi, half, uSig, ...
            openEdge, crit);
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

function t = local_tail(est, I, sigma, opt)
    [ny, nx] = size(I);
    e = est.gaussblur(I, sigma);  e(e <= 0) = eps('single');
    J = double(I) ./ double(e) * mean(e(:));
    p1 = est.acf1(reshape(J.', [], 1), opt.maxLag, nx);
    p2 = est.acf1(J(:),                opt.maxLag, ny);
    t  = mean([p1(opt.tailLags+1) p2(opt.tailLags+1)]);
end

function s = local_edge(s1, a1, s2, a2, lim)
% Width at which |tail| crosses LIM, linearly interpolated between two nodes.
    if a2 == a1, s = s1; return, end
    x = (lim - a1)/(a2 - a1);
    x = min(max(x, 0), 1);
    s = s1 + x*(s2 - s1);
end
