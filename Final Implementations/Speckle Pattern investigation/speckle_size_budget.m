%% SPECKLE_SIZE_BUDGET  Post-processing uncertainty, point by point.
%
% Combines the terms that the analysis itself contributes to the measured
% speckle diameter.
%
%   1. Flattening width. The sensitivity d(2a)/d(sigma) measured by
%      speckle_size_sensitivity.m, multiplied by the uncertainty on the
%      width itself. That uncertainty comes from the range of widths over
%      which the correlation tail stays inside its noise floor, measured per
%      point by speckle_size_sigma_band.m, with the fallback below used only
%      when no band is available.
%
%   2. Estimator choice. Spread between the 2D, average 1D and long vector
%      estimators, taken at the figure Hamarova et al. (2016) report.
%
% Output: speckle_size_budget.csv plus a printed table.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

clear; clc;

folder = fileparts(mfilename('fullpath'));
B = readtable(fullfile(folder, 'speckle_size_batch.csv'));
S = readtable(fullfile(folder, 'speckle_size_sensitivity.csv'));

% Measured widths of the admissible sigma band, when available. Each point
% then carries its own u(sigma) instead of the fallback below.
bandFile = fullfile(folder, 'speckle_size_sigma_band.csv');
if isfile(bandFile), W = readtable(bandFile); else, W = []; end

opt.uSigma_px      = 10;     % [px] fallback when no band was measured
opt.sigmaFloorFactor = 1;    % lowest admissible width, in speckle diameters
opt.uEstimator_pct = 1.5;    % [%]  spread between the three estimators
opt.coverage       = 2;      % k for the expanded uncertainty
opt.outFile        = fullfile(folder, 'speckle_size_budget.csv');

fid = fopen(opt.outFile, 'w');
fprintf(fid, ['point,file,twoa_px,u_flat_pct,u_estim_pct,' ...
              'u_total_pct,u_total_px,U_k2_px,u_sigma_px,' ...
              'u_sigma_source\n']);
fclose(fid);

fprintf('%5s %-10s %8s | %8s %8s | %9s %11s | %7s %s\n', ...
        'point', 'file', '2a [px]', 'flat', 'estim', ...
        'u total', 'U (k=2)', 'u(sig)', 'source');
fprintf('%s\n', repmat('-', 1, 92));

for i = 1:height(B)
    p    = B.point(i);
    twoa = B.twoa_1e_px(i);
    j    = find(S.point == p, 1);
    if isempty(j)
        warning('Point %d has no sensitivity entry, skipped.', p); continue
    end

    % The measured band comes from the tail criterion alone, which does not
    % exclude widths narrower than the speckle: over-flattening silences the
    % tail for the wrong reason, by removing the speckle correlation along
    % with the pedestal. The band is therefore clipped from below at
    % SIGMAFLOORFACTOR speckle diameters before the uncertainty is taken.
    uSig = opt.uSigma_px;  src = 'assumed';
    if ~isempty(W)
        w = find(W.point == p, 1);
        if ~isempty(w)
            sLo  = max(W.sigma_lo(w), opt.sigmaFloorFactor * twoa);
            sHi  = max(W.sigma_hi(w), sLo);
            uSig = (sHi - sLo)/2/sqrt(3);
            src  = 'measured';
        end
    end

    u1 = abs(S.d2a_dsigma_px(j)) * uSig / twoa * 100;               % [%]
    u2 = opt.uEstimator_pct;

    ut = sqrt(u1^2 + u2^2);

    fprintf(['%5d %-10s %8.3f | %7.2f%% %7.2f%% | ' ...
             '%8.2f%% %8.3f px | %7.1f %s\n'], ...
            p, B.file{i}, twoa, u1, u2, ut, ...
            opt.coverage*ut/100*twoa, uSig, src);

    fid = fopen(opt.outFile, 'a');
    fprintf(fid, '%d,%s,%.4f,%.3f,%.3f,%.3f,%.4f,%.4f,%.3f,%s\n', ...
            p, B.file{i}, twoa, u1, u2, ut, ut/100*twoa, ...
            opt.coverage*ut/100*twoa, uSig, src);
    fclose(fid);
end

fprintf('\nAssumptions: u(sigma) = %g px fallback, estimator %g%%, ', ...
        opt.uSigma_px, opt.uEstimator_pct);
fprintf('coverage factor k = %g.\n', opt.coverage);
fprintf('Written to %s\n', opt.outFile);
