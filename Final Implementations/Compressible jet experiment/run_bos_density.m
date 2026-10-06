%RUN_BOS_DENSITY  End-to-end: raw .txt  ->  density and density gradient.
%
%   Chains the corrected processing route:
%
%     read_bos_txt        raw text  ->  gridded U, V
%     bos_pipeline        steps 1-7: background removal, centreline, spline,
%                         straighten, crop, smooth, fold      ->  R.Dn_f
%     bos_optical_path    deflection -> projected optical path S
%     abel_invert_columns S -> (n_ref - n0), the LOCAL refractive index excess
%     Gladstone-Dale      (n_ref - n0) -> rho
%     finite difference   rho -> d(rho)/dr
%
%   NOTE ON THE ORDER OF THE LAST STEPS
%   The transverse integration in BOS_OPTICAL_PATH happens BEFORE the Abel
%   inversion. That is deliberate: BOS measures a deflection, which is the
%   gradient of the projected field, not the projected field itself. Feeding
%   the raw deflection to a plain Abel inversion violates the transform's
%   precondition (d(rho)/dy = (d(rho)/dr)(y/r) is not a function of r alone)
%   and on a synthetic Gaussian produced a sign-flipped profile near the axis
%   and 97% error at the peak. Integrating first gives S, which IS a genuine
%   plain-Abel projection, and recovers the density to 0.03%.
%
%   Because of that ordering, ABEL_INVERT_COLUMNS here returns the refractive
%   index excess DIRECTLY -- not a gradient. So the density follows from
%   Gladstone-Dale alone. Do NOT also call NFBOS_DENSITY_AXISYM: that function
%   belongs to the older route and would apply Eq. (4) to a quantity that is
%   not a gradient and integrate something already integrated.
%
%   UNITS
%   All lengths (x, y, displacements, W, l0, dr) share one unit -- mm for
%   these exports. K is in m^3/kg, which sets the density unit to kg/m^3.
%   The mixed units are fine here: K only converts a DIMENSIONLESS refractive
%   index excess into a density, and S/dr is dimensionless by construction.
%   The reported gradient is therefore kg/m^3 per mm; multiply by 1000 for
%   kg/m^3 per m.

clear; clc; %close all;

%% ======================= USER SETTINGS ==============================
% --- input data ---------------------------------------------------------
% DATA_SOURCE picks the folder DATAFILE is read from:
%   'cc'   the cross-correlation .csv exports in CC_RESULTS_DIR
%   'txt'  the original .txt exports in BOS_DATA_DIR
% Both are the same export format -- identical header, semicolon
% delimiter, grid and row count -- so READ_BOS_TXT reads either one and
% nothing downstream changes.
DATA_SOURCE    = 'cc';
% The exports sit in 'CC results' next to this script.
CC_RESULTS_DIR = fullfile(fileparts(mfilename('fullpath')), 'CC results');

DATAFILE = 'BOS_3cm_12x120001_4.csv';
%DATAFILE = 'BOS_3cm_5bar_12x12.txt';   % lives in BOS_DATA_DIR

% --- gas properties (AIR jet) ---------------------------------
K        = 2.26e-4;   % Gladstone-Dale constant, m^3/kg. Dry air, visible
                      % light (Schmidt et al. 2025, Sec. II.A). Use 1.58e-4
                      % only if the jet gas is argon.
rho_amb  = 1.2;       % ambient air density, kg/m^3
n0       = 1 + K*rho_amb;

% --- optical geometry ----------------------------------------------------
% l0 = distance from the jet axis to the BACKGROUND PATTERN, in mm.
% l0 sets the lever arm and scales the whole density result linearly.
l0 = 30;

% W = local width of the Schlieren object (jet diameter), mm.
%   'auto'  estimate per station from the extent of the folded signal
%   numeric scalar or [1 x Ns] vector to set it yourself
%
% Set to a fixed value rather than 'auto'. The automatic estimate saturated
% against the crop for about half the stations and raised a warning, but the
% saturation turns out not to matter: W enters only through the lever arm
%
%     l_eff = W/2 + l0

W_MODE = 4.4;
W_FALLBACK = 4.0;     % used where the automatic estimate fails
W_CAP_FRAC = 0.6;     % never let the auto estimate exceed this fraction of
                      % the crop diameter -- beyond that it is measuring the
                      % crop, not the jet

% --- processing options --------------------------------------------------
opts = struct();
opts.centerlineDetection = false;
opts.SplineFit = false;
opts.centerMethod = 'zerocross'; %opts.centerMethod = 'argmax'; %opts.centerMethod = 'centroid';
opts.nKnots       = 5;
opts.nMax         = 3;      % symmetric crop half-height, mm
opts.RadialSmoothing = false;
opts.smoothWin    = 5;      % radial smoothing window, points
opts.foldHalves   = true;
opts.abelLambda   = 'gcv';  % Tikhonov strength: 'gcv' uses GCV to pick it from the data;
                            % 'lcurve' or a fixed number (old default 1) also work

SIGN_FLIP = false;          % set true if the S diagnostic warns (see below)

% --- streamwise orientation ---------------------------------------------
% BOS_PIPELINE builds the arc-length coordinate s from the low-x end of the
% detected centreline, which is not necessarily the nozzle. Set this to put
% the nozzle at s = 0 with s increasing downstream.
%   'auto'  decide from the data: the nozzle end carries the strongest
%           density gradients, so whichever end has the larger mean |Dn_f|
%           is taken to be the nozzle
%   true    force a flip
%   false   leave the pipeline's orientation alone
FLIP_STREAMWISE = 'auto';

% Streamwise position of the nozzle in the pipeline's own s coordinate, used
% to shift the origin so the nozzle sits at s = 0. 'auto' finds the last
% station on the nozzle side still carrying signal above the quiescent floor.
% Override with a number if the automatic value looks wrong -- the printed
% diagnostic reports what it found.
NOZZLE_S = 38.37;

%% ======================= 1. READ ====================================
switch lower(DATA_SOURCE)
    case 'cc',  dataDir = CC_RESULTS_DIR;
    case 'txt', dataDir = bos_data_dir;
    otherwise
        error('run_bos_gradients_unfolded:badSource', ...
              'DATA_SOURCE must be ''cc'' or ''txt'', got ''%s''.', DATA_SOURCE);
end
dataPath = fullfile(dataDir, DATAFILE);
D = read_bos_txt(dataPath);

%% ======================= 2. PIPELINE (steps 1-7) ====================
R = bos_pipeline(D, opts);
dr = median(diff(R.n));
Ns = numel(R.s);

fprintf('\nmean jet tilt      : %.2f deg\n', R.info.meanTiltDeg);
fprintf('centreline residual: %.4f mm\n',   R.info.centerResidStd);

%% ======================= 2b. STREAMWISE ORIENTATION =================
% Applied here, before anything else consumes R, so every downstream array
% inherits the same orientation. This is a pure relabelling: each station
% holds a profile in r (or n), so reordering the stations cannot change any
% profile and introduces no sign issues.
%
% Two separate corrections are needed:
%   (a) DIRECTION - the pipeline builds s from the low-x end of the
%       centreline, which may put the nozzle at high s.
%   (b) ORIGIN - the measurement domain usually extends a little past the
%       nozzle into dead air. Flipping alone would put s = 0 in that dead
%       region rather than at the nozzle, so the origin is shifted too.
nEdge = max(5, round(0.2*Ns));
sigStart = mean(abs(R.Dn_f(:, 1:nEdge)),        'all', 'omitnan');
sigEnd   = mean(abs(R.Dn_f(:, end-nEdge+1:end)),'all', 'omitnan');

if ischar(FLIP_STREAMWISE) && strcmpi(FLIP_STREAMWISE,'auto')
    doFlip = sigEnd > sigStart;
    fprintf('nozzle side        : |Dn_f| start %.3e vs end %.3e -> nozzle at %s\n', ...
            sigStart, sigEnd, ternary_local(doFlip,'HIGH s','LOW s'));
else
    doFlip = logical(FLIP_STREAMWISE);
end

% --- locate the nozzle: last station (on the nozzle side) still carrying
%     signal well above the quiescent floor -------------------------------
pk = max(abs(R.Dn_f), [], 1);              % per-station peak
pk = movmean_local(pk, 31);
liveThr = 3 * median(pk(~isnan(pk)) , 'omitnan') * 0.5;   % half the median, x3 floor
live = pk > liveThr;
if ischar(NOZZLE_S) && strcmpi(NOZZLE_S,'auto')
    if doFlip
        idxN = find(live, 1, 'last');
    else
        idxN = find(live, 1, 'first');
    end
    if isempty(idxN)
        s_nozzle = ternary_local(doFlip, max(R.s), min(R.s));
        warning('run_bos_density:noNozzle', ...
            'Could not locate the nozzle automatically; using the domain edge.');
    else
        s_nozzle = R.s(idxN);
    end
    fprintf('nozzle position    : s = %.2f mm (auto)\n', s_nozzle);
else
    s_nozzle = NOZZLE_S;
    fprintf('nozzle position    : s = %.2f mm (user set)\n', s_nozzle);
end

if doFlip
    colFields = {'Ds','Dn','Ds_sm','Dn_sm','Ds_f','Dn_f'};
    for q = 1:numel(colFields)
        f = colFields{q};
        if isfield(R,f) && ~isempty(R.(f)), R.(f) = fliplr(R.(f)); end
    end
    R.theta = fliplr(R.theta);
    R.s     = s_nozzle - fliplr(R.s);   % nozzle -> 0, downstream positive
else
    R.s     = R.s - s_nozzle;           % offset only
end
fprintf('streamwise frame   : nozzle at s = 0, downstream positive, s = %.1f to %.1f mm\n', ...
        min(R.s), max(R.s));
fprintf('  (s < 0 is behind the nozzle: expect ambient density there)\n');

%% ======================= 3. LOCAL WIDTH W ===========================
if ischar(W_MODE) && strcmpi(W_MODE,'auto')
    % Radius at which the folded signal drops into the outer-edge noise,
    % doubled to give a diameter.
    %
    % The threshold is set from the noise floor rather than as a fraction of
    % the station peak. A pure fraction-of-peak rule (10%) let noisy
    % downstream stations run out to the crop boundary: 32.5% of stations
    % saturated at the full 11.97 mm crop, which is an artefact of the crop
    % rather than a measurement of the jet. W is also capped at W_CAP_FRAC of
    % the crop for the same reason.
    noise = median(abs(R.Dn_f(end-10:end, :)), 1, 'omitnan');  % outer-edge level
    noiseFloor = median(noise(~isnan(noise)));
    W = nan(1, Ns);
    for k = 1:Ns
        prof = abs(R.Dn_f(:,k));
        if all(isnan(prof)), continue; end
        pk = max(prof);
        if pk <= 0, continue; end
        thr = max(0.25*pk, 3*noiseFloor);      % whichever is the stricter test
        idx = find(prof > thr, 1, 'last');
        if ~isempty(idx), W(k) = 2*R.nHalf(idx); end
    end
    Wcap = W_CAP_FRAC * 2 * R.nHalf(end);
    satFrac = mean(W(~isnan(W)) >= Wcap);
    W(W > Wcap) = Wcap;
    W(isnan(W)) = W_FALLBACK;
    W = movmean_local(W, 51);          % smooth along the jet
    fprintf('W (auto)           : %.2f to %.2f mm (%.0f%% capped)\n', ...
            min(W), max(W), 100*satFrac);
    if satFrac > 0.25
        warning('run_bos_density:Wsaturated', ...
            ['%.0f%% of stations hit the W cap, so the automatic width is ' ...
             'being set by the crop rather than by the jet. Consider a ' ...
             'tighter opts.nMax, or set W_MODE to a fixed value.'], 100*satFrac);
    end
else
    W = W_MODE;
    if isscalar(W), fprintf('W (fixed)          : %.2f mm\n', W); end
end

%% ======================= 4. OPTICAL PATH S ==========================
oS = struct('dr', dr, 'signFlip', SIGN_FLIP);
[S, epsAng, l_eff] = bos_optical_path(R.Dn_f, R.nHalf, W, l0, n0, oS);

fprintf('lever arm l_eff    : %.2f to %.2f mm\n', min(l_eff), max(l_eff));
fprintf('peak |deflection|  : %.3e rad\n', max(abs(epsAng(:))));
fprintf('peak S             : %.3e mm\n',  max(S(:)));

%% ======================= 5. ABEL -> REFRACTIVE INDEX ================
% Returns (n_ref - n0) directly, because S is already the projected field.
[dn_ref, abelInfo] = abel_invert_columns(S, dr, opts.abelLambda);
fprintf('Abel lambda        : %s, lambdaRel = %.3g\n', abelInfo.mode, abelInfo.lambdaRel);

%% ======================= 6. DENSITY (Gladstone-Dale) ================
rho = rho_amb + dn_ref / K;                      % kg/m^3, on r >= 0

%% ======================= 7. RADIAL DENSITY GRADIENT ================
% Central differences in r, one-sided at the ends. drho_dr(1) is forced to
% zero: rho must be smooth through the axis of an axisymmetric flow, so its
% radial derivative vanishes there by symmetry.
drho_dr = nan(size(rho));
for k = 1:Ns
    v = rho(:,k);
    if any(isnan(v)), continue; end
    g = zeros(size(v));
    g(2:end-1) = (v(3:end) - v(1:end-2)) / (2*dr);
    g(1)   = 0;
    g(end) = (v(end) - v(end-1)) / dr;
    drho_dr(:,k) = g;
end

%% ======================= 7b. STREAMWISE DENSITY GRADIENT ================
% Central differences in s on the LOCAL density field, one-sided at the ends.
%
% This replaces the earlier route, which applied the 2D Eq. (4)/(5)
% conversion directly to R.Ds_sm. That algebra was correct on its own terms
% but returned a PATH-AVERAGED gradient -- Eq. (4) assumes d(rho)/ds is
% constant over the whole width W -- and it divided by W at every n,
% including the outer rows where the true path length through the jet is
% zero. Plotting it next to the Abel-inverted radial gradient therefore
% compared a projected quantity with a local one. Differencing rho, which
% is already local, keeps both gradients in the same frame and needs no
% extra assumption about W.
%
% Stations that failed the Abel inversion are NaN in rho, so a NaN
% propagates into its two neighbours here. That is deliberate: a difference
% across a missing station is not a measurement.
ds = median(diff(R.s));
drho_ds = nan(size(rho));
drho_ds(:, 2:end-1) = (rho(:, 3:end) - rho(:, 1:end-2)) / (2*ds);
drho_ds(:, 1)       = (rho(:, 2)     - rho(:, 1))       /    ds;
drho_ds(:, end)     = (rho(:, end)   - rho(:, end-1))   /    ds;

%% ======================= 8. MIRROR TO THE FULL PLANE ================
% rho is EVEN about the axis, so its streamwise derivative is EVEN too;
% the cross-stream gradient is ODD.
mirrorFull = @(A,par) [par*flipud(A(2:end,:)); A];
rho_full     = mirrorFull(rho,     +1);
drho_dn_full = mirrorFull(drho_dr, -1);
drho_ds_full = mirrorFull(drho_ds, +1);

good = ~all(isnan(rho),1);
fprintf('stations recovered : %d of %d\n', nnz(good), Ns);
fprintf('density range      : %.3f to %.3f kg/m^3\n', ...
        min(rho(:,good),[],'all'), max(rho(:,good),[],'all'));

%% ======================= 9. PLOTS ===================================
% ---- figure 1: input and centerline -----------------------------------
figure();
subplot(1,2,1);
imagesc(D.x, D.y, D.mag); axis xy image; colorbar; clim([0 0.02]);
title('Raw displacement magnitude [mm]'); ylabel('y [mm]');

subplot(1,2,2);
magb = hypot(D.U - R.bgU, D.V - R.bgV);
imagesc(D.x, D.y, magb); axis xy image; colorbar; clim([0 0.02]); hold on;
plot(R.centerX, R.centerY, 'w.', 'MarkerSize', 3);
xx = linspace(min(R.centerX), max(R.centerX), 400);
plot(xx, ppval(R.spline, xx), 'm-', 'LineWidth', 1.5);
% window the background offset was measured in; it must not touch the jet
bw = R.info.bgWindow;                % empty when background removal is off
if ~isempty(bw)
    rectangle('Position', [min(bw(:,1)) min(bw(:,2)) abs(diff(bw(:,1))) abs(diff(bw(:,2)))], ...
              'EdgeColor', 'r', 'LineWidth', 1.2);
end
title('Background removed (red: background window), with centerline and spline');
xlabel('x [mm]'); ylabel('y [mm]');

% ---- figure 2: straightened fields ------------------------------------
figure();
subplot(3,1,1);
imagesc(R.s, R.n, R.Ds_sm); axis xy; colorbar; clim([-0.02 0.02]);

title('D_s  (streamwise component, straightened)'); ylabel('n [mm]');

subplot(3,1,2);
imagesc(R.s, R.n, R.Dn_sm); axis xy; colorbar; clim([-0.02 0.02]);
title('D_n  (cross-stream component, straightened) - note the antisymmetry about n = 0');
xlabel('s [mm]'); ylabel('n [mm]');

R.Dmag_sm = hypot(R.Ds_sm, R.Dn_sm);
subplot(3,1,3);
imagesc(R.s, R.n, R.Dmag_sm); axis xy; colorbar; clim([0 0.02]);
title('D mag straightened');
xlabel('s [mm]'); ylabel('n [mm]');

% ---- figure 3: folded -------------------------------
R.Dmag_f = hypot(R.Ds_f, R.Dn_f);
% Mirror with the same rule as MIRRORFULL in step 8. Dropping row 1 (r = 0)
% from the reflected half is what keeps the axis from appearing twice; the
% earlier A(1:end-1,:) form duplicated it and shifted the whole lower half
% by one cell relative to R.n.
R.Dmag_f = [flipud(R.Dmag_f(2:end,:)); R.Dmag_f];

figure();
imagesc(R.s, R.n, R.Dmag_f); axis xy; colorbar; clim([0 0.02]);
title('D mag after being folded'); ylabel('r [mm]');

% --- Fig 4: the density gradient, full plane (the requested plot) --------
figure('Color','w','Position',[60 60 1150 420]);
imagesc(R.s, R.n, drho_dn_full); axis xy; colorbar; clim([-5 5])
try, colormap(bluewhitered_local()); catch, colormap(jet); end
title('Radial density gradient  \partial\rho/\partial n   [kg m^{-3} mm^{-1}]');
xlabel('s  [mm]'); ylabel('n  [mm]');

figure('Color','w','Position',[60 60 1150 420]);
imagesc(R.s, R.n, drho_ds_full); axis xy; colorbar; clim([-2 2])
try, colormap(bluewhitered_local()); catch, colormap(jet); end
title('Streamwise density gradient  \partial\rho/\partial s   [kg m^{-3} mm^{-1}]');
xlabel('s  [mm]'); ylabel('n  [mm]');

% --- Fig 5: density field for context ------------------------------------
figure('Color','w','Position',[60 60 1150 420]);
imagesc(R.s, R.n, rho_full); axis xy; colorbar; clim([1 3])
title('Density  \rho  [kg m^{-3}]');
xlabel('s  [mm]'); ylabel('n  [mm]');

% --- Fig 6: gradient magnitude, half plane -------------------------------
figure('Color','w','Position',[60 60 1150 420]);
imagesc(R.s, R.nHalf, abs(drho_dr)); axis xy; colorbar; clim([0 5])
title('|\partial\rho/\partial r|   [kg m^{-3} mm^{-1}]');
xlabel('s  [mm]'); ylabel('r  [mm]');

% --- Fig 7: sanity checks and sample profiles ----------------------------
figure('Color','w','Position',[60 60 1150 700]);

% pick four well-separated stations and label them by streamwise position
gi    = find(good);
fracs = [0.1 0.2 0.4 0.6 0.8];
stn   = gi(max(1, min(numel(gi), round(numel(gi)*fracs))));
lbl   = arrayfun(@(j) sprintf('s = %.1f mm', R.s(j)), stn, ...
                 'UniformOutput', false);

subplot(2,2,1);
plot(R.nHalf, S(:, stn), 'LineWidth', 1.1);
grid on; xlabel('r  [mm]'); ylabel('S  [mm]');
title('S');
legend(lbl, 'Location', 'northeast');

subplot(2,2,2);
plot(R.nHalf, rho(:, stn), 'LineWidth', 1.1); hold on;
xl = xlim; plot(xl, [rho_amb rho_amb], 'k--', 'LineWidth', 1); hold off;
grid on; xlabel('r  [mm]'); ylabel('\rho  [kg m^{-3}]');
title('Radial density profiles');
legend([lbl, {'ambient'}], 'Location', 'northeast');

subplot(2,2,3);
plot(R.nHalf, drho_dr(:, stn), 'LineWidth', 1.1); hold on;
xl = xlim; plot(xl, [0 0], 'k--', 'LineWidth', 1); hold off;
grid on; xlabel('r  [mm]'); ylabel('\partial\rho/\partial r  [kg m^{-3} mm^{-1}]');
title('Radial gradient profiles');
legend([lbl, {'zero'}], 'Location', 'northeast');

subplot(2,2,4);
plot(R.s(good), rho(1, good), 'LineWidth', 1.2); hold on;
xl = xlim; plot(xl, [rho_amb rho_amb], 'k--', 'LineWidth', 1);
yl = ylim;
for q = 1:numel(stn)
    plot([R.s(stn(q)) R.s(stn(q))], yl, ':', 'Color', [0.6 0.6 0.6]);
end
hold off; ylim(yl);
grid on; xlabel('s  [mm]'); ylabel('\rho on axis  [kg m^{-3}]');
title('Centreline density');
legend({'centreline \rho', 'ambient', 'sampled stations'}, ...
       'Location', 'northeast');

% --- Fig 8: regularisation parameter ------------------------------------
% Only drawn when lambda was chosen automatically. The GCV minimum is often
% shallow, so check that the chosen lambda sits in a clear dip and not on a
% long plateau.
if ~strcmp(abelInfo.mode, 'fixed')
    lr = abelInfo.lambdaGrid / abelInfo.lambdaScale;
    iC = abelInfo.iChosen;
    figure('Color','w','Position',[60 60 1150 420]);
    subplot(1,2,1);
    loglog(lr, abelInfo.gcv, 'LineWidth', 1.1); hold on;
    if isscalar(iC), loglog(lr(iC), abelInfo.gcv(iC), 'ro'); end
    hold off; grid on;
    xlabel('\lambda_{rel}'); ylabel('GCV');
    title('Generalised cross-validation');

    subplot(1,2,2);
    if isvector(abelInfo.resNorm)
        loglog(abelInfo.resNorm, abelInfo.semiNorm, 'LineWidth', 1.1); hold on;
        loglog(abelInfo.resNorm(iC), abelInfo.semiNorm(iC), 'ro'); hold off;
        title('L-curve, all stations');
    else
        semilogx(lr, sum(abelInfo.curvature, 2, 'omitnan'));
        title('L-curve curvature, summed over stations');
    end
    grid on; xlabel('||A f - S||'); ylabel('||L f||');
    sgtitle(sprintf('Abel regularisation: %s, \\lambda_{rel} = %.3g', ...
            abelInfo.mode, median(abelInfo.lambdaRel, 'omitnan')));
end


%% ======================= local helpers ==============================
function y = movmean_local(x, w)
    x = x(:).'; n = numel(x); y = x; h = floor(w/2);
    for i = 1:n
        a = max(1, i-h); b = min(n, i+h);
        y(i) = mean(x(a:b));
    end
end

function out = ternary_local(cond, a, b)
    if cond, out = a; else, out = b; end
end

function cmap = bluewhitered_local()
    m = 256; h = floor(m/2);
    up   = [linspace(0,1,h).', linspace(0,1,h).', ones(h,1)];
    down = [ones(m-h,1), linspace(1,0,m-h).', linspace(1,0,m-h).'];
    cmap = [up; down];
end
