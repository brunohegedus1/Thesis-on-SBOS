%RUN_BOS_GRADIENTS_UNFOLDED  Density gradients on the UNFOLDED (s,n) grid.
%
%   Plots d(rho)/dn and d(rho)/ds across the full straightened field, both
%   halves, BEFORE the fold and before the Abel inversion. Run it; three
%   figures and a printout: the whole field, a square ZOOM_SPAN-mm window
%   drawn at a true 1:1 aspect, and d(rho)/ds along the centreline.
%
%   WHAT THIS IS AND IS NOT
%   BOS measures a deflection, and a deflection is the LINE-OF-SIGHT
%   INTEGRAL of the refractive-index gradient. Only the Abel inversion makes
%   a quantity local, and the Abel inversion needs a folded half-profile. So
%   there is no local d(rho)/dn to be had at this stage, and this script does
%   not pretend otherwise.
%
%   GRADIENT_MODE selects which of the two thesis relations is applied. The
%   thesis equation numbers are set once, in EQ_PLANAR and EQ_VECTOR, and
%   the printout and every title read them from there. In the thesis
%   version of 2026-09-15 they are Eqs. 2.22 and 2.21.
%
%     'planar'  (default)  eps_x = (G*W/n0) * drho/dx     [kg/m^3 per mm]
%               inverted to drho/dx = n0*eps_x/(G*W). The thesis derives
%               this from the vector relation for the case where the
%               transverse gradients do NOT vary along the line of sight,
%               i.e. a 2D PLANAR object of thickness W.
%
%     'vector'             eps_vec = (G/n0) * int grad(rho) dz  [kg/m^3]
%               inverted to int grad(rho) dz = (n0/G)*eps_vec. No W, no
%               assumption about how the gradient varies along the ray.
%
%   The old names 'eq2.18' and 'eq2.17' are still accepted.
%
%   READ THE ASSUMPTION IN THE PLANAR RELATION. This jet is AXISYMMETRIC, so
%   a ray crosses a gradient that rotates and changes sign along its own
%   path, and the premise of the planar relation is violated. What it
%   returns here is therefore the EQUIVALENT PLANAR gradient: the uniform
%   drho/dx that would have bent the ray by as much as the real field did.
%   It is the right quantity to plot for a schlieren-like view of the
%   structure and it carries the right units, but it understates the true
%   peak and smears the shock cells over the chord. Recovering the local
%   value is what the Abel-route equations (2.23 and 2.24 in the current
%   thesis) exist for, and that is what RUN_BOS_DENSITY does.
%
%   So: compare the magnitudes here against RUN_BOS_DENSITY's drho_dr and
%   drho_ds expecting this script to read LOW. Compare the structure and
%   the cell spacing expecting agreement.
%
%   WHY LOOK AT IT UNFOLDED
%   Because this is where the fold can be checked before it is committed to.
%   For an axisymmetric jet on a correctly located axis:
%
%       d(rho)/dn  is ODD  about n = 0   (it must vanish on the axis)
%       d(rho)/ds  is EVEN about n = 0
%
%   Both follow from rho(r) being a function of |n| alone. The printed parity
%   diagnostic measures how far the data departs from that. It is the most
%   direct test of the centreline fit there is: a mislocated or tilted axis
%   breaks the parity immediately, and one RMS number captures it. Once the
%   halves are folded, that information is averaged away and cannot be
%   recovered.
%
%   THE CONVERSION, AND WHERE IT COMES FROM
%   Both relations live in mainmatter/02_TheoreticalBackground.tex:
%
%     vector  (Eq. 2.21)  \label{equation:epsilon-nabla rho}
%                   eps_vec = (G/n0) int grad(rho) dz
%     planar  (Eq. 2.22)  \label{equation: Displacement and Density Gradient}
%                   eps_x = (G*W/n0) drho/dx
%
%   solved for the gradient in each case. The planar relation is the vector
%   one with grad(rho) pulled out of the integral, which the thesis states
%   requires the transverse gradients to be independent of the line of
%   sight.
%
%   NOTATION. The thesis writes the Gladstone-Dale coefficient G; this code
%   and RUN_BOS_DENSITY write K. Same quantity, 2.26e-4 m^3/kg for dry air.
%   The thesis line-of-sight element is z; here s is the STREAMWISE
%   coordinate of the straightened grid, so the planar relation's x maps to n or to s
%   depending on which displacement component is used.
%
%   The deflection itself comes from the measured displacement exactly as in
%   BOS_OPTICAL_PATH, so the two scripts cannot drift apart:
%
%       l_eff = (W + 2*l0) / Ccorr           lever arm, Ccorr = 2
%       eps   = D / l_eff                    deflection angle [rad]
%
%   giving, per component,
%
%       int (drho/dn) dz = n0 * eps_n / G,   eps_n = R.Dn_sm / l_eff
%       int (drho/ds) dz = n0 * eps_s / G,   eps_s = R.Ds_sm / l_eff
%
%   The same l_eff serves both: the lever arm is a property of the optical
%   layout, not of the direction the ray was bent in.
%
%   The vector relation gives a vector, so its magnitude is reported too:
%
%       | int grad(rho) dz | = (n0/G) * |eps| = (n0/G) * |D| / l_eff
%
%   Note this is |int grad(rho) dz| and NOT int |grad(rho)| dz. The two
%   differ wherever the gradient changes direction along the ray, which it
%   does across every shock in the cell. The magnitude map is the honest
%   schlieren-like image of the field; it is not a path-integrated |grad|.
%
%   CAVEAT CARRIED OVER FROM RUN_BOS_DENSITY
%   Ccorr = 2 is the far-field limit. l0 = 30 mm against a ~4 mm object is
%   not obviously far field, and l0 scales this result linearly. Treat the
%   MAGNITUDES as provisional; the STRUCTURE and the parity check do not
%   depend on either constant.
%
%   Settings below are kept in the same order and with the same names as
%   RUN_BOS_DENSITY so the two can be held in step by eye.

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

DATAFILE = 'BOS_3cm_12x120001_5.csv';
% Previous .txt input:   'BOS_3cm_8bar_12x12_M0.txt'  (DATA_SOURCE = 'txt')

% Dataset nomenclature. Titles name the dataset by letter, not by file name.
% A file that is not in this list (a .txt input, for example) gets no
% letter. The file names carry no pressure, so check which run each letter
% is before you quote a result against a tap setting.
DATASET_FILES = {'BOS_3cm_12x120001.csv',   'BOS_3cm_12x120001_1.csv', ...
                 'BOS_3cm_12x120001_2.csv', 'BOS_3cm_12x120001_3.csv', ...
                 'BOS_3cm_12x120001_4.csv', 'BOS_3cm_12x120001_5.csv'};
DATASET_TAGS  = {'a', 'b', 'c', 'd', 'e', 'f'};

GRADIENT_MODE = 'planar';   % 'planar'  planar relation, /(G*W)  [kg/m^3 per mm]
                            % 'vector'  vector integral form     [kg/m^3]

% Thesis equation numbers, used in the printout and in every title. Update
% these two lines if chapter 2 is renumbered; nothing else hardcodes them.
EQ_PLANAR = 'Eq. 2.22';     % \label{equation: Displacement and Density Gradient}
EQ_VECTOR = 'Eq. 2.21';     % \label{equation:epsilon-nabla rho}

% --- gas properties (AIR jet) -------------------------------------------
K        = 2.26e-4;   % Gladstone-Dale constant, m^3/kg (dry air, visible)
rho_amb  = 1.2;       % ambient air density, kg/m^3
n0       = 1 + K*rho_amb;

% --- optical geometry ---------------------------------------------------
l0 = 30;              % jet axis to background pattern, mm
Ccorr = 2;            % near-field correction; 2 is the far-field limit,
                      % matching BOS_OPTICAL_PATH

% W = local width of the Schlieren object (jet diameter), mm. Scalar or a
% [1 x Ns] vector. The 'auto' estimate of RUN_BOS_DENSITY is deliberately
% NOT carried over: W enters only through l_eff = (W + 2*l0)/Ccorr, so the
% whole plausible range of W moves the answer by ~6%, well under the l0
% uncertainty above. A fixed value keeps this script short and honest.
W = 4.0;

% --- export of figure 1, which is Figure 6.9 of the thesis ---------------
% The thesis includes this figure with \includegraphics[width=\textwidth],
% so it is exported exactly \textwidth wide and nothing is rescaled on the
% page. Text set at the document's \footnotesize therefore prints at
% \footnotesize. TUDELFT-REPORT loads book at 10 pt, where \footnotesize is
% 8 pt, and report.log gives \textwidth = 448.1309 pt = 6.20 in (TeX points,
% 1/72.27 in). The file name follows the dataset tag, so dataset f writes
% over 'Dataset f density gradients.jpg'.
EXPORT_FIG1   = true;
EXPORT_DIR    = 'C:\Users\bruno\Documents\GitHub\Thesis-on-SBOS\figures';
FONT_PT       = 8;                  % \footnotesize of the 10 pt class
FIG_WIDTH_IN  = 448.1309/72.27;     % \textwidth
FIG_HEIGHT_IN = 4.40;               % keeps the three stacked panels legible
EXPORT_DPI    = 600;

% --- processing options (must match RUN_BOS_DENSITY) --------------------
opts = struct();
opts.centerlineDetection = false;
opts.SplineFit    = false;
opts.centerMethod = 'zerocross';
opts.nKnots       = 5;
opts.nMax         = 5;      % symmetric crop half-height, mm
opts.RadialSmoothing = false;
opts.smoothWin    = 5;
opts.foldHalves   = true;   % kept TRUE: the folded field is used only to
                            % find the nozzle end, so that the s axis here
                            % matches RUN_BOS_DENSITY exactly. Nothing
                            % plotted below comes from it.

SIGN_FLIP = false;          % same convention as RUN_BOS_DENSITY

FLIP_STREAMWISE = 'auto';
NOZZLE_S = 38.37;


% --- square-window figure ------------------------------------------------
% Figure 2 draws a window ZOOM_SPAN mm on BOTH axes with a true 1:1 aspect,
% so the shock-cell geometry is undistorted. n is centred on the axis, so
% the window is n = +-ZOOM_SPAN/2 and s = ZOOM_S0 to ZOOM_S0 + ZOOM_SPAN.
ZOOM_SPAN = 6;        % mm, same extent on both axes
ZOOM_S0   = 0;        % mm, left edge of the window (0 = the nozzle)

% --- centreline figure ----------------------------------------------------
% Figure 3 plots the gradient along n = 0. A single pixel row is noisy, so
% it is also averaged over |n| <= CL_BAND. That average is safe for
% d(rho)/ds: it is EVEN about the axis and varies slowly across it. Set
% CL_BAND = 0 to plot the axis row only.
CL_BAND = 0.25;       % mm, half-width of the band averaged about the axis
% Figure 3 also marks the cell features: compression peaks, expansion
% troughs, density minima and maxima, and the end of the first cell. They
% are DETECTED on a copy of the band mean smoothed over CL_SMOOTH mm and
% drawn on the unsmoothed curve. Over 0.2 to 0.8 mm the first-cell length
% moved by 0.07 mm at most on the six CC datasets.
CL_SMOOTH = 0.4;      % mm, smoothing window used only for feature detection
%% ======================= 1. READ ====================================
switch lower(DATA_SOURCE)
    case 'cc',  dataDir = CC_RESULTS_DIR;
    case 'txt', dataDir = bos_data_dir;
    otherwise
        error('run_bos_gradients_unfolded:badSource', ...
              'DATA_SOURCE must be ''cc'' or ''txt'', got ''%s''.', DATA_SOURCE);
end
dataPath = fullfile(dataDir, DATAFILE);
if ~isfile(dataPath)
    % Name what IS there, so a typo in DATAFILE is a one-look fix.
    avail = [dir(fullfile(dataDir, '*.csv')); dir(fullfile(dataDir, '*.txt'))];
    error('run_bos_gradients_unfolded:noFile', ...
          '%s not found in\n  %s\nAvailable:\n  %s', DATAFILE, dataDir, ...
          strjoin({avail.name}, '\n  '));
end
iTag = find(strcmpi(DATASET_FILES, DATAFILE), 1);
if isempty(iTag), dsTag = ''; else, dsTag = sprintf('Dataset %s', DATASET_TAGS{iTag}); end
fprintf('input              : %s  (%s)%s\n', DATAFILE, upper(DATA_SOURCE), ...
        tern(isempty(dsTag), '', ['  ->  ' dsTag]));
D = read_bos_txt(dataPath);

%% ======================= 2. PIPELINE (steps 1-7) ====================
R  = bos_pipeline(D, opts);
dn = median(diff(R.n));
ds = median(diff(R.s));
Ns = numel(R.s);

fprintf('\nmean jet tilt      : %.2f deg\n', R.info.meanTiltDeg);
fprintf('centreline residual: %.4f mm\n',   R.info.centerResidStd);

%% ======================= 2b. STREAMWISE ORIENTATION =================
% Copied from RUN_BOS_DENSITY so both scripts put the nozzle at s = 0 with
% s increasing downstream. A pure relabelling of the station axis: each
% column holds a profile in n, so reordering columns cannot alter a profile.
nEdge    = max(5, round(0.2*Ns));
sigStart = mean(abs(R.Dn_f(:, 1:nEdge)),         'all', 'omitnan');
sigEnd   = mean(abs(R.Dn_f(:, end-nEdge+1:end)), 'all', 'omitnan');

if ischar(FLIP_STREAMWISE) && strcmpi(FLIP_STREAMWISE,'auto')
    doFlip = sigEnd > sigStart;
    fprintf('nozzle side        : |Dn_f| start %.3e vs end %.3e -> nozzle at %s\n', ...
            sigStart, sigEnd, tern(doFlip,'HIGH s','LOW s'));
else
    doFlip = logical(FLIP_STREAMWISE);
end

if ischar(NOZZLE_S) && strcmpi(NOZZLE_S,'auto')
    pk = movmean_l(max(abs(R.Dn_f), [], 1), 31);
    live = pk > 1.5*median(pk(~isnan(pk)), 'omitnan');
    idxN = find(live, 1, tern(doFlip,'last','first'));
    if isempty(idxN)
        s_nozzle = tern(doFlip, max(R.s), min(R.s));
        warning('run_bos_gradients_unfolded:noNozzle', ...
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
    for f = {'Ds','Dn','Ds_sm','Dn_sm','Ds_f','Dn_f'}
        if isfield(R,f{1}) && ~isempty(R.(f{1})), R.(f{1}) = fliplr(R.(f{1})); end
    end
    R.theta = fliplr(R.theta);
    R.s     = s_nozzle - fliplr(R.s);
    % Reversing the s axis reverses the direction a POSITIVE streamwise
    % displacement points in, so the streamwise components change sign.
    % Reordering the columns alone kept the old sign: d(rho)/ds came out
    % positive through the lip expansion, where density must fall. The
    % cross-stream components do not change, since n is not reversed.
    for f = {'Ds','Ds_sm','Ds_f'}
        if isfield(R,f{1}) && ~isempty(R.(f{1})), R.(f{1}) = -R.(f{1}); end
    end
else
    R.s     = R.s - s_nozzle;
end
fprintf('streamwise frame   : nozzle at s = 0, s = %.1f to %.1f mm\n', ...
        min(R.s), max(R.s));

%% ======================= 3. UNFOLDED GRADIENTS ======================
% The two straightened, smoothed components, on the FULL n grid. No fold,
% no Abel inversion. R.Dn_sm is the cross-stream shift, R.Ds_sm the
% streamwise one, both in mm.
Dn = R.Dn_sm;
Ds = R.Ds_sm;
if SIGN_FLIP, Dn = -Dn;  Ds = -Ds; end

if isscalar(W), Wv = W*ones(1,Ns); else, Wv = W(:).'; end
if numel(Wv) ~= Ns
    error('run_bos_gradients_unfolded:badW', ...
          'W must be scalar or have one entry per station (%d).', Ns);
end
l_eff = (Wv + 2*l0) / Ccorr;                      % [1 x Ns], mm

% The vector relation (EQ_VECTOR) inverted, per component. Units kg/m^3 (i.e. kg/m^3*mm of
% column density excess, per mm of n or s).
Lrep      = repmat(l_eff, numel(R.n), 1);
eps_n     = Dn ./ Lrep;                  % deflection angle [rad]
eps_s     = Ds ./ Lrep;
grad_n_pr = n0 * eps_n / K;
grad_s_pr = n0 * eps_s / K;

switch lower(GRADIENT_MODE)
    case {'planar','eq2.18','eq218','pathaveraged'}
        % Planar relation:  eps_x = (G*W/n0) drho/dx  ->  drho/dx = n0 eps/(G*W)
        Wrep   = repmat(Wv, numel(R.n), 1);
        grad_n = grad_n_pr ./ Wrep;
        grad_s = grad_s_pr ./ Wrep;
        uLab   = 'kg m^{-3} mm^{-1}';
        eqTx   = EQ_PLANAR;
        % \epsilon, not \varepsilon: MATLAB's TeX interpreter has no
        % \varepsilon, and one unknown symbol blanks the whole title.
        modeTx    = [eqTx ', equivalent planar gradient (n_0\epsilon/GW)'];
        modePlain = [eqTx ', equivalent planar gradient n0*eps/(G*W)'];
    case {'vector','eq2.17','eq217','projected'}
        % Vector relation:  int grad(rho) dz = (n0/G) eps_vec
        grad_n = grad_n_pr;  grad_s = grad_s_pr;
        uLab   = 'kg m^{-3}';
        eqTx   = EQ_VECTOR;
        modeTx    = [eqTx ', \int\nabla\rho dz (n_0\epsilon/G)'];
        modePlain = [eqTx ', int grad(rho) dz = n0*eps/G'];
    otherwise
        error('run_bos_gradients_unfolded:badMode', ...
              'GRADIENT_MODE must be ''planar'' or ''vector''.');
end

% Magnitude of the vector relation. See the header: this is
% |int grad(rho) dz|, not int |grad(rho)| dz.
grad_mag = hypot(grad_n, grad_s);

fprintf('\nlever arm l_eff    : %.2f to %.2f mm\n', min(l_eff), max(l_eff));
fprintf('peak |deflection|  : %.3e rad\n', max(hypot(eps_n(:), eps_s(:))));
fprintf('gradient mode      : %s [%s]\n', modePlain, uLab);
fprintf('\ngradient from thesis %s,  G = K = %.3e m^3/kg\n', eqTx, K);
fprintf('  peak |drho/dn|   : %7.3f   RMS %7.4f\n', ...
        max(abs(grad_n(:))), rms_l(grad_n));
fprintf('  peak |drho/ds|   : %7.3f   RMS %7.4f\n', ...
        max(abs(grad_s(:))), rms_l(grad_s));
fprintf('  peak |grad rho|  : %7.3f   RMS %7.4f   [%s]\n', ...
        max(grad_mag(:)), rms_l(grad_mag), uLab);

%% ======================= 4. PARITY DIAGNOSTIC =======================
% Pair each n = +a with n = -a. R.n contains an exact zero and the crop is
% symmetric, so the pairing is exact; the code still takes the shorter side
% in case a future crop is not.
[~, i0] = min(abs(R.n));
nUp = numel(R.n) - i0;  nLo = i0 - 1;  nH = min(nUp, nLo);
iU  = i0 + (1:nH);                       % n > 0, increasing
iL  = i0 - (1:nH);                       % n < 0, mirrored to match iU

% Split each field into its even and odd parts about n = 0. One part is the
% parity the physics demands and the fold keeps; the other is contamination
% the fold discards. Report the ratio of the two -- NOT the raw sum over the
% field RMS, which double-counts and makes a 45% contamination read as 82%.
gnU = grad_n(iU,:);  gnL = grad_n(iL,:);
gnKeep = 0.5*(gnU - gnL);                % ODD  part: what d(rho)/dn should be
gnDrop = 0.5*(gnU + gnL);                % EVEN part: contamination

gsU = grad_s(iU,:);  gsL = grad_s(iL,:);
gsKeep = 0.5*(gsU + gsL);                % EVEN part: what d(rho)/ds should be
gsDrop = 0.5*(gsU - gsL);                % ODD  part: contamination

contamN = rms_l(gnDrop)/rms_l(gnKeep);
contamS = rms_l(gsDrop)/rms_l(gsKeep);
keptN   = rms_l(gnKeep)/rms_l([gnU; gnL]);
keptS   = rms_l(gsKeep)/rms_l([gsU; gsL]);

fprintf('\nparity about n = 0.  For an axisymmetric jet on a correct axis the\n');
fprintf('wrong-parity part is zero; step 7 discards whatever is left of it.\n');
fprintf('  %-26s %14s %14s\n', '', 'contamination', 'RMS kept');
fprintf('  %-26s %13.1f%% %13.1f%%\n', 'drho/dn (should be ODD)', ...
        100*contamN, 100*keptN);
fprintf('  %-26s %13.1f%% %13.1f%%\n', 'drho/ds (should be EVEN)', ...
        100*contamS, 100*keptS);
% Repeat over the jet core only. The obvious objection to the numbers above
% is that they are set by noise in the quiescent outer field, where the
% right-parity part is near zero and the ratio blows up. Restricting to
% |n| < CORE_N and to stations still carrying signal tests that directly.
% If the core-only figure is as large as the whole-field one, or larger, the
% contamination is in the jet, not in the background.
CORE_N = 2.5;                                   % mm
rowCore = R.n(iU) < CORE_N;
pkS     = max(abs(R.Dn_f), [], 1);
colCore = pkS > 0.25*max(pkS);
cc = @(U,L,sg) rms_l(0.5*(U(rowCore,colCore) + sg*L(rowCore,colCore)));
contamNc = cc(gnU,gnL,+1)/cc(gnU,gnL,-1);
contamSc = cc(gsU,gsL,-1)/cc(gsU,gsL,+1);
fprintf('  %-26s %13.1f%% %13s   (|n| < %.1f mm, %d of %d stations)\n', ...
        'same, jet core only', 100*contamNc, '-', CORE_N, nnz(colCore), Ns);
fprintf('  %-26s %13.1f%% %13s\n', '', 100*contamSc, '');
fprintf(['  contamination = RMS(wrong parity) / RMS(right parity)\n' ...
         '  RMS kept      = fraction of the field RMS the fold retains\n']);
fprintf('  Centreline residual for reference: %.4f mm\n', R.info.centerResidStd);
fprintf(['  These are a direct test of the axis: a mislocated or tilted\n' ...
         '  centreline turns odd signal into even signal and shows up here\n' ...
         '  before the fold hides it.\n']);

%% ======================= 5. FIGURES =================================
cmap = bwr_l();
cN = clim_l(grad_n);  cS = clim_l(grad_s);

% ---- figure 1: the two gradient maps, unfolded ------------------------
figure('Color','w','Position',[60 60 1250 880]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(R.s, R.n, grad_n); axis xy; colormap(gca, cmap); clim([-cN cN]);
cb = colorbar; cb.Label.String = uLab;
hold on; yline(0, 'k--', 'LineWidth', 0.8); xline(0, 'k:', 'LineWidth', 0.8);
title(sprintf('\\partial\\rho/\\partial n'));
ylabel('n [mm]');

nexttile;
imagesc(R.s, R.n, grad_s); axis xy; colormap(gca, cmap); clim([-cS cS]);
cb = colorbar; cb.Label.String = uLab;
hold on; yline(0, 'k--', 'LineWidth', 0.8); xline(0, 'k:', 'LineWidth', 0.8);
title('\partial\rho/\partial s');
ylabel('n [mm]');

% Magnitude of the vector relation, and the closest thing here to a
% conventional schlieren image. Sequential map, since it cannot be negative.
nexttile;
imagesc(R.s, R.n, grad_mag); axis xy; colormap(gca, parula);
clim([0 clim_l(grad_mag)]);
cb = colorbar; cb.Label.String = uLab;
hold on; yline(0, 'w--', 'LineWidth', 0.8); xline(0, 'w:', 'LineWidth', 0.8);
title(sprintf('|\\nabla\\rho|, vector magnitude'));
xlabel('s [mm]'); ylabel('n [mm]');

% No overall title: the thesis caption carries that text, and a title inside
% the image would repeat it.

% Size and fonts for the page, then export. Every text object is set to
% FONT_PT, titles included, so nothing is left at a MATLAB default that
% would print larger than \footnotesize.
fig1 = gcf;
set(findall(fig1, '-property', 'FontSize'), 'FontSize', FONT_PT);
set(fig1, 'Units', 'inches', 'PaperPositionMode', 'auto');
fig1.Position(3:4) = [FIG_WIDTH_IN FIG_HEIGHT_IN];

if EXPORT_FIG1
    if ~isfolder(EXPORT_DIR)
        error('run_bos_gradients_unfolded:noExportDir', ...
              'EXPORT_DIR does not exist:\n  %s', EXPORT_DIR);
    end
    fig1Name = fullfile(EXPORT_DIR, ...
        sprintf('%s density gradients.jpg', tern(isempty(dsTag), 'Dataset', dsTag)));
    [wIn, hIn] = export_at_width(fig1, fig1Name, FIG_WIDTH_IN, EXPORT_DPI);
    fprintf('\nfigure 1 exported   : %s\n', fig1Name);
    fprintf('                      %.3f x %.3f in at %d dpi, %g pt text\n', ...
            wIn, hIn, EXPORT_DPI, FONT_PT);
    fprintf('                      \\textwidth is %.3f in, so the page scale is %.3f\n', ...
            FIG_WIDTH_IN, FIG_WIDTH_IN/wIn);
end


% ---- figure 2: square window, 1:1 aspect ------------------------------
% Figure 1 stretches 38 mm of s across the same width as 10 mm of n, about
% 4:1, so every cell is drawn wider than it is. This window is ZOOM_SPAN mm
% on BOTH axes and uses AXIS IMAGE, so a millimetre reads the same in either
% direction and the cell geometry is true: the barrel shock leaves the lip
% at an angle you can measure off the page.
sWin = [ZOOM_S0, ZOOM_S0 + ZOOM_SPAN];
nWin = [-ZOOM_SPAN/2, ZOOM_SPAN/2];
si = R.s >= sWin(1) & R.s <= sWin(2);
ni = R.n >= nWin(1) & R.n <= nWin(2);
if nnz(si) < 2 || nnz(ni) < 2
    warning('run_bos_gradients_unfolded:emptyZoom', ...
        ['The %.1f mm window at s = %.1f mm holds no data (s runs %.1f to ' ...
         '%.1f mm). Adjust ZOOM_S0 or ZOOM_SPAN.'], ...
         ZOOM_SPAN, ZOOM_S0, min(R.s), max(R.s));
else
    sZ = R.s(si);  nZ = R.n(ni);
    gnZ = grad_n(ni,si);  gsZ = grad_s(ni,si);  gmZ = grad_mag(ni,si);
    % Colour limits from the WINDOW, not the whole field. Reusing the
    % full-field limits would flatten this panel, since the near-lip
    % gradients are the strongest in the jet.
    czN = clim_l(gnZ);  czS = clim_l(gsZ);  czM = clim_l(gmZ);

    figure('Color','w','Position',[90 90 1400 560]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile;
    imagesc(sZ, nZ, gnZ); axis xy image; colormap(gca, cmap); clim([-czN czN]);
    cb = colorbar('southoutside'); cb.Label.String = uLab;
    hold on; yline(0,'k--','LineWidth',0.8);
    title('\partial\rho/\partial n');
    xlabel('s [mm]'); ylabel('n [mm]');

    nexttile;
    imagesc(sZ, nZ, gsZ); axis xy image; colormap(gca, cmap); clim([-czS czS]);
    cb = colorbar('southoutside'); cb.Label.String = uLab;
    hold on; yline(0,'k--','LineWidth',0.8);
    title('\partial\rho/\partial s');
    xlabel('s [mm]');

    nexttile;
    imagesc(sZ, nZ, gmZ); axis xy image; colormap(gca, parula); clim([0 czM]);
    cb = colorbar('southoutside'); cb.Label.String = uLab;
    hold on; yline(0,'w--','LineWidth',0.8);
    title('|\nabla\rho|');
    xlabel('s [mm]');

    sgtitle(tagTitle(dsTag, sprintf(['%s, %.1f x %.1f mm at true 1:1 aspect ' ...
                     '(s = %.1f to %.1f mm)'], ...
            eqTx, ZOOM_SPAN, ZOOM_SPAN, sWin(1), sWin(2))), 'FontWeight','bold');

    fprintf('\nsquare window: s = %.1f to %.1f mm, n = %+.1f to %+.1f mm\n', ...
            sWin(1), sWin(2), nWin(1), nWin(2));
    fprintf('  %d x %d pixels, peak |grad rho| in window %.3f [%s]\n', ...
            nnz(ni), nnz(si), max(gmZ(:)), uLab);
end

% ---- figure 3: centreline density gradient ----------------------------
% On the axis only d(rho)/ds carries signal. d(rho)/dn is ODD about n = 0,
% so it must vanish there. It is not plotted; its RMS on the axis row is
% printed below as a check on the axis position.
bandRows = abs(R.n) <= CL_BAND + 1e-9;
gsAxis   = grad_s(i0,:);                         % the axis row alone
gsBand   = mean(grad_s(bandRows,:), 1, 'omitnan');
gnAxis   = grad_n(i0,:);

% Cell features. Along the axis each cell reads:
%   expansion trough -> density MINIMUM (zero crossing, - to +) ->
%   compression peak (the shock) -> density MAXIMUM (zero crossing, + to -)
% The density maximum after the first compression peak ends the first cell,
% measured from the nozzle at s = 0.
F = cl_features(R.s, gsBand, CL_SMOOTH);

figure('Color','w','Position',[120 120 1250 460]);
hold on; grid on; box on;
hA = plot(R.s, gsAxis, '-', 'Color',[0.65 0.75 0.95], 'LineWidth',0.8);
hB = plot(R.s, gsBand, '-', 'Color',[0.00 0.25 0.70], 'LineWidth',1.6);
yline(0, 'k-', 'HandleVisibility','off');
xline(0, 'k:', 'nozzle', 'LabelVerticalAlignment','bottom', ...
      'HandleVisibility','off');

colC   = [0.85 0.00 0.00];      % compression peak   (red)
colE   = [0.00 0.00 0.00];      % expansion trough   (black)
colMin = [0.00 0.60 0.00];      % density minimum    (green)
colMax = [0.85 0.45 0.00];      % density maximum    (orange)
hF = gobjects(0);  lgF = {};
% Peaks and troughs sit on the plotted (unsmoothed) band mean; crossings
% sit on the zero line, where they are defined.
iC = F.idx(F.lobes(F.lobes(:,1) > 0, 4));
iE = F.idx(F.lobes(F.lobes(:,1) < 0, 4));
sMin = F.cross(F.cross(:,2) > 0, 1);
sMax = F.cross(F.cross(:,2) < 0, 1);
if ~isempty(iC)
    hF(end+1) = plot(R.s(iC), gsBand(iC), '^', 'LineStyle','none', 'MarkerSize',7, ...
                     'MarkerFaceColor',colC, 'MarkerEdgeColor',colC);
    lgF{end+1} = 'compression peak';
end
if ~isempty(iE)
    hF(end+1) = plot(R.s(iE), gsBand(iE), 'v', 'LineStyle','none', 'MarkerSize',7, ...
                     'MarkerFaceColor',colE, 'MarkerEdgeColor',colE);
    lgF{end+1} = 'expansion trough';
end
if ~isempty(sMin)
    hF(end+1) = plot(sMin, zeros(size(sMin)), 'o', 'LineStyle','none', 'MarkerSize',6, ...
                     'MarkerFaceColor',colMin, 'MarkerEdgeColor',colMin);
    lgF{end+1} = 'density minimum';
end
if ~isempty(sMax)
    hF(end+1) = plot(sMax, zeros(size(sMax)), 's', 'LineStyle','none', 'MarkerSize',6, ...
                     'MarkerFaceColor',colMax, 'MarkerEdgeColor',colMax);
    lgF{end+1} = 'density maximum';
end
if isfinite(F.Lmax)
    hF(end+1) = xline(F.Lmax, '--', 'Color',colMax, 'LineWidth',1.2);
    lgF{end+1} = sprintf('end of first cell, L_{max} = %.2f mm', F.Lmax);
end

xlim([min(R.s) max(R.s)]);
xlabel('s [mm]');
ylabel(sprintf('\\partial\\rho/\\partial s on the axis  [%s]', uLab));
legend([hA hB hF], ...
       [{'\partial\rho/\partial s, axis row (n = 0)', ...
         sprintf('\\partial\\rho/\\partial s, mean over |n| \\leq %.2f mm (%d rows)', ...
                 CL_BAND, nnz(bandRows))}, lgF], ...
       'Location','northeast', 'Box','off', 'NumColumns', 2);
%title(tagTitle(dsTag, sprintf('centreline density gradient, %s', eqTx)), ...
%      'FontWeight','bold');

[gMax, iMax] = max(gsBand);  [gMin, iMin] = min(gsBand);
fprintf('\ncentreline, mean over |n| <= %.2f mm (%d rows):\n', CL_BAND, nnz(bandRows));
fprintf('  max drho/ds  %+.3f at s = %.2f mm   (density rising, compression)\n', ...
        gMax, R.s(iMax));
fprintf('  min drho/ds  %+.3f at s = %.2f mm   (density falling, expansion)\n', ...
        gMin, R.s(iMin));
fprintf('  RMS of drho/dn on the axis row: %.4f  (zero for a perfect axis)\n', ...
        rms_l(gnAxis));
fprintf('\nfirst shock cell (features detected over a %.2f mm smoothing window):\n', CL_SMOOTH);
fprintf('  L_max  %5.2f mm   nozzle to the first density maximum after the shock\n', F.Lmax);
fprintf('  L_min  %5.2f mm   between the first two density minima\n', F.Lmin);
fprintf('  L_avg  %5.2f mm   (L_max + L_min)/2\n', (F.Lmax + F.Lmin)/2);
fprintf('  L_C    %5.2f mm   between the first two compression peaks\n', F.LC);

%% ======================= local helpers ==============================
function out = tern(c,a,b), if c, out = a; else, out = b; end, end

function [wIn, hIn] = export_at_width(fig, fname, targetWin, dpi)
%Export FIG so the written file is exactly TARGETWIN inches wide.
%   EXPORTGRAPHICS crops to the tight content box, so the file comes out
%   narrower than the figure and LaTeX would then scale it up at
%   \includegraphics[width=\textwidth], enlarging the text with it. The
%   figure width is therefore adjusted until the cropped file matches.
%   Font sizes are absolute points and do not change with the figure size.
%   The crop removes a nearly constant margin, so the width correction is
%   additive and converges in one or two passes. The probe is written at the
%   export resolution, because a coarser one quantises the crop.
    tmp = [tempname '.png'];
    cleaner = onCleanup(@() delete_if_present(tmp)); %#ok<NASGU>
    for it = 1:4
        exportgraphics(fig, tmp, 'Resolution', dpi);
        info = imfinfo(tmp);
        wNow = info.Width / dpi;
        if abs(wNow - targetWin) < 0.005, break; end
        fig.Position(3) = fig.Position(3) + (targetWin - wNow);
    end
    exportgraphics(fig, fname, 'Resolution', dpi);
    info = imfinfo(fname);
    wIn = info.Width / dpi;  hIn = info.Height / dpi;
end

function delete_if_present(f)
    if isfile(f), delete(f); end
end

function t = tagTitle(tag, t)
%Prefix a title with the dataset tag, or capitalise it when there is none.
    if isempty(tag)
        t(1) = upper(t(1));
    else
        t = [tag ': ' t];
    end
end

function y = movmean_l(x, w)
    x = x(:).'; n = numel(x); y = x; h = floor(w/2);
    for i = 1:n
        y(i) = mean(x(max(1,i-h):min(n,i+h)), 'omitnan');
    end
end

function v = rms_l(A)
    a = A(:); a = a(isfinite(a));
    if isempty(a), v = NaN; else, v = sqrt(mean(a.^2)); end
end

function F = cl_features(s, g, winMM)
%Shock-cell features of the centreline d(rho)/ds.
%   Smooths G over WINMM mm and marks LOBES where the smoothed signal stays
%   beyond a noise threshold: +1 = compression, -1 = expansion. Each lobe
%   keeps its extremum. One zero crossing is placed between each pair of
%   neighbouring lobes: +1 = density minimum (- to +), -1 = density maximum.
%
%   F.lobes  [sign iStart iEnd iExtremum], indices into F.s
%   F.idx    maps an index into F.s back to an index into S
%   F.cross  [s, type]
%   F.Lmax   nozzle (s = 0) to the first density maximum after the first
%            compression peak
%   F.Lmin   between the density minima before the first and second
%            compression peaks
%   F.LC     between the first two compression peaks
    s = s(:).';  g = g(:).';
    w   = max(3, round(winMM/median(diff(s))));
    gsm = movmean(g, w, 'omitnan');
    idx = find(s >= 0 & isfinite(gsm));
    s   = s(idx);  gs = gsm(idx);
    F = struct('s',s, 'gs',gs, 'idx',idx, 'thr',NaN, 'lobes',zeros(0,4), ...
               'cross',zeros(0,2), 'C1',NaN, 'Lmax',NaN, 'Lmin',NaN, 'LC',NaN);
    if numel(s) < 10, return; end

    % Threshold: 3 sigma of the far field, and never below 8% of the
    % near-field peak, so a quiet far field cannot promote noise to a lobe.
    far  = s >= 0.8*max(s);
    near = s <= 15;  if ~any(near), near = true(size(s)); end
    F.thr = max(3*std(gs(far)), 0.08*max(abs(gs(near))));

    lab = zeros(size(gs));  lab(gs > F.thr) = 1;  lab(gs < -F.thr) = -1;
    lobes = zeros(0,4);  i = 1;  n = numel(gs);
    while i <= n
        if lab(i) == 0, i = i + 1; continue; end
        j = i;
        while j < n && lab(j+1) == lab(i), j = j + 1; end
        if s(j) - s(i) >= 0.2                 % shorter than 0.2 mm is noise
            seg = i:j;
            if lab(i) > 0, [~,m] = max(gs(seg)); else, [~,m] = min(gs(seg)); end
            lobes(end+1,:) = [lab(i) i j seg(m)]; %#ok<AGROW>
        end
        i = j + 1;
    end
    if isempty(lobes), return; end

    % Neighbouring lobes of the same sign are one lobe split by noise.
    L = lobes(1,:);
    for q = 2:size(lobes,1)
        if lobes(q,1) == L(end,1)
            if abs(gs(lobes(q,4))) > abs(gs(L(end,4))), L(end,4) = lobes(q,4); end
            L(end,3) = lobes(q,3);
        else
            L(end+1,:) = lobes(q,:); %#ok<AGROW>
        end
    end
    F.lobes = L;

    cross = zeros(0,2);
    for q = 1:size(L,1)-1
        seg = L(q,4):L(q+1,4);
        sc  = find(diff(sign(gs(seg))) ~= 0);
        if isempty(sc), continue; end
        loc = zeros(size(sc));
        for t = 1:numel(sc)
            a = seg(sc(t));  b = a + 1;
            loc(t) = s(a) - gs(a)*(s(b)-s(a))/(gs(b)-gs(a));
        end
        cross(end+1,:) = [median(loc), L(q+1,1)]; %#ok<AGROW>
    end
    F.cross = cross;

    Cs = s(L(L(:,1) > 0, 4));
    if isempty(Cs), return; end
    F.C1 = Cs(1);
    mx = cross(cross(:,2) < 0 & cross(:,1) > Cs(1), 1);
    if ~isempty(mx), F.Lmax = mx(1); end
    if numel(Cs) >= 2
        F.LC = Cs(2) - Cs(1);
        mins = cross(cross(:,2) > 0, 1);
        m1 = mins(mins < Cs(1));
        m2 = mins(mins > Cs(1) & mins < Cs(2));
        if ~isempty(m1) && ~isempty(m2), F.Lmin = m2(end) - m1(end); end
    end
end

function c = clim_l(A)
%Robust symmetric colour limit: the 99th percentile of |A|, so a few
%outlying pixels cannot wash the map out.
    a = abs(A(:)); a = a(isfinite(a));
    if isempty(a), c = 1; else, c = max(quantile(a, 0.99), eps); end
end

function cmap = bwr_l()
    m = 256; h = floor(m/2);
    up   = [linspace(0,1,h).', linspace(0,1,h).', ones(h,1)];
    down = [ones(m-h,1), linspace(1,0,m-h).', linspace(1,0,m-h).'];
    cmap = [up; down];
end
