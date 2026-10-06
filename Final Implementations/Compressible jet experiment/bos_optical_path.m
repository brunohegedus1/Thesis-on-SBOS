function [S, eps, l_eff] = bos_optical_path(Dn_f, nHalf, W, l0, n0, opts)
%BOS_OPTICAL_PATH  Projected optical path S from measured BOS displacements.
%
%   [S, EPS, L_EFF] = BOS_OPTICAL_PATH(DN_F, NHALF, W, L0, N0)
%   [...] = BOS_OPTICAL_PATH(..., OPTS)
%
%   This is the step that must come BETWEEN bos_pipeline and the Abel
%   inversion when the measurement is a DEFLECTION (as BOS always is).
%
%   WHY THIS EXISTS
%   A BOS camera measures how far a ray was bent, and a ray is bent only by
%   the refractive-index gradient TRANSVERSE to its own path. For an
%   axisymmetric field that introduces a projection factor y/r, so the
%   measured deflection obeys
%
%       eps(y) = 2y * integral_y^R  (dn_ref/dr) dr / sqrt(r^2 - y^2)
%
%   In order to apply Abel's inversion the identity  eps = d(S)/dy is needed, where S is the
%   projected optical path. S is a genuine plain-Abel projection of the
%   refractive-index excess, so once you
%   have S the existing chord-length matrix in ABEL_INVERT_COLUMNS is the
%   correct operator. Integrating first also avoids building the deflection
%   kernel itself, which is singular on axis (its first row is all zeros).
%
%   INPUTS
%     Dn_f  : [Nr x Ns] FOLDED, PROJECTED cross-stream displacement, i.e.
%             R.Dn_f from bos_pipeline
%             Row 1 = axis (r = 0). Units of length (mm for these exports).
%     nHalf : [Nr x 1] radial coordinate, R.nHalf. Same length units.
%     W     : local width of the Schlieren object. Scalar or [1 x Ns].
%     l0    : object-to-background distance (scalar), same length units.
%     n0    : ambient refractive index (dimensionless), e.g. 1.000271 for air.
%
%   OPTS (optional)
%     opts.dr       : radial spacing, default median(diff(nHalf)).
%     opts.signFlip : multiply the deflection by -1 before integrating,
%                     default false. See the sign check below.
%     opts.quiet    : suppress the diagnostic printout, default false.
%
%   OUTPUTS
%     S     : [Nr x Ns] projected optical path excess, ready to hand
%             straight to ABEL_INVERT_COLUMNS. Scaled so that the Abel
%             inverse of S is (n_ref - n0) directly.
%     eps   : [Nr x Ns] deflection angle in radians, the intermediate.
%     l_eff : [1 x Ns] lever arm used, (W + 2*l0)/Ccorr.
%
%   SIGN CONVENTION -- CHECK THIS
%   Whether a positive Dn_f means the ray bent toward or away from the axis
%   depends on your camera orientation and on how the straightening in
%   bos_pipeline defined the normal direction. Getting it backwards flips
%   the sign of S and yields a jet LESS dense than ambient. Since a
%   compressed jet core must be denser, this function warns when S comes out
%   negative on axis. If it does, set opts.signFlip = true.
%
%   NEXT STEP
%       S      = bos_optical_path(R.Dn_f, R.nHalf, W, l0, n0);
%       dn_ref = abel_invert_columns(S, median(diff(R.n)), 'gcv');
%       rho    = (n0 - 1)/K + dn_ref/K;

    if nargin < 6 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'dr')       || isempty(opts.dr),       opts.dr = median(diff(nHalf)); end
    if ~isfield(opts,'signFlip') || isempty(opts.signFlip), opts.signFlip = false; end
    if ~isfield(opts,'quiet')    || isempty(opts.quiet),    opts.quiet = false; end

    [Nr, Ns] = size(Dn_f);
    if numel(nHalf) ~= Nr
        error('bos_optical_path:sizeMismatch', ...
            'nHalf has %d entries but Dn_f has %d rows.', numel(nHalf), Nr);
    end
    if isscalar(W)
        W = W*ones(1,Ns);
    elseif numel(W) ~= Ns
        error('bos_optical_path:badW', ...
            'W must be scalar or have one entry per station (%d), got %d.', Ns, numel(W));
    end
    W = W(:).';

    % --- lever arm: apparent shift = l_eff * deflection angle -------------
    % Ccorr is the near-field correction of Eqs. (15)-(16); the remaining
    % (W + 2*l0) is the geometric arm from Eq. (4). Note the path-length W
    % is NOT included here -- the Abel step downstream handles the
    % integration through the object.
    Ccorr = 2;%./ (1 + W./l0);
    l_eff = (W + 2*l0) ./ Ccorr;

    % --- step 1: deflection angle -----------------------------------------
    eps = Dn_f ./ repmat(l_eff, Nr, 1);
    if opts.signFlip, eps = -eps; end

    % --- step 2: integrate transversely, inward from the outer edge -------
    % S(y) = -integral_y^R eps dy'. S vanishes outside the jet, so the
    % outermost row anchors the integration at zero.
    S = nan(Nr, Ns);
    dr = opts.dr;
    for k = 1:Ns
        e = eps(:,k);
        if any(isnan(e)), continue; end
        v = zeros(Nr,1);
        for i = Nr-1:-1:1
            v(i) = v(i+1) - 0.5*(e(i) + e(i+1))*dr;
        end
        S(:,k) = v;
    end

    % scale so that abel_invert_columns(S) returns (n_ref - n0) directly
    S = n0 * S;

    % --- sign diagnostic ---------------------------------------------------
    good = ~all(isnan(S),1);
    axisVals = S(1, good);
    if ~isempty(axisVals)
        fracNeg = mean(axisVals < 0);
        if ~opts.quiet
            fprintf('bos_optical_path: %d of %d stations integrated\n', nnz(good), Ns);
            fprintf('  S on axis: median %+.4e, %.0f%% of stations negative\n', ...
                median(axisVals), 100*fracNeg);
        end
        if fracNeg > 0.5
            warning('bos_optical_path:negativeS', ...
                ['S is negative on axis for %.0f%% of stations, implying a jet core ' ...
                 'LESS dense than ambient. For a compressed jet this usually means ' ...
                 'the displacement sign convention is reversed -- try opts.signFlip = true.'], ...
                 100*fracNeg);
        end
    end
end
