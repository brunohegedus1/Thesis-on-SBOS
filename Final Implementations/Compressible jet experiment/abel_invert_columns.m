function [f, info] = abel_invert_columns(F, dr, lambda, opts)
%ABEL_INVERT_COLUMNS  Regularised Abel inversion of each column of F.
%
%   f = ABEL_INVERT_COLUMNS(F, DR)
%   f = ABEL_INVERT_COLUMNS(F, DR, LAMBDA)
%   f = ABEL_INVERT_COLUMNS(F, DR, LAMBDA, OPTS)
%   [f, INFO] = ABEL_INVERT_COLUMNS(...)
%
%   F    : [N x Ns] projected (line-of-sight integrated) profiles. Row 1
%          is the symmetry axis (r = 0), row N the outer edge. Each COLUMN
%          is one independent axial station.
%   DR   : radial grid spacing.
%   LAMBDA : how the Tikhonov (second-difference) strength is chosen.
%          'gcv'            (default) one lambda for the whole field, at
%                           the minimum of the generalised cross-validation
%                           function summed over all valid stations
%          'lcurve'         one lambda for the whole field, at the point
%                           of maximum curvature of the L-curve. See the
%                           caveat below before using it.
%          'lcurve-station' a separate L-curve corner for every station
%          numeric scalar   fixed strength, relative to norm(A)/norm(L),
%                           as in earlier versions. 0 gives a plain,
%                           unregularised onion-peel.
%
%   OPTS (struct, all optional, used by the automatic modes only)
%     opts.nLambda  : points on the log-spaced lambda grid, default 200
%     opts.relRange : [min max] of the grid, relative to norm(A)/norm(L),
%                     default [1e-4 1e3]
%
%   Returns f : [N x Ns] the LOCAL (radially resolved) quantity.
%
%   INFO struct
%     info.mode        : 'fixed' | 'gcv' | 'lcurve' | 'lcurve-station'
%     info.lambda      : lambda used. Scalar, or [1 x Ns] (NaN on stations
%                        that were not inverted) for 'lcurve-station'
%     info.lambdaRel   : the same, divided by info.lambdaScale, so it can
%                        be compared with the old fixed setting of 1
%     info.lambdaScale : norm(A,'fro')/norm(L,'fro')
%     info.lambdaGrid  : [nLambda x 1] grid the criterion was evaluated on
%     info.resNorm     : ||A f - F||, [nLambda x 1] summed over stations,
%                        or [nLambda x nGood] for 'lcurve-station'
%     info.semiNorm    : ||L f||, same size as resNorm
%     info.gcv         : [nLambda x 1] GCV function
%     info.curvature   : curvature of the log-log L-curve, same size as
%                        resNorm
%     info.iChosen     : index into lambdaGrid of the chosen lambda
%     info.stations    : column indices of F that were inverted
%
%   METHOD
%   Onion-peeling: the medium is modelled as N concentric annuli of
%   constant value. A ray with impact parameter y_i = (i-1)*dr crosses
%   annulus j over a chord of length
%       2*( sqrt(r_j^2 - y_i^2) - sqrt(r_{j-1}^2 - y_i^2) )
%   giving an upper-triangular system A*f = F (Dasch 1992).
%
%   Onion-peeling is exact for noise-free input but strongly noise-
%   amplifying, because each annulus is recovered by subtracting the
%   contributions of all the annuli outside it. It is therefore solved in
%   Tikhonov form (Daun et al. 2006),
%       minimise  ||A f - F||^2 + lambda^2 ||L f||^2
%   with L the second-difference operator, as the least-squares problem
%   [A; lambda*L] f = [F; 0], by economy QR at each lambda.
%
%   CHOOSING LAMBDA: GENERALISED CROSS-VALIDATION
%   GCV (Golub, Heath & Wahba 1979) picks the lambda whose solution best
%   predicts a data point left out of the fit, without needing an estimate
%   of the noise level. For one station
%       G(lambda) = ||A f - F||^2 / (N - trace(H))^2
%   where H = A (A'A + lambda^2 L'L)^-1 A' maps the data onto the fitted
%   data. With [A; lambda*L] = Q*R and Q1 the first N rows of Q, H equals
%   Q1*Q1', so trace(H) = ||Q1||_F^2 costs nothing extra. Every station
%   shares A and L, so the whole field is one problem with the residuals
%   summed over stations and the same trace for each. One shared lambda
%   also means streamwise changes in the result reflect the data rather
%   than a station-to-station change in smoothing.
%
%   GCV assumes uncorrelated noise. The S that BOS_OPTICAL_PATH produces
%   is an integral of the deflection, so its noise is correlated along r.
%   A synthetic test with noise added to the deflection BEFORE that
%   integration still put GCV within 5% of the lowest achievable error.
%
%   WHY NOT THE L-CURVE BY DEFAULT
%   The L-curve (Hansen & O'Leary 1993; recommended for flame tomography
%   by Akesson & Daun 2008) takes the corner of log||A f - F|| against
%   log||L f||. That corner exists when the residual levels off at the
%   noise floor for small lambda. Here A is square and invertible, so the
%   residual keeps falling to zero instead, and onion-peeling is only
%   mildly ill-posed. The curve then has no sharp corner. On the 8 bar
%   data its curvature has several weak peaks, and in synthetic Gaussian
%   tests the maximum-curvature lambda gave up to about 3x the lowest
%   achievable error, where GCV stayed within 30%. The L-curve modes
%   remain available for comparison.
%
%   A chosen lambda at either end of the grid means the optimum may lie
%   outside it: the function warns, and you should widen opts.relRange.
%
%   Columns containing any NaN are returned as NaN: an Abel inversion
%   needs a complete radial profile, and silently interpolating across a
%   gap would smear an unknown amount of signal inward.
%
%   REFERENCES
%   C. J. Dasch, Appl. Opt. 31(8), 1146-1152 (1992).
%   K. J. Daun, K. A. Thomson, F. Liu, G. J. Smallwood, Appl. Opt. 45(19),
%       4638-4646 (2006).
%   G. H. Golub, M. Heath, G. Wahba, Technometrics 21(2), 215-223 (1979).
%   P. C. Hansen, D. P. O'Leary, SIAM J. Sci. Comput. 14(6), 1487-1503
%       (1993).
%   E. O. Akesson, K. J. Daun, Appl. Opt. 47(3), 407-416 (2008).

    if nargin < 3 || isempty(lambda), lambda = 'gcv'; end
    if nargin < 4 || isempty(opts),   opts = struct(); end
    if ~isfield(opts,'nLambda')  || isempty(opts.nLambda),  opts.nLambda  = 200; end
    if ~isfield(opts,'relRange') || isempty(opts.relRange), opts.relRange = [1e-4 1e3]; end

    [N, Ns] = size(F);
    if N < 3
        error('abel_invert_columns:tooFewRows', ...
            'Need at least 3 radial points, got %d.', N);
    end

    % --- forward (onion-peel) matrix -----------------------------------
    A = zeros(N, N);
    for i = 1:N
        yi = (i-1)*dr;
        for j = i:N
            outer = sqrt(max((j*dr)^2     - yi^2, 0));
            inner = sqrt(max(((j-1)*dr)^2 - yi^2, 0));
            A(i,j) = 2*(outer - inner);
        end
    end

    % --- smoothness penalty: second difference --------------------------
    L = full(spdiags(repmat([1 -2 1], N-2, 1), 0:2, N-2, N));
    lamScale = norm(A,'fro') / norm(L,'fro');

    f = nan(N, Ns);
    good = ~any(isnan(F), 1);
    info = struct('mode', '', 'lambda', NaN, 'lambdaRel', NaN, ...
                  'lambdaScale', lamScale, 'lambdaGrid', [], ...
                  'resNorm', [], 'semiNorm', [], 'gcv', [], ...
                  'curvature', [], 'iChosen', [], 'stations', find(good));
    if ~any(good), return; end
    Fg = F(:, good);

    % --- fixed strength (old behaviour) ---------------------------------
    if isnumeric(lambda)
        info.mode      = 'fixed';
        info.lambda    = lambda * lamScale;
        info.lambdaRel = lambda;
        f(:, good) = local_solve(A, L, info.lambda, Fg);
        return
    end

    mode = lower(lambda);
    if ~any(strcmp(mode, {'gcv', 'lcurve', 'lcurve-station'}))
        error('abel_invert_columns:badLambda', ...
            'LAMBDA must be numeric, ''gcv'', ''lcurve'' or ''lcurve-station'', got ''%s''.', lambda);
    end
    info.mode = mode;

    % --- evaluate residual, seminorm and trace(H) on the grid -------------
    lamGrid = lamScale * logspace(log10(opts.relRange(1)), ...
                                  log10(opts.relRange(2)), opts.nLambda).';
    nL   = numel(lamGrid);
    nG   = size(Fg, 2);
    trH  = zeros(nL, 1);
    if strcmp(mode, 'lcurve-station')
        res = zeros(nL, nG);
        sem = zeros(nL, nG);
        for q = 1:nL
            [fq, trH(q)] = local_solve(A, L, lamGrid(q), Fg);
            res(q,:) = sqrt(sum((A*fq - Fg).^2, 1));
            sem(q,:) = sqrt(sum((L*fq).^2, 1));
        end
        resG = sqrt(sum(res.^2, 2));
        semG = sqrt(sum(sem.^2, 2));
    else
        % Field-wide norms only. The solution is f = M*F with M an N x N
        % matrix, so summed over stations
        %     ||(A*M - I) F||_F^2 = trace(P*G*P'),  P = A*M - I
        %     ||L*M*F||_F^2       = trace(Z*G*Z'),  Z = L*M
        % with the Gram matrix G = F*F' formed once. Each lambda then costs
        % O(N^3) however many stations there are.
        G    = Fg * Fg.';
        resG = zeros(nL, 1);
        semG = zeros(nL, 1);
        for q = 1:nL
            [M, trH(q)] = local_solve(A, L, lamGrid(q), eye(N));
            P = A*M - eye(N);
            Z = L*M;
            resG(q) = sqrt(max(sum((P*G).*P, 'all'), 0));
            semG(q) = sqrt(max(sum((Z*G).*Z, 'all'), 0));
        end
    end
    info.lambdaGrid = lamGrid;
    info.gcv        = resG.^2 ./ (N - trH).^2;

    switch mode
        case 'gcv'
            [~, iC] = min(info.gcv);
            [~, kap] = local_corner(lamGrid, resG, semG);
        case 'lcurve'
            [iC, kap] = local_corner(lamGrid, resG, semG);
        case 'lcurve-station'
            iC  = zeros(1, nG);
            kap = zeros(nL, nG);
            for k = 1:nG
                [iC(k), kap(:,k)] = local_corner(lamGrid, res(:,k), sem(:,k));
            end
    end

    nEdge = nnz(iC <= 2 | iC >= nL-1);
    if nEdge > 0
        warning('abel_invert_columns:atGridEdge', ...
            ['%d of %d chosen lambda values lie at the edge of the grid ' ...
             '(lambdaRel %.3g to %.3g). Widen opts.relRange.'], ...
             nEdge, numel(iC), opts.relRange(1), opts.relRange(2));
    end
    info.curvature = kap;
    info.iChosen   = iC;

    if strcmp(mode, 'lcurve-station')
        % stations sharing a lambda share one factorisation
        fg = zeros(N, nG);
        for q = unique(iC)
            sel = (iC == q);
            fg(:, sel) = local_solve(A, L, lamGrid(q), Fg(:, sel));
        end
        f(:, good) = fg;
        lamAll = nan(1, Ns);
        lamAll(good)   = lamGrid(iC);
        info.lambda    = lamAll;
        info.lambdaRel = lamAll / lamScale;
        info.resNorm   = res;
        info.semiNorm  = sem;
    else
        info.lambda    = lamGrid(iC);
        info.lambdaRel = lamGrid(iC) / lamScale;
        info.resNorm   = resG;
        info.semiNorm  = semG;
        f(:, good) = local_solve(A, L, info.lambda, Fg);
    end
end

% ======================================================================
function [X, trH] = local_solve(A, L, lam, B)
%Tikhonov solution for every column of B at one lambda.
%   Economy QR of the stacked matrix [A; lam*L]. The stacked right-hand
%   side is [B; 0], so only the first N rows of Q, Q1, contribute to
%   Q'*[B; 0]. The influence matrix is Q1*Q1', so its trace is ||Q1||_F^2.
    N = size(A, 1);
    if lam == 0
        X = A \ B;
        trH = N;
        return
    end
    [Q, R] = qr([A; lam*L], 0);
    Q1  = Q(1:N, :);
    X   = R \ (Q1.' * B);
    trH = sum(Q1.^2, 'all');
end

% ======================================================================
function [iC, kappa] = local_corner(lam, rho, eta)
%Index of maximum curvature of the log-log L-curve (log rho, log eta).
%   Parametrised by t = log(lambda). As lambda grows a textbook L-curve
%   runs down its vertical arm and turns right along its horizontal arm, a
%   counter-clockwise turn, so the corner is the maximum of the signed
%   curvature
%       kappa = (x' y'' - x'' y') / (x'^2 + y'^2)^(3/2).
    t = log(lam(:));
    x = log(max(rho(:), realmin));
    y = log(max(eta(:), realmin));
    dx  = gradient(x, t);   dy  = gradient(y, t);
    ddx = gradient(dx, t);  ddy = gradient(dy, t);
    kappa = (dx.*ddy - ddx.*dy) ./ (dx.^2 + dy.^2).^1.5;
    kappa([1 end]) = NaN;          % one-sided differences at the ends
    [~, iC] = max(kappa);
end
