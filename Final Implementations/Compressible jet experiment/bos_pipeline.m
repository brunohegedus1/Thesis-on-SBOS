function R = bos_pipeline(D, opts)
%BOS_PIPELINE  Straighten a tilted/curved jet and fold it about its axis.
%
%   R = BOS_PIPELINE(D)
%   R = BOS_PIPELINE(D, OPTS)
%
%   D is the struct returned by READ_BOS_TXT (fields x, y, U, V).
%   Executes the seven processing steps:
%
%     1. remove background offset          (vector offset, see note below)
%     2. detect the jet centerline
%     3. fit a smooth spline to it
%     4. resample into spline-aligned coordinates (straighten)
%     5. crop symmetrically about the centerline
%     6. smooth the radial profiles
%     7. optionally fold upper/lower halves
%
%   There is deliberately NO Abel inversion here. A BOS displacement is a
%   deflection, the transverse derivative of the projected field, so it is
%   not a plain Abel projection and must not be inverted directly. The
%   density route integrates it first: see RUN_BOS_DENSITY, which passes
%   R.Dn_f to BOS_OPTICAL_PATH and then ABEL_INVERT_COLUMNS.
%
%   IMPORTANT -- STEP 1 IS A VECTOR SUBTRACTION.
%   The background offset in these exports is a near-uniform displacement
%   VECTOR (in an earlier dataset roughly (-0.0107, -0.0135) mm), not a
%   scalar pedestal on the magnitude. Its length happens to equal the
%   median magnitude, so subtracting it from |d| looks plausible but is
%   wrong: it corrupts the direction of every vector. This function
%   therefore estimates and subtracts an offset from U and V separately.
%
%   STEP 1 USES A FIXED WINDOW OF QUIESCENT AIR.
%   The offset is the median of U and V inside a rectangle that must lie
%   entirely outside the jet, set by opts.bgWindow and computed by
%   BOS_BACKGROUND_OFFSET, which PLOT_DISPLACEMENT_HISTOGRAM also uses. The
%   default is the window below the jet, towards the lower-left of the
%   valid mask, chosen by inspection of the 3 cm exports. It replaces an
%   automatic estimate
%   that excluded a band around a provisional centreline. That estimate
%   also discarded every column too weak to locate the axis in, which are
%   the quietest columns of all. Check the window whenever the field of
%   view, the nozzle position or the jet direction changes.
%
%   IMPORTANT -- STEP 2 DEFAULT IS NOT "MAXIMUM DISPLACEMENT".
%   For an axisymmetric jet the cross-stream BOS displacement is an ODD
%   function of the distance from the axis: it passes through zero ON the
%   axis and peaks on either side. The magnitude therefore has a local
%   MINIMUM on the axis (visible as the dark centres of the shock-diamond
%   rings), so a per-column maximum lands on one flank, not the axis, and
%   flips sides erratically. On this dataset:
%       'argmax'    scatter about a straight fit ~2.79 mm
%       'centroid'  scatter ~0.70 mm
%       'zerocross' scatter ~0.24 mm   <-- default
%   'argmax' and 'centroid' are still selectable via opts.centerMethod.
%
%   OPTS (struct, all optional)
%    -- step 1 --
%     opts.bgWindow : [x1 y1; x2 y2], two opposite corners of the
%                     background rectangle, in the units of D.x and D.y.
%                     Default [] = the default of BOS_BACKGROUND_OFFSET,
%                     [3.244 14.9421; -19.27 4.536].
%     opts.bgMode   : 'median' (default) | 'mean' | 'none'
%    -- step 2 --
%     opts.centerMethod : 'zerocross' (default) | 'centroid' | 'argmax'
%     opts.centerBand   : half-width of the search band about the running
%                         estimate, default 3
%     opts.centerIter   : refinement iterations, default 4
%     opts.qualMinPeak  : columns whose peak |d| after background removal
%                         is below this are ignored as too weak. Default
%                         [] = auto (40th percentile of column peaks).
%    -- step 3 --
%     opts.nKnots    : interior knots for the least-squares spline,
%                      default 5. Keep this SMALL. A smoothing spline
%                      driven only by a residual target will satisfy the
%                      residual while oscillating violently between
%                      knots: on this dataset that produced slopes of
%                      +/-65 (true slope -0.28) and an arc length of
%                      145 mm across a 38 mm span. Knot count is the
%                      stable control.
%     opts.clipSigma : robust outlier rejection threshold, default 2.5
%     opts.clipIter  : rejection iterations, default 5
%    -- steps 4/5 --
%     opts.ds    : streamwise spacing of the straightened grid,
%                  default = median grid spacing
%     opts.dn    : cross-stream spacing, default = median grid spacing
%     opts.nMax  : half-height of the symmetric crop, default 6
%    -- step 6 --
%     opts.smoothWin : moving-average window in POINTS along n, applied
%                      only across the radial direction, default 5.
%                      Use 1 to disable.
%    -- step 7 --
%     opts.foldHalves : true (default) to average the two halves.
%
%   OUTPUT struct R
%     R.s, R.n        : straightened coordinates (n includes exactly 0)
%     R.Ds, R.Dn      : straightened streamwise / cross-stream shifts
%     R.Ds_sm,R.Dn_sm : after radial smoothing
%     R.nHalf         : radial coordinate of the folded half (>= 0)
%     R.Ds_f, R.Dn_f  : folded profiles (upper half only if opts.foldHalves
%                       is false). Still PROJECTED, line-of-sight integrated.
%     R.centerX, R.centerY : detected centerline points used for the fit
%     R.spline        : the fitted spline (MATLAB spline struct)
%     R.theta         : local tangent angle at each station [rad]
%     R.bgU, R.bgV    : background offsets removed
%     R.info          : diagnostics, including bgWindow (the window used, []
%                       if opts.bgMode is 'none') and bgPoints (the number
%                       of valid vectors the offset came from)

    if nargin < 2 || isempty(opts), opts = struct(); end
    dxg = median(diff(D.x));  dyg = median(diff(D.y));
    if ~isfield(opts,'bgWindow'), opts.bgWindow = []; end
    opts = local_default(opts,'bgMode','median');
    opts = local_default(opts,'centerMethod','zerocross');
    opts = local_default(opts,'centerBand',3);
    opts = local_default(opts,'centerIter',4);
    opts = local_default(opts,'qualMinPeak',[]);
    opts = local_default(opts,'nKnots',5);
    opts = local_default(opts,'clipSigma',2.5);
    opts = local_default(opts,'clipIter',5);
    opts = local_default(opts,'ds',median([dxg dyg]));
    opts = local_default(opts,'dn',median([dxg dyg]));
    opts = local_default(opts,'nMax',6);
    opts = local_default(opts,'smoothWin',5);
    opts = local_default(opts,'foldHalves',true);

    x = D.x(:).';  y = D.y(:);
    U = D.U;       V = D.V;
    [Ny, Nx] = size(U);

    %% ---------------- STEP 1: background offset (vector) --------------
    % Median of U and V over a fixed window of quiescent air.
    if strcmpi(opts.bgMode,'none')
        bgU = 0;  bgV = 0;  bgN = 0;  bgWin = [];
        fprintf('step 1: background offset removal disabled\n');
    else
        [bgU, bgV, bgN, bgWin] = bos_background_offset(x, y, U, V, ...
                                     opts.bgWindow, opts.bgMode);
        xw = sort(bgWin(:,1));  yw = sort(bgWin(:,2));
        fprintf(['step 1: background offset removed  (u0 = %+.5f, v0 = %+.5f) ' ...
                 'from %d vectors in x [%.2f, %.2f], y [%.2f, %.2f]\n'], ...
                 bgU, bgV, bgN, xw(1), xw(2), yw(1), yw(2));
    end
    Ub = U - bgU;  Vb = V - bgV;

    %% ---------------- STEP 2: centerline detection --------------------
    if opts.centerlineDetection
        cen = local_centerline(Ub, Vb, x, y, opts.centerMethod, opts);
        okc = ~isnan(cen);
        if nnz(okc) < 20
            error('bos_pipeline:noCenterline', ...
                ['Centerline detection found only %d usable columns. Check that ' ...
                 'the jet runs roughly left-to-right and that opts.qualMinPeak ' ...
                 'is not too strict.'], nnz(okc));
        end
        cx = x(okc).';  cy = cen(okc).';
        fprintf('step 2: centerline from %d columns (method: %s)\n', numel(cx), opts.centerMethod);
    else

        %cx = [17.2876, 13.343,-9.648, -1.77, -5.96, -13.24, -19.502]; %Points obtained by visual inspection of 3cm 8bar 12x12
        %cy = [13.56, 14.666, 20.88, 18.5796, 19.78, 21.849, 23.50];

        cx = [17.81, 17.2416,-9.648, -1.77, -5.96, -13.24, -19.502]; %Points obtained by visual inspection of 3cm 5bar 12x12
        cy = [13.31, 13.377, 20.88, 18.5796, 19.78, 21.849, 23.50];


    end
    
     %% ---------------- STEP 3: robust smooth spline --------------------
     if opts.SplineFit
        [sp, keep] = local_fit_spline(cx, cy, opts);
        cx = cx(keep);  cy = cy(keep);
        resid = cy - ppval(sp, cx);
        fprintf('step 3: spline fit, %d pts kept, residual sd = %.4f\n', numel(cx), std(resid));
    else
    %% ---------------- STEP 3: linear regression ----------------------
    % Robust line fit
    
        p = local_robust_polyfit(cx, cy);
        
        x0 = min(cx);
        breaks = [x0 max(cx)];
        
        y0 = polyval(p,x0);
        
        sp = mkpp(breaks,[0 0 p(1) y0]);
    
        keep = true(size(cx));
    
        resid = cy - polyval(p, cx);
    
        fprintf('step 3: linear regression, residual sd = %.4f\n', std(resid));
    end 
    %% ------------ STEPS 4 & 5: straighten + symmetric crop ------------
    [s_grid, n_grid, Ds, Dn, theta] = local_straighten(Ub, Vb, x, y, sp, ...
                                        min(cx), max(cx), opts);
    fprintf('step 4/5: straightened onto %d x %d grid (s: %.2f..%.2f, n: +/-%.2f)\n', ...
        numel(n_grid), numel(s_grid), s_grid(1), s_grid(end), opts.nMax);

    %% ---------------- STEP 6: light radial smoothing ------------------
    if opts.RadialSmoothing
        Ds_sm = local_smooth_radial(Ds, opts.smoothWin);
        Dn_sm = local_smooth_radial(Dn, opts.smoothWin);
        fprintf('step 6: radial smoothing, window = %d points\n', opts.smoothWin);
    else
        Ds_sm = Ds;
        Dn_sm = Dn;
        fprintf('step 6: radial smoothing, deactivated \n');
    end

    %% ---------------- STEP 7: fold the two halves ---------------------
    cIdx = find(n_grid == 0, 1);
    if isempty(cIdx)
        [~, cIdx] = min(abs(n_grid));
    end
    M = min(cIdx-1, numel(n_grid)-cIdx);
    nHalf = n_grid(cIdx:cIdx+M);
    if opts.foldHalves
        % Dn is ODD about the axis, Ds is EVEN -- fold with the right parity
        Dn_f = local_fold(Dn_sm, cIdx, M, true);
        Ds_f = local_fold(Ds_sm, cIdx, M, false);
        fprintf('step 7: folded halves (Dn antisymmetric, Ds symmetric)\n');
    else
        Dn_f = Dn_sm(cIdx:cIdx+M, :);
        Ds_f = Ds_sm(cIdx:cIdx+M, :);
        fprintf('step 7: folding disabled, using upper half only\n');
    end

    %% ---------------- pack ---------------------------------------------
    R = struct();
    R.s = s_grid;  R.n = n_grid;  R.nHalf = nHalf;
    R.Ds = Ds;     R.Dn = Dn;
    R.Ds_sm = Ds_sm;  R.Dn_sm = Dn_sm;
    R.Ds_f = Ds_f;    R.Dn_f = Dn_f;
    R.centerX = cx;  R.centerY = cy;
    R.spline = sp;   R.theta = theta;
    R.bgU = bgU;     R.bgV = bgV;
    R.info = struct('centerResidStd', std(resid), ...
                    'meanTiltDeg', mean(theta)*180/pi, ...
                    'nStations', numel(s_grid), ...
                    'bgWindow', bgWin, ...
                    'bgPoints', bgN, ...
                    'opts', opts);
end

% ======================================================================
function s = local_default(s,f,v)
    if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end


% ======================================================================
function cen = local_centerline(U, V, x, y, method, opts)
%Per-column estimate of the jet axis position.
    [Ny, Nx] = size(U);
    Mag = hypot(U, V);

    % column quality: reject columns whose signal is too weak to localise
    colPeak = max(Mag, [], 1);   % max already ignores NaN
    if isempty(opts.qualMinPeak)
        qmin = local_prctile(colPeak(~isnan(colPeak)), 40);
    else
        qmin = opts.qualMinPeak;
    end
    qual = colPeak > qmin;

    % Provisional straight-line estimate from a magnitude-weighted centroid.
    thr = local_prctile(Mag(~isnan(Mag)), 90);
    Wt  = Mag - thr;  Wt(isnan(Wt) | Wt<0) = 0;
    cen = nan(1, Nx);
    for j = 1:Nx
        w = Wt(:,j);
        if sum(w) > 0, cen(j) = sum(y.*w)/sum(w); end
    end
    ok = ~isnan(cen);
    p  = local_robust_polyfit(x(ok), cen(ok));

    % Always converge the CENTROID estimate first, whatever the requested
    % method. 'zerocross' picks the sign change nearest the running guess,
    % so seeding it from a crude guess lets it lock onto the wrong crossing
    % and then re-derive the search angle from its own bad output -- a
    % feedback loop that converges confidently to the wrong axis. On this
    % dataset that produced a -8.8 deg axis instead of the true -15.9 deg.
    % Converging the centroid first removes that failure mode.
    for it = 1:max(opts.centerIter,3)
        pred = polyval(p, x);
        cenC = nan(1, Nx);
        for j = 1:Nx
            if ~qual(j), continue; end
            band = abs(y - pred(j)) < opts.centerBand;
            if nnz(band) < 8, continue; end
            w = Wt(:,j);  w(~band) = 0;
            if sum(w) > 0, cenC(j) = sum(y.*w)/sum(w); end
        end
        okc = ~isnan(cenC);
        if nnz(okc) < 10, break; end
        p = local_robust_polyfit(x(okc), cenC(okc));
        cen = cenC;
    end

    % Now refine with the requested method, seeded from the converged fit.
    for it = 1:opts.centerIter
        pred = polyval(p, x);
        cen  = nan(1, Nx);
        for j = 1:Nx
            if ~qual(j), continue; end
            band = abs(y - pred(j)) < opts.centerBand;
            if nnz(band) < 8, continue; end
            switch lower(method)
                case 'argmax'
                    col = Mag(:,j);  col(~band) = NaN;
                    if all(isnan(col)), continue; end
                    [~, k] = max(col);
                    cen(j) = y(k);

                case 'centroid'
                    w = Wt(:,j);  w(~band) = 0;
                    if sum(w) > 0, cen(j) = sum(y.*w)/sum(w); end

                case 'zerocross'
                    % cross-stream component is odd about the axis: find
                    % where it changes sign, nearest the running estimate
                    th = atan(p(1));
                    Dn = -U(:,j)*sin(th) + V(:,j)*cos(th);
                    yy = y(band);  dd = Dn(band);
                    g  = ~isnan(dd);
                    yy = yy(g);  dd = dd(g);
                    if numel(dd) < 8, continue; end
                    sc = find(diff(sign(dd)) ~= 0);
                    if isempty(sc), continue; end
                    zc = nan(numel(sc),1);
                    for q = 1:numel(sc)
                        i1 = sc(q);  f1 = dd(i1);  f2 = dd(i1+1);
                        if f2 ~= f1
                            zc(q) = yy(i1) - f1*(yy(i1+1)-yy(i1))/(f2-f1);
                        end
                    end
                    zc = zc(~isnan(zc));
                    if isempty(zc), continue; end
                    [dmin, k] = min(abs(zc - pred(j)));
                    % ignore crossings implausibly far from the estimate
                    if dmin < opts.centerBand
                        cen(j) = zc(k);
                    end

                otherwise
                    error('bos_pipeline:badMethod', ...
                        'Unknown centerMethod "%s".', method);
            end
        end
        ok = ~isnan(cen);
        if nnz(ok) < 10, break; end
        p = local_robust_polyfit(x(ok), cen(ok));
    end
end

% ======================================================================
function p = local_robust_polyfit(x, v)
%Straight-line fit with iterative sigma clipping (portable, no toolbox).
    x = x(:);  v = v(:);
    keep = true(size(x));
    p = polyfit(x, v, 1);
    for it = 1:5
        r  = v - polyval(p, x);
        sd = std(r(keep));
        if sd <= 0, break; end
        nk = abs(r) < 2.5*sd;
        if nnz(nk) < 10 || isequal(nk, keep), break; end
        keep = nk;
        p = polyfit(x(keep), v(keep), 1);
    end
end

% ======================================================================
function q = local_prctile(v, pct)
%Percentile without the Statistics Toolbox.
    v = sort(v(:));
    if isempty(v), q = NaN; return; end
    if numel(v) == 1, q = v; return; end
    pos = 1 + (numel(v)-1)*pct/100;
    lo  = floor(pos);  hi = ceil(pos);
    if lo == hi
        q = v(lo);
    else
        q = v(lo) + (pos-lo)*(v(hi)-v(lo));
    end
end

% ======================================================================
function [sp, keep] = local_fit_spline(cx, cy, opts)
%Least-squares cubic spline with FEW interior knots, plus sigma clipping.
    keep = true(size(cx));
    sp = [];
    for it = 1:opts.clipIter
        xk = cx(keep);  yk = cy(keep);
        if numel(xk) < opts.nKnots + 6
            error('bos_pipeline:tooFewPts', ...
                'Not enough centerline points (%d) for %d knots.', ...
                numel(xk), opts.nKnots);
        end
        kn = linspace(min(xk), max(xk), opts.nKnots+2);
        kn = kn(2:end-1);
        if exist('spap2','file') == 2
            spf = spap2(augknt([min(xk) kn max(xk)], 4), 4, xk, yk);
            sp  = fn2fm(spf, 'pp');
        else
            sp = local_lsq_cubic_spline(xk, yk, kn);
        end
        r  = cy - ppval(sp, cx);
        sd = std(r(keep));
        newKeep = abs(r) < opts.clipSigma*sd;
        if isequal(newKeep, keep), break; end
        keep = newKeep;
    end
end

% ======================================================================
function pp = local_lsq_cubic_spline(xk, yk, interiorKnots)
%Least-squares cubic spline via a truncated-power basis (no toolbox needed).
%   The basis is built on a CENTRED AND SCALED abscissa. The raw truncated
%   power basis is badly conditioned when x is large and off-centre -- on
%   this dataset cond(A) was 1.3e6 raw versus 3.6e4 scaled -- and the
%   resulting coefficient error showed up as a visibly oscillating
%   centreline fit. The pp is converted back to the original x at the end.
    xk = xk(:);  yk = yk(:);
    x0 = mean(xk);
    xsc = std(xk);
    if xsc <= 0, xsc = 1; end
    xs = (xk - x0)/xsc;
    t  = (interiorKnots(:).' - x0)/xsc;

    A  = [ones(numel(xs),1), xs, xs.^2, xs.^3];
    for q = 1:numel(t)
        A = [A, max(xs - t(q), 0).^3]; %#ok<AGROW>
    end
    c = A \ yk;

    % piecewise-polynomial in the SCALED variable
    brs = [min(xs), t, max(xs)];
    nb  = numel(brs)-1;
    coefs = zeros(nb,4);
    for b = 1:nb
        a0 = c(1); a1 = c(2); a2 = c(3); a3 = c(4);
        for q = 1:numel(t)
            if brs(b) >= t(q)   % this truncated term is active on the piece
                tq = t(q);  cq = c(4+q);
                a0 = a0 + cq*(-tq)^3;
                a1 = a1 + cq*3*tq^2;
                a2 = a2 + cq*(-3*tq);
                a3 = a3 + cq;
            end
        end
        L = brs(b);          % re-centre the piece on its left breakpoint
        d0 = a0 + a1*L + a2*L^2 + a3*L^3;
        d1 = a1 + 2*a2*L + 3*a3*L^2;
        d2 = a2 + 3*a3*L;
        d3 = a3;
        % undo the x-scaling so the pp is valid in the ORIGINAL variable:
        % a local offset dx_orig corresponds to dx_scaled = dx_orig/xsc
        coefs(b,:) = [d3/xsc^3, d2/xsc^2, d1/xsc, d0];
    end
    br = brs*xsc + x0;
    pp = mkpp(br, coefs);
end

% ======================================================================
function [s_grid, n_grid, Ds, Dn, theta_s] = local_straighten(U, V, x, y, sp, xlo, xhi, opts)
%Resample onto a spline-aligned grid, rotating vectors by the LOCAL angle.
    % arc length along the spline
    xf   = linspace(xlo, xhi, 4000).';
    dsp  = fnder_pp(sp);
    dyf  = ppval(dsp, xf);
    integ= sqrt(1 + dyf.^2);
    sOf  = [0; cumsum(0.5*(integ(2:end)+integ(1:end-1)).*diff(xf))];
    sTot = sOf(end);

    s_grid = (0:opts.ds:sTot);
    x_s    = interp1(sOf, xf, s_grid, 'linear', 'extrap');
    y_s    = ppval(sp, x_s);
    theta_s= atan(ppval(dsp, x_s));       % local tangent angle

    % symmetric cross-stream grid that contains EXACTLY zero
    M      = floor(opts.nMax/opts.dn);
    n_grid = (-M:M).'*opts.dn;

    [S_TH, ~] = meshgrid(theta_s, n_grid); %#ok<ASGLU>
    TH = repmat(theta_s(:).', numel(n_grid), 1);
    XS = repmat(x_s(:).',     numel(n_grid), 1);
    YS = repmat(y_s(:).',     numel(n_grid), 1);
    NN = repmat(n_grid,       1, numel(s_grid));
    %TH = -15*pi/180;
    
    % unit normal is (-sin(theta), cos(theta))
    Xq = XS - NN.*sin(TH);
    Yq = YS + NN.*cos(TH);

    [XX, YY] = meshgrid(x, y);
    Uq = interp2(XX, YY, U, Xq, Yq, 'linear', NaN);
    Vq = interp2(XX, YY, V, Xq, Yq, 'linear', NaN);

    % rotate the VECTOR COMPONENTS into the local tangential/normal frame.
    % This must use the local theta, not one global angle -- that is the
    % whole point of following a curved spline.
    Ds =  Uq.*cos(TH) + Vq.*sin(TH);
    Dn = -Uq.*sin(TH) + Vq.*cos(TH);
end

% ======================================================================
function dpp = fnder_pp(pp)
%Derivative of a pp-form spline (avoids requiring the Curve Fitting Toolbox).
    [br, cf, L, K, dd] = unmkpp(pp); %#ok<ASGLU>
    if K <= 1
        dpp = mkpp(br, zeros(L,1));  return
    end
    newc = cf(:,1:K-1) .* repmat((K-1):-1:1, L, 1);
    dpp  = mkpp(br, newc);
end

% ======================================================================
function Q = local_smooth_radial(Q, win)
%NaN-aware moving average along the radial (first) dimension only.
    if win <= 1, return; end
    nanm  = isnan(Q);
    filled= Q;  filled(nanm) = 0;
    k     = ones(win,1)/win;
    num   = conv2(filled,        k, 'same');
    den   = conv2(double(~nanm), k, 'same');
    Q     = num ./ den;
    Q(den == 0) = NaN;
end

% ======================================================================
function F = local_fold(Q, cIdx, M, antisym)
%Average the two halves about the axis with the correct parity.
    pos = Q(cIdx:cIdx+M, :);
    neg = Q(cIdx:-1:cIdx-M, :);
    if antisym, neg = -neg; end
    F = nan(size(pos));
    both = ~isnan(pos) & ~isnan(neg);
    onlyP= ~isnan(pos) &  isnan(neg);
    onlyN=  isnan(pos) & ~isnan(neg);
    F(both) = 0.5*(pos(both) + neg(both));
    F(onlyP)= pos(onlyP);
    F(onlyN)= neg(onlyN);
end
