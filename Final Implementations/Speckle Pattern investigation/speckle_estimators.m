function f = speckle_estimators()
%SPECKLE_ESTIMATORS  Handles to the speckle-size routines shared by
% speckle_size.m and validate_speckle_size.m, so that the validation
% exercises the same code that produces the reported numbers.
%
%   f = speckle_estimators();
%   [acf, lagx, lagy]     = f.acf2(I, maxLag)                2D normalized ACF
%   [acf, lag]            = f.acf1(v, maxLag, segLen)        1D normalized ACF
%   [aMean, aStd, nUsed]  = f.avg1d(A, maxLag, thr)          average 1D method
%   w                     = f.widths(lag, prof, thr)         threshold crossings
%   a                     = f.cross(lag, p, level)           one crossing
%   B                     = f.gaussblur(A, sigma)            Gaussian low pass
%
% See speckle_size.m for the method and its reference.

    f.acf2      = @local_acf2;
    f.acf1      = @local_acf1;
    f.avg1d     = @local_avg1d;
    f.widths    = @local_widths;
    f.cross     = @local_cross;
    f.gaussblur = @local_gaussblur;
end

function [acf, lagx, lagy] = local_acf2(I, maxLag)
% 2D normalized autocorrelation of intensity via Wiener-Khinchin.
% Zero padding avoids circular wrap-around; normalizing by the zero-lag value
% of the mean-subtracted field reproduces eq. (1) with r_I(0,0) = 1. The
% overlap correction below removes the triangular bias of the FFT estimator,
% which would otherwise scale r_I at lag d by (1-d/n) in each direction.
    [ny, nx] = size(I);
    Ic = single(I) - mean(I(:));
    Ny = 2*ny;  Nx = 2*nx;                       % zero padded size
    F  = fft2(Ic, Ny, Nx);
    C  = real(ifft2(F .* conj(F)));
    C  = C / C(1,1);
    C  = fftshift(C);
    cy = Ny/2 + 1;  cx = Nx/2 + 1;
    L  = min([maxLag, ny-1, nx-1]);
    acf  = double(C(cy-L:cy+L, cx-L:cx+L));
    lagx = -L:L;  lagy = -L:L;

    wx = 1 - abs(lagx)/nx;
    wy = 1 - abs(lagy)/ny;
    acf = acf ./ (wy(:) * wx);
end

function [acf, lag] = local_acf1(v, maxLag, segLen)
% 1D normalized autocorrelation of a (long) intensity vector, lags 0..maxLag.
%
% SEGLEN is the length of the contiguous pieces the vector is made of; for
% the long-vector method of Hamarova et al. that is one image row (column),
% not the whole vector. It matters because the pairs that straddle a joint
% between two rows are uncorrelated: at lag d a fraction d/segLen of the
% pairs is spurious, which pulls r_I down by exactly (1 - d/segLen). The
% same factor also removes the triangular bias of the FFT estimator, so one
% correction covers both. Defaults to the whole vector.
    v = double(v(:));
    v = v - mean(v);
    N = numel(v);
    if nargin < 3 || isempty(segLen), segLen = N; end
    M = 2^nextpow2(2*N);
    F = fft(v, M);
    C = real(ifft(F .* conj(F)));
    C = C / C(1);
    L = min(maxLag, N-1);
    lag = 0:L;
    acf = C(1:L+1).' ./ (1 - lag/segLen);
end

function [aMean, aStd, nUsed] = local_avg1d(A, maxLag, thr)
% Average 1D correlation method, eq. (4)-(5): one ACF per row of A, the
% threshold crossing read from each, then averaged over the rows.
    A = double(A);
    [m, n] = size(A);
    A = A - mean(A, 2);                          % per-profile mean removal
    M = 2^nextpow2(2*n);
    F = fft(A, M, 2);
    C = real(ifft(F .* conj(F), [], 2));
    C = C ./ C(:,1);                             % r_I = 1 at zero lag
    L = min(maxLag, n-1);
    lag = 0:L;
    C = C(:, 1:L+1) ./ (1 - lag/n);              % overlap correction

    aMean = nan(1, numel(thr));
    aStd  = nan(1, numel(thr));
    nUsed = zeros(1, numel(thr));
    for t = 1:numel(thr)
        a = nan(m, 1);
        for i = 1:m
            a(i) = local_cross(lag, C(i,:), thr(t));
        end
        ok = ~isnan(a);
        nUsed(t) = sum(ok);
        aMean(t) = mean(a(ok));
        aStd(t)  = std(a(ok));
    end
end

function w = local_widths(lag, prof, thr)
% Threshold crossings of an ACF profile, positive lags only.
    pos  = lag >= 0;
    lagp = lag(pos);  p = prof(pos);
    w = nan(1, numel(thr));
    for t = 1:numel(thr)
        w(t) = local_cross(lagp, p, thr(t));
    end
end

function a = local_cross(lag, p, level)
% First down-crossing of LEVEL, linearly interpolated for sub-pixel accuracy.
    idx = find(p < level, 1, 'first');
    if isempty(idx) || idx == 1
        a = NaN;  return
    end
    p1 = p(idx-1);  p2 = p(idx);
    a  = lag(idx-1) + (p1 - level) / (p1 - p2) * (lag(idx) - lag(idx-1));
end

function B = local_gaussblur(A, sigma)
% Separable Gaussian low pass, evaluated in the frequency domain so that a
% wide kernel costs no more than a narrow one. The image is first padded by
% replication over 3*sigma, which keeps the circular convolution from
% wrapping the opposite edge into the estimate.
    r = ceil(3*sigma);
    P = local_padreplicate(double(A), r);
    [m, n] = size(P);

    gy = exp(-((-floor(m/2):ceil(m/2)-1).^2) / (2*sigma^2)).';
    gx = exp(-((-floor(n/2):ceil(n/2)-1).^2) / (2*sigma^2));
    gy = ifftshift(gy / sum(gy));
    gx = ifftshift(gx / sum(gx));

    P = real(ifft2(fft2(P) .* (fft(gy) * fft(gx))));
    B = single(P(r+1:end-r, r+1:end-r));
end

function P = local_padreplicate(A, r)
% Replicate padding without the Image Processing Toolbox.
    [m, n] = size(A);
    ri = [ones(1,r), 1:m, m*ones(1,r)];
    ci = [ones(1,r), 1:n, n*ones(1,r)];
    P  = A(ri, ci);
end
