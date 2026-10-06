function [u0, v0, nWin, window, sel] = bos_background_offset(X, Y, U, V, window, mode)
%BOS_BACKGROUND_OFFSET  Uniform background displacement vector from a window.
%
%   [U0, V0] = BOS_BACKGROUND_OFFSET(X, Y, U, V)
%   [U0, V0, NWIN, WINDOW, SEL] = BOS_BACKGROUND_OFFSET(X, Y, U, V, WINDOW, MODE)
%
%   Estimates the near-uniform displacement VECTOR that a small shift
%   between the reference and measurement images adds to every vector, as
%   the median (or mean) of U and V inside a rectangle of quiescent air.
%   U and V are treated separately: subtracting one scalar from |d| would
%   corrupt the direction of every vector.
%
%   Subtract the result yourself:  U - U0,  V - V0.
%
%   X, Y   : coordinates of the vectors. Either arrays the same size as U,
%            or a row x [1 x Nx] and a column y [Ny x 1] for gridded
%            U, V [Ny x Nx]. Same length unit as WINDOW.
%   U, V   : displacement components. NaN marks invalid vectors.
%   WINDOW : [x1 y1; x2 y2], two opposite corners in either order. Empty
%            or omitted gives the default window below the jet in the 3 cm
%            exports, [3.244 14.9421; -19.27 4.536] mm. The window must lie
%            entirely outside the jet: check it whenever the field of view,
%            the nozzle position or the jet direction changes.
%   MODE   : 'median' (default) | 'mean'
%
%   U0, V0 : the offset, in the units of U and V
%   NWIN   : number of valid vectors inside the window
%   WINDOW : the window actually used
%   SEL    : logical mask of the vectors used, the size of U
%
%   Used by BOS_PIPELINE (step 1) and PLOT_DISPLACEMENT_HISTOGRAM, so both
%   remove the same background.

    if nargin < 5 || isempty(window), window = [3.244 14.9421; -19.27 4.536]; end
    if nargin < 6 || isempty(mode),   mode = 'median'; end
    if ~isequal(size(window), [2 2])
        error('bos_background_offset:badWindow', ...
            'WINDOW must be [x1 y1; x2 y2], got a %dx%d array.', size(window,1), size(window,2));
    end

    xw = sort(window(:,1));  yw = sort(window(:,2));
    inX = X >= xw(1) & X <= xw(2);
    inY = Y >= yw(1) & Y <= yw(2);
    sel = (inX & inY) & ~isnan(U) & ~isnan(V);     % implicit expansion for a grid
    nWin = nnz(sel);
    if nWin == 0
        error('bos_background_offset:emptyWindow', ...
            ['No valid vectors inside the background window x [%.3f, %.3f], ' ...
             'y [%.3f, %.3f]. The data covers x [%.3f, %.3f], y [%.3f, %.3f]. ' ...
             'Choose a window of quiescent air inside the mask.'], ...
             xw(1), xw(2), yw(1), yw(2), min(X(:)), max(X(:)), min(Y(:)), max(Y(:)));
    end

    switch lower(mode)
        case 'median', u0 = median(U(sel));  v0 = median(V(sel));
        case 'mean',   u0 = mean(U(sel));    v0 = mean(V(sel));
        otherwise
            error('bos_background_offset:badMode', ...
                'MODE must be ''median'' or ''mean'', got ''%s''.', mode);
    end
end
