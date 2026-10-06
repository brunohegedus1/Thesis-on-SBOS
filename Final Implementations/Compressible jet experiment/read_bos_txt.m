function D = read_bos_txt(filename, opts)
%READ_BOS_TXT  Read a semicolon-delimited BOS/PIV displacement export onto a grid.
%
%   D = READ_BOS_TXT(FILENAME)
%   D = READ_BOS_TXT(FILENAME, OPTS)
%
%   Reads a text export with a one-line header and semicolon-separated
%   columns of the form
%
%     x [mm];y [mm];x-displacement [mm];y-displacement [mm];Displacement [mm];...
%
%   and reshapes the scattered rows onto the regular (x,y) grid they
%   describe. Points whose x- AND y-displacement are both exactly zero
%   are treated as MASKED/INVALID (not as genuine zero displacement) and
%   are returned as NaN -- in these exports the zero-filled border is
%   padding outside the interrogated region, not measured data.
%
%   OPTS (struct, all fields optional)
%     opts.delimiter  : field delimiter, default ';'
%     opts.headerLines: header lines to skip, default 1
%     opts.colX       : column index of x,  default 1
%     opts.colY       : column index of y,  default 2
%     opts.colU       : column index of x-displacement, default 3
%     opts.colV       : column index of y-displacement, default 4
%     opts.maskZeros  : treat (u==0 & v==0) as invalid, default true
%
%   OUTPUT struct D with fields
%     D.x  : [1 x Nx] unique x coordinates (ascending)
%     D.y  : [Ny x 1] unique y coordinates (ascending)
%     D.U  : [Ny x Nx] x-displacement, NaN where invalid
%     D.V  : [Ny x Nx] y-displacement, NaN where invalid
%     D.mag: [Ny x Nx] displacement magnitude, NaN where invalid
%
%   The displacement units are whatever the file uses (mm in this export);
%   all downstream functions are unit-agnostic as long as the coordinate
%   and displacement units are used consistently.

    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = local_default(opts, 'delimiter',   ';');
    opts = local_default(opts, 'headerLines', 1);
    opts = local_default(opts, 'colX', 1);
    opts = local_default(opts, 'colY', 2);
    opts = local_default(opts, 'colU', 3);
    opts = local_default(opts, 'colV', 4);
    opts = local_default(opts, 'maskZeros', true);

    if ~isfile(filename)
        error('read_bos_txt:notFound', 'File not found: %s', filename);
    end

    % readmatrix (R2019a+) handles CRLF and the text header transparently;
    % fall back to textscan on older MATLAB and on Octave.
    if exist('readmatrix','file') == 2
        raw = readmatrix(filename, 'Delimiter', opts.delimiter, ...
                         'NumHeaderLines', opts.headerLines);
    else
        raw = local_read_fallback(filename, opts.delimiter, opts.headerLines);
    end

    nCols = size(raw,2);
    need  = max([opts.colX opts.colY opts.colU opts.colV]);
    if nCols < need
        error('read_bos_txt:tooFewColumns', ...
            'File has %d columns but column %d was requested. Check opts.col* and the delimiter.', ...
            nCols, need);
    end

    xr = raw(:,opts.colX);  yr = raw(:,opts.colY);
    ur = raw(:,opts.colU);  vr = raw(:,opts.colV);

    keep = ~(isnan(xr) | isnan(yr));    % drop any trailing/blank rows
    xr=xr(keep); yr=yr(keep); ur=ur(keep); vr=vr(keep);

    % Unique coordinates define the grid. Round to a tolerance derived from
    % the data itself so floating-point noise in the export does not create
    % spurious extra grid lines.
    xs = sort(unique(xr));  ys = sort(unique(yr));
    tolx = 0.01*median(diff(xs));  toly = 0.01*median(diff(ys));
    xu = local_uniquetol(xs, tolx);
    yu = local_uniquetol(ys, toly);
    Nx = numel(xu);  Ny = numel(yu);

    % Map every row to its grid cell. Done with a nearest-neighbour lookup
    % rather than a pairwise distance matrix: for this dataset the latter
    % would be 941 x 701986 elements (several GB) and exhaust memory.
    ix = interp1(xu, 1:Nx, xr, 'nearest', 'extrap');
    iy = interp1(yu, 1:Ny, yr, 'nearest', 'extrap');
    ix = min(max(round(ix), 1), Nx);
    iy = min(max(round(iy), 1), Ny);

    U = nan(Ny, Nx);  V = nan(Ny, Nx);
    U(sub2ind([Ny Nx], iy(:), ix(:))) = ur;
    V(sub2ind([Ny Nx], iy(:), ix(:))) = vr;

    if opts.maskZeros
        bad = (U==0 & V==0);
        U(bad) = NaN;  V(bad) = NaN;
    end
    
    D.x = xu(:).';
    D.y = yu(:);
    D.U = U;
    D.V = V;
    D.mag = hypot(U, V);

    fprintf('read_bos_txt: %d x %d grid (%d rows), %.1f%% valid\n', ...
        Ny, Nx, numel(xr), 100*mean(~isnan(U(:))));
end

% ======================================================================
function raw = local_read_fallback(filename, delim, nHeader)
%Portable numeric reader for older MATLAB / Octave (no readmatrix).
    fid = fopen(filename, 'r');
    if fid < 0
        error('read_bos_txt:cannotOpen', 'Could not open %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));
    for k = 1:nHeader
        fgetl(fid);
    end
    % probe the first data line to count columns
    pos  = ftell(fid);
    line = fgetl(fid);
    if ~ischar(line)
        error('read_bos_txt:noData', 'No data rows found in %s', filename);
    end
    nCol = numel(strsplit(strtrim(line), delim));
    fseek(fid, pos, 'bof');

    fmt = repmat('%f', 1, nCol);
    C   = textscan(fid, fmt, 'Delimiter', delim, ...
                   'CollectOutput', true, 'EndOfLine', '\r\n');
    raw = C{1};
    if isempty(raw)
        fseek(fid, pos, 'bof');
        C   = textscan(fid, fmt, 'Delimiter', delim, 'CollectOutput', true);
        raw = C{1};
    end
end

% ======================================================================
function s = local_default(s, f, v)
    if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

% ======================================================================
function u = local_uniquetol(v, tol)
    v = sort(v(:));
    u = v(1);
    for k = 2:numel(v)
        if v(k) - u(end) > tol
            u(end+1,1) = v(k); %#ok<AGROW>
        end
    end
end
