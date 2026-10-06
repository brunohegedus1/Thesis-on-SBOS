function [I, pitch_um] = read_intensity_map(f)
%READ_INTENSITY_MAP  Read a DaVis intensity export into an image matrix.
%
%   [I, pitch_um] = read_intensity_map(f)
%
% Handles both export formats produced so far:
%
%   x [pixel];y [pixel];Average [counts];Number of data points [samples]
%   x [mm];y [mm];Average [];Number of data points [samples]
%
% The two differ only in the coordinate columns, so the grid is recovered
% from the data rather than from the units: x varies fastest, therefore the
% row length is the number of samples before y first changes. PITCH_UM is
% the sample spacing in micrometres when the file carries millimetres, and
% NaN when it carries pixel indices.
%
% I is returned as single, with rows along y and columns along x. The
% vertical direction follows the order of the file, which for a descending
% y column means the first row is the top of the frame. A flip does not
% affect an autocorrelation, so no attempt is made to normalize it.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

    fid = fopen(f, 'r');
    assert(fid > 0, 'Cannot open %s', f);
    header = fgetl(fid);
    head = textscan(fid, '%f%f%*f%*[^\n]', 20000, 'Delimiter', ';');
    fclose(fid);

    x = head{1};  y = head{2};
    assert(numel(x) > 2, 'File %s holds no data rows.', f);

    nx = find(y ~= y(1), 1, 'first') - 1;
    assert(~isempty(nx) && nx > 1, ...
        ['Could not determine the row length of %s. The first %d rows all ' ...
         'share one y value, so either the grid is wider than that or the ' ...
         'file is not ordered with x varying fastest.'], f, numel(y));

    % The coordinates are written with six significant figures, so the
    % spacing wobbles in the last digit. The tolerance below accepts that
    % rounding while still rejecting a genuinely irregular or unsorted grid.
    dx = diff(x(1:nx));
    assert(all(dx > 0) && (max(dx) - min(dx)) < 0.05 * median(dx), ...
        'The x coordinates of %s are not a regular increasing grid.', f);

    if contains(lower(header), '[mm]')
        pitch_um = median(dx) * 1000;
    else
        pitch_um = NaN;                      % pixel indices, no physical scale
    end

    fid = fopen(f, 'r');
    C = textscan(fid, '%*f%*f%f%*[^\n]', 'Delimiter', ';', 'HeaderLines', 1);
    fclose(fid);
    v = C{1};

    ny = numel(v)/nx;
    assert(ny == round(ny), ...
        ['Grid is incomplete in %s: %d samples do not divide into rows of ' ...
         '%d.'], f, numel(v), nx);
    I = single(reshape(v, nx, ny)).';        % rows = y, columns = x
end
