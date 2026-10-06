function sort_results_file(f)
%SORT_RESULTS_FILE  Order the rows of a results CSV by their point number.
%
% The batch scripts append reprocessed points to the end of the file, so a
% run that redoes a single point leaves it out of order. This restores the
% ordering. The first column must hold the point number and the first line
% must be the header.
%
% Bruno Hegedus, MSc thesis - speckle size investigation.

    if ~isfile(f), return, end
    txt = strsplit(fileread(f), newline);
    txt = txt(~cellfun(@isempty, strtrim(txt)));
    if numel(txt) < 3, return, end

    body = txt(2:end);
    p = nan(1, numel(body));
    for i = 1:numel(body)
        v = sscanf(body{i}, '%d', 1);
        if ~isempty(v), p(i) = v; end
    end
    [~, ord] = sort(p);

    fid = fopen(f, 'w');
    fprintf(fid, '%s\n', strtrim(txt{1}));
    for i = ord, fprintf(fid, '%s\n', strtrim(body{i})); end
    fclose(fid);
end
