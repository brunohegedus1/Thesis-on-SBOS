function d = bos_data_dir()
%BOS_DATA_DIR  Folder holding the raw BOS .txt exports.
%
%   Single place to change if the data moves again. Every script that reads
%   an export goes through here:
%
%       D = read_bos_txt(fullfile(bos_data_dir, 'BOS_3cm_8bar_12x12.txt'));
%
%   The exports live outside the code folder because they are ~45 MB each
%   and do not belong in a repository.

    % The .txt exports sit in BOS_data next to this function.
    d = fullfile(fileparts(mfilename('fullpath')), 'BOS_data');

    if ~isfolder(d)
        error('bos_data_dir:missing', ...
            ['BOS data folder not found:\n    %s\n' ...
             'Edit BOS_DATA_DIR if the exports have moved.'], d);
    end
end
