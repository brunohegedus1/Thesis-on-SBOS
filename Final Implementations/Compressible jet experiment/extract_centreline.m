% Runs run_bos_density.m (unchanged logic, plots stripped) for several files
% and saves the centreline density of each.
src  = fileparts(mfilename('fullpath'));           % RUN_BOS_DENSITY sits here
if ~isfile(fullfile(src, 'run_bos_density.m'))
    src = 'C:\Users\bruno\Documents\MSC THESIS\BOS jet Post processing';
end
here = fileparts(mfilename('fullpath'));
addpath(src);
txt = fileread(fullfile(src, 'run_bos_density.m'));
iPlot = strfind(txt, '%% ======================= 9. PLOTS');
iHelp = strfind(txt, '%% ======================= local helpers');
body  = [txt(1:iPlot-1) txt(iHelp:end)];
body  = strrep(body, 'clear; clc; %close all;', '');

files = {'BOS_3cm_12x120001.csv', 'BOS_3cm_12x120001_1.csv', 'BOS_3cm_12x120001_4.csv'};
out = struct();
for i = 1:numel(files)
    b = regexprep(body, '^DATAFILE = ''[^'']*'';', ...
                  sprintf('DATAFILE = ''%s'';', files{i}), 'once', 'lineanchors');
    tmp = fullfile(here, sprintf('tmp_run_%d.m', i));
    fid = fopen(tmp, 'w'); fwrite(fid, b); fclose(fid);
    run(tmp);
    delete(tmp);
    out(i).file = files{i};
    out(i).s = R.s;
    out(i).rho_axis = rho(1,:);
    out(i).rho_r025 = mean(rho(R.nHalf <= 0.25, :), 1);
    out(i).lambda = abelInfo.lambdaRel;
    clearvars -except src here body files out i
end
save(fullfile(here, 'centreline.mat'), 'out');
T = [];
for i = 1:numel(out)
    fprintf('%s: s %.2f..%.2f, rho axis %.3f..%.3f\n', out(i).file, ...
        min(out(i).s), max(out(i).s), min(out(i).rho_axis), max(out(i).rho_axis));
end
