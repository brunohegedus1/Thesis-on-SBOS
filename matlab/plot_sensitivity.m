% PLOT_SENSITIVITY  Measured sensitivity against signed defocus distance.
%
%   Plots S against l for one campaign of the angled glass experiments,
%   with error bars sized from the run-to-run scatter of the measured
%   displacement.
%
%   Input files (relative to this script):
%       ../data/angled_glass_experimental.csv
%           Points, sheet_ref, campaign, f#, f, l, df, S, M, type,
%           delta (px), StdDev (px)
%       ../data/angled_glass_errors.csv
%           Points, sheet_ref, df/S/M/delta errors in %
%
%   Numbering. Points 1-20 follow the tables in Chapter 5. The sheet_ref
%   column keeps the original workbook labels (P6-P14 and P16-P16.10) so a
%   row can still be traced back to the spreadsheet.
%
%   Campaigns. Both used the angled glass plate, but they answer different
%   questions and share f = 200 mm, so they must not be pooled:
%
%     'slider-sweep'   Points 10-20. The in-focus plane is fixed and the
%                      object rides a linear slider, so M is constant at
%                      0.746 and only l changes. Point 10 sits in the
%                      in-focus plane (l = 0), where epsilon and therefore
%                      S are undefined; its S field is empty and it is
%                      dropped before plotting.
%
%     'setup-matrix'   Points 1-9. Each point is a different optical
%                      layout, so M changes from point to point.
%
%   Sign convention. The CSV stores the magnitude of the defocus distance
%   in l and its sign in type (-1 in front of the object, +1 behind it),
%   following Chapter 4. The signed value is rebuilt here so that negative
%   defocus lands on the left half of the plot.
%
%   Grouping. Set groupBy to 'defocus' to separate the branches by the
%   sign of l, which is the useful split when one f-stop is used
%   throughout, or to 'fstop' to colour by f-number.
%
%   Choice of error bar. Two quantities are available and they measure
%   different things:
%
%     'scatter' (default)  StdDev/delta, the repeatability of the
%                          displacement measurement. This is an
%                          uncertainty, which is what an error bar should
%                          show.
%
%     'bias'               The S error column, equal to
%                          |S_exp - l*M_th| / (l*M_th), the deviation from
%                          the theoretical prediction. This is an
%                          agreement metric, not an uncertainty. Plotting
%                          it as a bar double-counts the gap already
%                          visible between the point and theory, and it
%                          understates points that happen to land close to
%                          theory.
%
%   Use 'bias' only to reproduce the original spreadsheet chart.

clear; close all; clc

%% Configuration
scriptDir   = fileparts(mfilename('fullpath'));
dataDir     = fullfile(scriptDir, '..', 'data');
campaign    = 'slider-sweep';   % 'slider-sweep' (points 10-20) or 'setup-matrix' (1-9)
focalLength = 0.200;            % [m] lens to plot; use 0.105 for the short lens
groupBy     = 'defocus';        % 'defocus' or 'fstop'
errorSource = 'scatter';        % 'scatter' or 'bias'
saveFigure  = false;
outFile     = fullfile(scriptDir, ...
              sprintf('sensitivity_%s_f%03.0fmm.pdf', campaign, focalLength*1e3));

%% Load and merge
expData = loadCsv(fullfile(dataDir, 'angled_glass_experimental.csv'));
errData = loadCsv(fullfile(dataDir, 'angled_glass_errors.csv'));

% sheet_ref appears in both files. Drop the duplicate so the join does not
% produce suffixed columns.
errData.sheet_ref = [];

T = innerjoin(expData, errData, 'Keys', 'Points');

% Keep the requested campaign, then the requested lens. The two campaigns
% share f = 200 mm, so filtering on focal length alone would mix them.
T = T(T.campaign == string(campaign), :);
T = T(abs(T.f - focalLength) < 1e-9, :);

% Points sitting in the in-focus plane have no defined S and drop out.
T = T(~isnan(T.S), :);

if isempty(T)
    error('No points for campaign ''%s'' at f = %g m.', campaign, focalLength);
end

%% Derived quantities
lSigned = T.l .* T.type;                                % [m] negative = in front

% Order by position so the console summary reads left to right like the plot.
[lSigned, ord] = sort(lSigned);
T = T(ord, :);

relScatter = T.("StdDev (px)") ./ T.("delta (px)");     % [-] repeatability
relBias    = T.("S error (%)") / 100;                   % [-] deviation from theory

switch lower(errorSource)
    case 'scatter'
        relErr   = relScatter;
        barLabel = 'error bars: scatter of the delta displacement';
    case 'bias'
        relErr   = relBias;
        barLabel = 'error bars: deviation from theory';
    otherwise
        error('errorSource must be ''scatter'' or ''bias''.');
end

sErrAbs = T.S .* relErr;                                % [m] bar half-length

% Console summary so the plotted numbers can be checked against the sheet.
summary = table(T.Points, string(T.sheet_ref), T.("f#"), lSigned, T.S, ...
                100*relScatter, 100*relBias, sErrAbs, ...
    'VariableNames', {'Point','SheetRef','fStop','l_signed_m','S_m', ...
                      'scatter_pct','bias_pct','bar_half_m'});
disp(summary);
fprintf('Campaign: %s   bars: %s\n\n', campaign, errorSource);

%% Series definition
switch lower(groupBy)
    case 'defocus'
        keys    = [-1 1];
        inGroup = @(k) sign(lSigned) == k;
        labels  = {'negative defocus (in front)', 'positive defocus (behind)'};
        colours = {[0.259 0.522 0.957], [0.957 0.427 0.063]};
    case 'fstop'
        keys    = sort(unique(T.("f#")), 'descend')';
        inGroup = @(k) T.("f#") == k;
        labels  = arrayfun(@(k) sprintf('f# = %d', k), keys, 'UniformOutput', false);
        colours = arrayfun(@fstopColour, keys, 'UniformOutput', false);
    otherwise
        error('groupBy must be ''defocus'' or ''fstop''.');
end

%% Axes ranges, chosen per campaign so each fills the frame
switch campaign
    case 'slider-sweep'
        xLim = [-0.06 0.06];  xTick = -0.05:0.01:0.05;
        yLim = [0 0.045];     yTick = 0:0.005:0.045;
    otherwise
        xLim = [-0.20 0.25];  xTick = -0.20:0.05:0.25;
        yLim = [0 0.14];      yTick = 0:0.02:0.14;
end

%% Plot
figure('Color', 'w', 'Position', [100 100 780 440]);
hold on; grid on; box on

xline(0, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');

h = gobjects(numel(keys), 1);
for k = 1:numel(keys)
    m = inGroup(keys(k));
    if ~any(m), continue; end

    h(k) = errorbar(lSigned(m), T.S(m), sErrAbs(m), 'o', ...
                    'LineStyle',       'none', ...
                    'Color',           [0.15 0.15 0.15], ...  % bars
                    'MarkerFaceColor', colours{k}, ...
                    'MarkerEdgeColor', colours{k}, ...
                    'MarkerSize',      6, ...
                    'LineWidth',       0.9, ...
                    'CapSize',         6, ...
                    'DisplayName',     labels{k});
end
h = h(isgraphics(h));

xlabel('l [m]');
ylabel('S [m]');
title({sprintf('Sensitivity for f = %g mm  (%s)', focalLength*1e3, campaign), barLabel});

xlim(xLim);  xticks(xTick);
ylim(yLim);  yticks(yTick);
ytickformat('%.4f');

legend(h, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
set(gca, 'FontSize', 10, 'GridAlpha', 0.15, 'Layer', 'top');
set(get(gca, 'Title'), 'FontSize', 11);

hold off

if saveFigure
    exportgraphics(gcf, outFile, 'ContentType', 'vector');
    fprintf('Saved %s\n', outFile);
end

%% Helpers
function T = loadCsv(file)
% Reads a campaign CSV, keeping the header text verbatim so that columns
% such as "delta (px)" and "S error (%)" can be addressed by their real
% names.
    T = readtable(file, 'VariableNamingRule', 'preserve');
end

function c = fstopColour(fStop)
% Colours chosen to match the original spreadsheet chart.
    switch fStop
        case 16, c = [0.957 0.427 0.063];   % orange
        case 11, c = [0.259 0.522 0.957];   % blue
        case 4,  c = [0.180 0.663 0.314];   % green
        otherwise, c = [0.40 0.40 0.40];    % grey fallback
    end
end
