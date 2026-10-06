% PLOT_SENSITIVITY  Measured sensitivity against signed defocus distance.
%
%   Plots S against l for one campaign of the angled glass experiments,
%   with error bars sized from the run-to-run scatter of the measured
%   displacement.
%
%   Input files (relative to this script):
%       data/angled_glass_experimental.csv
%           Points, sheet_ref, campaign, f#, f, l, df, S, M, type,
%           delta (px), StdDev (px)
%       data/angled_glass_errors.csv
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
dataDir     = fullfile(scriptDir, 'data');          % the CSVs sit in this folder
campaign    = 'setup-matrix';   % 'slider-sweep' (points 10-20) or 'setup-matrix' (1-9)
focalLength = 0.200;            % [m] lens to plot; use 0.105 for the short lens
groupBy     = 'defocus';        % 'defocus' or 'fstop'
errorSource = 'scatter';        % 'scatter' or 'bias'
saveFigure  = true;
outFormat   = 'jpg';            % 'jpg' (raster) or 'pdf' (vector)

% Text size on the page. The figure is drawn at the size it occupies in the
% report and every label is set to fontSize, so the exported file needs no
% rescaling by LaTeX and its letters print at exactly fontSize points.
% textWidth is the class text width (a4paper, hscale = 0.75, so
% 0.75*210 = 157.5 mm) and widthFrac the subfigure width used in Chapter 6,
% where both panels sit in 0.49\textwidth boxes.
fontSize    = 8;                % pt, all text in the figure
textWidth   = 15.75;            % cm, \textwidth of the report class
widthFrac   = 0.49;             % fraction of \textwidth the panel occupies
figW        = widthFrac*textWidth;
figH        = 6.00;             % cm, free to choose, the panels sit side by side
printDpi    = 600;              % raster resolution of the exported jpg

% File names as they are included in the report.
switch campaign
    case 'slider-sweep', baseName = 'S vs l slider experiment';
    otherwise,           baseName = 'S vs l experimental';
end
outFile     = fullfile(scriptDir, sprintf('%s.%s', baseName, outFormat));

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
    % Tick spacing is wider than the data would allow, because at 8 pt on a
    % panel 7.7 cm wide the former labels collided: eleven labels left about
    % 7 mm each for text roughly 8.5 mm long.
    case 'slider-sweep'
        xLim = [-0.06 0.06];  xTick = -0.05:0.025:0.05;
        yLim = [0 0.045];     yTick = 0:0.01:0.04;
    otherwise
        xLim = [-0.20 0.25];  xTick = -0.20:0.10:0.20;
        yLim = [0 0.14];      yTick = 0:0.02:0.14;
end

%% Plot
figure('Color', 'w', 'Units', 'centimeters', ...
       'Position', [4 4 figW figH], 'PaperPositionMode', 'auto');
hold on; grid on; box on

xline(0, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');

h = gobjects(numel(keys), 1);
for k = 1:numel(keys)
    m = inGroup(keys(k));
    if ~any(m), continue; end

    % Marker, bar and cap sizes are in points, so they do not shrink with
    % the figure. They are set for the printed panel size, where the former
    % values, chosen for a figure about three times wider, covered the data.
    h(k) = errorbar(lSigned(m), T.S(m), sErrAbs(m), 'o', ...
                    'LineStyle',       'none', ...
                    'Color',           [0.15 0.15 0.15], ...  % bars
                    'MarkerFaceColor', colours{k}, ...
                    'MarkerEdgeColor', colours{k}, ...
                    'MarkerSize',      2.8, ...
                    'LineWidth',       0.8, ...
                    'CapSize',         3, ...
                    'DisplayName',     labels{k});
end
h = h(isgraphics(h));

xlabel('l [m]');
ylabel('S [m]');
% Titles as they appear in the thesis figures. barLabel is still printed
% to the console so the error bar choice stays traceable.
switch campaign
    case 'slider-sweep'
        % Shorter than the former "(linear slider sweep)": at 8 pt on a
        % panel 7.7 cm wide the longer title ran past the axes.
        figTitle = sprintf('Sensitivity for f= %gmm, slider sweep', focalLength*1e3);
    otherwise
        figTitle = sprintf('Sensitivity comparison for f= %gmm', focalLength*1e3);
end
title(figTitle);
fprintf('%s\n', barLabel);

xlim(xLim);  xticks(xTick);  xtickangle(0);   % keep labels horizontal
ylim(yLim);  yticks(yTick);
ytickformat('%.4f');

% The defocus labels are long, so they are stacked; the f-number labels are
% short enough to share one row. Side by side, the defocus pair is about
% 9 cm wide and would run off a panel of 7.7 cm.
switch lower(groupBy)
    case 'defocus', legCols = 1;
    otherwise,      legCols = numel(h);
end
legend(h, 'Location', 'southoutside', 'NumColumns', legCols, 'Box', 'off');
% One size for every letter in the figure. MATLAB draws axis labels and the
% title larger than the axes font size by default, so both multipliers are
% set to 1 before the font size is applied to all text objects.
set(gca, 'FontSize', fontSize, 'GridAlpha', 0.15, 'Layer', 'top', ...
         'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
ax = gca;
ax.Toolbar.Visible = 'off';        % keep the axes toolbar out of exported images
set(get(gca, 'Title'), 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);

hold off

if saveFigure
    % PRINT, not EXPORTGRAPHICS: the latter crops the margins, which makes
    % the file narrower than the panel. LaTeX would then stretch it back to
    % 0.49\textwidth and the letters would print larger than fontSize.
    if strcmpi(outFormat, 'pdf')
        set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [figW figH], ...
                 'PaperPosition', [0 0 figW figH]);
        print(gcf, outFile, '-dpdf', '-vector');
    else
        print(gcf, outFile, '-djpeg95', sprintf('-r%d', printDpi));
    end
    fprintf('Saved %s  (%.2f x %.2f cm, %d pt text)\n', outFile, figW, figH, fontSize);
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
