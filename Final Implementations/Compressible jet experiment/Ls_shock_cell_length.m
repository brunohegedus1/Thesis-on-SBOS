% VISUAL_LS
%
%   Shock cell lengths of the underexpanded air jet, measured from the
%   centreline density gradient, fitted with two relations.
%
%   Emden (1899), as given by van Hinsberg & Rosgen (2014) Eq. 25:
%
%       Ls/Djet = k*sqrt(NPR - NPRc)                     (through the origin)
%
%   Hartmann & Lazarus (1941), their Eq. 27, the same form with an offset:
%
%       Ls/Djet = a*sqrt(NPR - NPRc) + b
%
%   NPRc is the choking pressure ratio, van Hinsberg Eq. 26:
%
%       NPRc = ((gamma+1)/2)^(gamma/(gamma-1))
%
%   Two cell lengths are used, and each relation is fitted to one of them:
%     Lvis    first cell, nozzle to the first density maximum (L_max)
%             -> Hartmann-Lazarus fit
%     Ls_avg  average length of all cells AFTER the first
%             -> Emden fit
%   Emden's own coefficient came from a weighted multi-cell average, and
%   Hartmann & Lazarus gave a = 0.98, b = 0.30 for the first cell (Powell
%   2010, Int J Aeroacoust 9:207-236; Hartmann & Lazarus 1941, pp. 35-36).
%   Both coefficients are fitted to the measurements here. Three literature
%   curves are drawn for reference only: Emden's k = 0.88 for air, Hartmann
%   & Lazarus's first-cell a = 0.98, b = 0.30, and Hartmann & Lazarus's
%   cells after the first, which they gave in Emden's form with k = 1.06.
%
%   pbar is the stagnation pressure p0 in bar absolute, estimated from the
%   Prandtl-Meyer fit, so NPR = p0/p_amb.

clear; clc; close all

%% Measurements ------------------------------------------------------------
Djet  = 3.5;                                    % nozzle exit diameter, mm
gamma = 1.4;                                    % air
p_amb = 1.013;                                  % ambient pressure, bar

%pbar = (3:8).'+ p_amb;                                 % tap setting, bar gauge
%pbar = [ 2.082 2.461 3.067 3.218 3.597 3.74].';
%pbar = [ 2.25 2.613 2.887 3.072 3.387 3.925].';
pbar = [ 2.15 2.353 2.647 3.072 3.387 3.925].';
%Lvis = [2.16 3.20 4.17 4.80 5.00 5.40].';       % visually inspected shock cell length, mm
%Lrho_max = [1.84 3.36 4.05 3.59 4 4.65].';       % density data maxima shock cell length, mm
%Lrho_min = [2.3 2.9 3.31 4.74 5.02 4.97].';       % density data minimum shock cell length, mm
%Lvis = (Lrho_min + Lrho_max)/2;
Lvis = [2.13 3.14 3.9 4.48 4.95 5.0].';          % first cell (L_max), mm

% Average length of all cells after the first, datasets a to f, mm. From
% the centreline d(rho)/ds of run_bos_gradients_unfolded.m: common
% least-squares slope of compression peaks, density minima and density
% maxima against cell index, first cell excluded, 0.4 mm detection window.
% 4 to 6 cells per dataset; changes by at most 0.02 mm for 0.2-0.8 mm windows.
Ls_avg = [1.853 2.657 3.358 3.912 4.520 4.806].';

%% Reference data: Emden, digitised ----------------------------------------
% The CSV is already in the EMDEN COORDINATE: column 1 is sqrt(NPR - NPRc),
% column 2 is Ls/Djet. It therefore overlays the right-hand panel as-is, and
% needs NPR = x^2 + NPRc to be drawn on the left.
%
% These points are PLOTTED ONLY. They take no part in either fit, which
% stays a fit to this rig's own measurements -- mixing digitised literature
% points into it would make the reported coefficients a blend of two
% experiments and the residuals meaningless.
% Emden_data.csv sits next to this script. Set EMDEN_CSV to it to draw the
% points; left empty, the figure shows this rig's measurements only, as in
% the thesis.
EMDEN_CSV = '';   % fullfile(fileparts(mfilename('fullpath')), 'Emden_data.csv')
if isfile(EMDEN_CSV)
    Eref   = readmatrix(EMDEN_CSV);
    xE_ref = Eref(:,1);   yE_ref = Eref(:,2);
else
    xE_ref = [];  yE_ref = [];
    warning('Ls_shock_cell_length:noEmdenCsv', ...
            '%s not found; the reference points are skipped.', EMDEN_CSV);
end
emdCol = [0 0 0];

%% Literature coefficients, drawn for reference -----------------------------
kE_lit  = 0.88;                                 % Emden, air (Powell 2010)
aHL_lit = 0.98;  bHL_lit = 0.30;                % Hartmann & Lazarus (1941), first cell
kHL_lit = 1.06;                                 % Hartmann & Lazarus (1941), cells
                                                % after the first, Emden form

%% Derived -----------------------------------------------------------------
NPRc = ((gamma+1)/2)^(gamma/(gamma-1));
NPR  = (pbar)/p_amb;
x    = sqrt(NPR - NPRc);                        % Emden abscissa
y    = Lvis/Djet;                               % first cell, non-dimensional
yA   = Ls_avg/Djet;                             % cells after the first
Mj   = sqrt((NPR.^((gamma-1)/gamma) - 1)*2/(gamma-1));

%% Fit 1: Emden, one parameter, through the origin, to Ls_avg --------------
kE = sum(yA.*x)/sum(x.^2);
rE = yA - kE*x;

%% Fit 2: Hartmann-Lazarus, two parameters, to the first cell --------------
P   = polyfit(x, y, 1);
aHL = P(1);  bHL = P(2);
rHL = y - (aHL*x + bHL);

%% Report ------------------------------------------------------------------
fprintf('gamma = %.2f, Djet = %.2f mm, NPRc = %.4f\n\n', gamma, Djet, NPRc);
fprintf('%8s %7s %6s | %8s %8s %9s | %9s %9s %8s\n', ...
        'p0[bar]','NPR','Mj','Ls1[mm]','Ls1/D','Hart-Laz','Lsavg[mm]','Lsavg/D','Emden');
for i = 1:numel(pbar)
    fprintf('%8.3f %7.3f %6.2f | %8.2f %8.3f %9.3f | %9.3f %9.3f %8.3f\n', ...
        pbar(i), NPR(i), Mj(i), Lvis(i), y(i), aHL*x(i)+bHL, Ls_avg(i), yA(i), kE*x(i));
end

fprintf('\nEmden, cells after the first   Ls/Djet = %.3f*sqrt(NPR - %.3f)\n', kE, NPRc);
fprintf('  rms residual   %.4f  (%.3f mm)   max %.4f\n', ...
        sqrt(mean(rE.^2)), sqrt(mean(rE.^2))*Djet, max(abs(rE)));
fprintf('Hartmann-Lazarus, first cell   Ls/Djet = %.3f*sqrt(NPR - %.3f) %+.3f\n', aHL, NPRc, bHL);
fprintf('  rms residual   %.4f  (%.3f mm)   max %.4f\n', ...
        sqrt(mean(rHL.^2)), sqrt(mean(rHL.^2))*Djet, max(abs(rHL)));

rHLlit = y - (aHL_lit*x + bHL_lit);
fprintf('Hartmann-Lazarus, literature a = %.2f, b = %+.2f  (first cell, not fitted)\n', ...
        aHL_lit, bHL_lit);
fprintf('  rms residual   %.4f  (%.3f mm)   max %.4f\n', ...
        sqrt(mean(rHLlit.^2)), sqrt(mean(rHLlit.^2))*Djet, max(abs(rHLlit)));

rHLk = yA - kHL_lit*x;
fprintf('Hartmann-Lazarus, Emden form k = %.2f  (cells after the first, not fitted)\n', kHL_lit);
fprintf('  rms residual   %.4f  (%.3f mm)   max %.4f\n', ...
        sqrt(mean(rHLk.^2)), sqrt(mean(rHLk.^2))*Djet, max(abs(rHLk)));

fprintf('\nThe Hartmann-Lazarus fit reaches zero at NPR = %.2f. Below that it is\n', ...
        NPRc + (bHL/aHL)^2);
fprintf('unphysical: it is an empirical fit, valid only over the measured range\n');
fprintf('NPR = %.2f to %.2f.\n', min(NPR), max(NPR));

%% Plot --------------------------------------------------------------------
NPRp = linspace(NPRc, 10.5, 400).';
xp   = sqrt(NPRp - NPRc);

% Marker styles shared by both figures
mk1   = {'o', 'MarkerSize',5, 'MarkerFaceColor',[0.15 0.15 0.15], ...
         'MarkerEdgeColor','k', 'LineStyle','none'};             % first cell
mkAvg = {'s', 'MarkerSize',6, 'MarkerFaceColor','w', ...
         'MarkerEdgeColor','k', 'LineWidth',1.1, 'LineStyle','none'}; % Ls_avg
lgFit = {sprintf('Emden, k = %.2f', kE_lit), ...
         sprintf('Hartmann-Lazarus, k = %.2f, b = %+.2f', aHL_lit, bHL_lit), ...
         sprintf('Hartmann-Lazarus, Emden form, k = %.2f', kHL_lit), ...
         sprintf('Emden fitted to L_{s,avg}, k = %.3f', kE), ...
         sprintf('Hartmann-Lazarus fitted to L_s, k = %.3f, b = %+.3f', aHL, bHL), ...
         'L_s, first cell', ...
         'L_{s,avg}, cells after the first'};

figure('Color','w','Position',[100 100 1050 680]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% --- left: physical coordinate. Both models are curves here ---------------
nexttile([2 1]); hold on
plot(NPRp, kE_lit*xp,    ':',  'LineWidth',1.3, 'Color',[0.45 0.45 0.45]);
plot(NPRp, aHL_lit*xp + bHL_lit, '-.', 'LineWidth',1.3, 'Color',[0.45 0.45 0.45]);
plot(NPRp, kHL_lit*xp,   '--', 'LineWidth',1.3, 'Color',[0.45 0.45 0.45]);
plot(NPRp, kE*xp,        '-',  'LineWidth',1.7);
plot(NPRp, aHL*xp + bHL, '--', 'LineWidth',1.7);
plot(NPR,  y,  mk1{:});
plot(NPR,  yA, mkAvg{:});
if ~isempty(xE_ref)
    plot(xE_ref.^2 + NPRc, yE_ref, '*', 'MarkerSize',7, ...
         'LineWidth',1.1, 'Color',emdCol, 'LineStyle','none');
end
yline(0, 'k-', 'HandleVisibility','off');
xline(NPRc, 'k:', 'NPR_c', 'LabelVerticalAlignment','bottom', ...
      'HandleVisibility','off');
xlabel('NPR = p_0 / p_{amb}'); ylabel('L_s / d_N');
lg = lgFit;
if ~isempty(xE_ref)
    lg{end+1} = sprintf('Emden data, digitised (%d pts)', numel(xE_ref));
end
legend(lg, 'Location','northwest','Box','off','FontSize',8);
grid on; xlim([1.5 10.5]); ylim([-0.3 2.6]);
title('Against NPR');

% --- right: Emden coordinate. Both models become straight lines, so the
%     only difference between them is whether the line passes through the
%     origin. The offset in Hartmann-Lazarus is the intercept -------------
nexttile([2 1]); hold on
plot(xp, kE_lit*xp,    ':',  'LineWidth',1.3, 'Color',[0.45 0.45 0.45]);
plot(xp, aHL_lit*xp + bHL_lit, '-.', 'LineWidth',1.3, 'Color',[0.45 0.45 0.45]);
plot(xp, kHL_lit*xp,   '--', 'LineWidth',1.3, 'Color',[0.45 0.45 0.45]);
plot(xp, kE*xp,        '-',  'LineWidth',1.7);
plot(xp, aHL*xp + bHL, '--', 'LineWidth',1.7);
plot(x,  y,  mk1{:});
plot(x,  yA, mkAvg{:});
if ~isempty(xE_ref)
    plot(xE_ref, yE_ref, '*', 'MarkerSize',7, ...
         'LineWidth',1.1, 'Color',emdCol, 'LineStyle','none');
end
yline(0, 'k-'); plot(0, bHL, 'd', 'MarkerSize',7, ...
     'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k');
text(0.06, bHL, sprintf('  b = %+.3f', bHL), 'FontSize',8, ...
     'VerticalAlignment','middle');
xlabel('(NPR - NPR_c)^{1/2}'); ylabel('L_s / d_N');
grid on; xlim([0 3]); ylim([-0.6 2.6]);
title('Emden coordinate');

nexttile([1 2]); hold on
plot(NPR, rE*Djet,  'o-',  'LineWidth',1.2, 'MarkerFaceColor','auto');
plot(NPR, rHL*Djet, 's--', 'LineWidth',1.2, 'MarkerFaceColor','auto');
yline(0,'k:');
xlabel('NPR'); ylabel('residual [mm]'); grid on; xlim([1.5 10.5]);
legend({'Emden, to L_{s,avg}','Hartmann-Lazarus, to L_s'}, 'Location','best', ...
       'Box','off', 'Orientation','horizontal');
title('Fit residuals','FontWeight','normal');
sgtitle('Shock cell length of an underexpanded air jet, D_{jet} = 3.5 mm','FontWeight','bold');

exportgraphics(gcf, 'visual_Ls.png', 'Resolution', 150);
fprintf('\nfigure written to visual_Ls.png\n');

%% Plot 2: Emden coordinate only --------------------------------------------
% The same relations and data as the right-hand panel above, on their own.
% Drawn after the export above, so visual_Ls.png still captures figure 1.
% Colours are set explicitly so they do not depend on MATLAB's colour order.
colEfit  = [0.93 0.69 0.13];                    % Emden, fitted
colHLfit = [0.49 0.18 0.56];                    % Hartmann-Lazarus, fitted
colLit   = [0.45 0.45 0.45];                    % literature curves

% This figure is the one used in the report. It is sized to the text width
% of the class (a4paper with hscale = 0.75, so 0.75*210 = 157.5 mm) and all
% of its text is set to 8 pt below, so that the exported image, included at
% width=\textwidth, shows 8 pt letters on the page.
FIG_W = 15.75;                                  % cm, the text width
FIG_H = 10.50;                                  % cm. The former aspect was
                                                % much wider, but at 8 pt the
                                                % legend then covered the data
figure('Color','w','Units','centimeters', ...
       'Position',[4 4 FIG_W FIG_H], 'PaperPositionMode','auto');
hold on; grid on; box on;
plot(xp, kE_lit*xp,              ':',  'LineWidth',1.3, 'Color',colLit);
plot(xp, aHL_lit*xp + bHL_lit,   '-.', 'LineWidth',1.3, 'Color',colLit);
plot(xp, kHL_lit*xp,             '--', 'LineWidth',1.3, 'Color',colLit);
plot(xp, kE*xp,                  '-',  'LineWidth',1.7, 'Color',colEfit);
plot(xp, aHL*xp + bHL,           '--', 'LineWidth',1.7, 'Color',colHLfit);
plot(x,  y,  mk1{:});
plot(x,  yA, mkAvg{:});
lg2 = lgFit;
if ~isempty(xE_ref)
    plot(xE_ref, yE_ref, '*', 'MarkerSize',7, 'LineWidth',1.1, ...
         'Color',emdCol, 'LineStyle','none');
    lg2{end+1} = sprintf('Emden data, digitised (%d pts)', numel(xE_ref));
end
xlabel('(NPR - NPR_c)^{1/2}'); ylabel('L_s / d_N');
xlim([0 3]); ylim([0 2.6]);
% Below the axes in two columns: at 8 pt the longest entry is wider than
% half the panel, so inside the axes it would sit over the measurements.
legend(lg2, 'Location','southoutside', 'NumColumns',2, 'Box','off', 'FontSize',8);
title('Shock cell length relationships');

% Every letter at 8 pt: tick labels, axis labels and title included. MATLAB
% draws the axis labels and the title larger than the axes font size by
% default, so both multipliers are set to 1 first.
set(gca, 'FontSize',8, 'LabelFontSizeMultiplier',1, 'TitleFontSizeMultiplier',1);
set(findall(gcf,'-property','FontSize'), 'FontSize',8);

% Printed, not exported with EXPORTGRAPHICS: the latter crops the margins,
% which would shrink the image and enlarge the text once LaTeX scales it
% back up to the text width.
print(gcf, 'shock cell lengths.jpg', '-djpeg90', '-r300');
fprintf('figure written to shock cell lengths.jpg (%.2f x %.2f cm, 8 pt text)\n', ...
        FIG_W, FIG_H);
