% Centreline density normalised by the fully expanded density, stacked in the
% style of Panda & Seasholtz (1999), Fig. 4.
here = fileparts(mfilename('fullpath'));
load(fullfile(here, 'centreline.mat'), 'out');

gamma   = 1.4;
rho_amb = 1.2;          % same value as run_bos_density.m
Djet    = 3.5;          % nozzle diameter [mm]
Mj      = [1.082 1.177 1.435];      % Table 6.9, datasets a, b, e
NPR     = [2.088 2.353 3.344];
tag     = {'a', 'b', 'e'};
rho_j   = rho_amb * NPR.^((gamma-1)/gamma);
offset  = 0.8;          % vertical shift between curves (eight minor divisions of 0.1)
xMax    = 10;

fig = figure('Color', 'w', 'Position', [100 100 760 820]);
ax = axes(fig); hold(ax, 'on');
nC = numel(out);
for i = 1:nC
    k   = nC - i;                           % highest Mj at the bottom, as in the reference
    x   = out(i).s / Djet;
    y   = out(i).rho_axis / rho_j(i) + k*offset;
    keep = x >= 0 & x <= xMax & ~isnan(y);
    x = x(keep); y = y(keep);
    [x, ord] = sort(x); y = y(ord);
    plot(ax, x, y, 'k-', 'LineWidth', 1.0);
    mk = 1:6:numel(x);
    plot(ax, x(mk), y(mk), 'k.', 'MarkerSize', 7);
    yline(ax, 1 + k*offset, 'k--', 'LineWidth', 0.9);
    text(ax, xMax - 0.1, 1 + k*offset + 0.42, ...
        sprintf('%s)  M_j = %.3f,  NPR = %.3f', tag{i}, Mj(i), NPR(i)), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 1);
end
xlim(ax, [0 xMax]);
ylim(ax, [0.3 1 + (nC-1)*offset + 0.6]);
ax.YTick = 0.4:0.2:1.4;   % labels apply to the bottom curve only
ax.YAxis.MinorTickValues = 0.3:0.1:ax.YLim(2);
ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
ax.Box = 'on'; ax.FontSize = 12; ax.TickDir = 'in';
xlabel(ax, 'Downstream distance, x/D');
ylabel(ax, 'Time-averaged centreline density, \rho/\rho_j');
title(ax, {'SBOS centreline density, datasets a, b, e', ...
    sprintf('Each curve shifted by %.1f; dashed lines mark \\rho/\\rho_j = 1', offset)}, ...
    'FontWeight', 'normal', 'FontSize', 11);
exportgraphics(fig, fullfile(here, 'centreline_density_panda_style.png'), 'Resolution', 200);

% print a few numbers for checking
for i = 1:nC
    x = out(i).s / Djet; r = out(i).rho_axis / rho_j(i);
    m = x >= 0.5 & x <= 6 & ~isnan(r);
    fprintf('%s: rho_j = %.3f, mean rho/rho_j (x/D 0.5..6) = %.3f, min %.3f, max %.3f, far field (x/D>9) = %.3f (ambient = %.3f)\n', ...
        tag{i}, rho_j(i), mean(r(m)), min(r(m)), max(r(m)), mean(r(x > 9), 'omitnan'), rho_amb/rho_j(i));
end

%% ---- Dataset e against Panda & Seasholtz (1999), Fig. 4, M_j = 1.43 ----
% 'Pandas data.csv' holds the digitised M_j = 1.43 curve, already unshifted:
% column 1 is x/D, column 2 is rho/rho_j.
P = readmatrix(fullfile(here, 'Pandas data.csv'), 'NumHeaderLines', 1);
[~, iu] = unique(P(:,1), 'stable');      % the digitised file repeats x/D = 7
P = P(iu, :);

iE = find(strcmp(tag, 'e'));
xE = out(iE).s / Djet;
yE = out(iE).rho_axis / rho_j(iE);
keep = xE >= 0 & xE <= xMax & ~isnan(yE);
[xE, ord] = sort(xE(keep)); yE = yE(keep); yE = yE(ord);

fig2 = figure('Color', 'w', 'Position', [120 120 900 420]);
ax2 = axes(fig2); hold(ax2, 'on');
plot(ax2, xE, yE, 'k-', 'LineWidth', 1.0);
plot(ax2, P(:,1), P(:,2), 'o--', 'Color', [0.80 0.15 0.10], 'MarkerSize', 5, ...
     'MarkerFaceColor', [0.80 0.15 0.10], 'LineWidth', 1.0);
yline(ax2, 1, 'k:', 'LineWidth', 0.9);
xlim(ax2, [0 xMax]); ylim(ax2, [0.3 1.6]);
ax2.Box = 'on'; ax2.FontSize = 12; ax2.XMinorTick = 'on'; ax2.YMinorTick = 'on';
xlabel(ax2, 'Downstream distance, x/D');
ylabel(ax2, 'Centreline density, \rho/\rho_j');
legend(ax2, {sprintf('SBOS, dataset e (M_j = %.3f)', Mj(iE)), ...
             'Panda & Seasholtz (1999), M_j = 1.43', '\rho/\rho_j = 1'}, ...
       'Location', 'northeast');
exportgraphics(fig2, fullfile(here, 'centreline_density_e_vs_panda.png'), 'Resolution', 200);
