%SURFMAP_SHEET  One printed sheet of surface maps: rows = f-stops, columns = quantities.
function surfmap_sheet(panels, l, f, fnum, S, smin_raw, isoS, isoSmin, W, axh, FS, outfile, RES)
% One sheet: rows = f-stops, columns = quantities, one colourbar per column.
nc = numel(panels); nr = numel(fnum);

left = 1.65; right = 0.35; top = 0.80; bottom = 2.85;
colgap = 0.95; rowgap = 0.62;
axw = (W - left - right - (nc-1)*colgap) / nc;
H   = bottom + nr*axh + (nr-1)*rowgap + top;

fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w', ...
             'PaperUnits', 'centimeters', 'PaperSize', [W H], 'Visible', 'off');

for c = 1:nc
    x0 = left + (c-1)*(axw + colgap);
    for r = 1:nr
        y0 = bottom + (nr-r)*(axh + rowgap);
        ax = axes(fig, 'Units', 'centimeters', 'Position', [x0 y0 axw axh], ...
                  'FontSize', FS, 'Layer', 'top', 'Box', 'on', 'TickDir', 'out');
        hold(ax, 'on');
        % one black line per fill level, as in the original maps
        contourf(ax, l, f, panels(c).data(:,:,r).', panels(c).nlev, ...
                 'LineColor', 'k', 'LineWidth', 0.2);
        contour(ax, l, f, S(:,:,r).', [isoS isoS], 'w', 'LineWidth', 1.2);
        contour(ax, l, f, smin_raw(:,:,r).', [isoSmin isoSmin], 'r', 'LineWidth', 1.2);
        clim(ax, panels(c).clim);
        xlim(ax, [min(l) max(l)]); ylim(ax, [min(f) max(f)]);
        xticks(ax, 0:0.05:0.15); yticks(ax, 0.05:0.05:0.2);
        if r < nr, xticklabels(ax, []); else, xlabel(ax, '$|l|$ [m]', 'FontSize', FS, 'Interpreter', 'latex'); end
        if c > 1
            yticklabels(ax, []);
        else
            % the row header rides on the y label, so no extra left margin is needed
            ylabel(ax, {sprintf('$f_{\\#} = %d$', fnum(r)), '$f$ [m]'}, ...
                   'FontSize', FS, 'Interpreter', 'latex');
        end
        if r == 1, title(ax, panels(c).title, 'FontSize', FS+1, 'FontWeight', 'normal', ...
                         'Interpreter', 'latex'); end
        % shared colourbar: drawn once per column, under the bottom row
        if r == nr
            cb = colorbar(ax, 'southoutside', 'Units', 'centimeters', 'FontSize', FS);
            cb.Position = [x0 1.15 axw 0.30];
            cb.Label.String = panels(c).cbar;
            cb.Label.FontSize = FS;
            ax.Units = 'centimeters';          % colorbar can shrink the axes
            ax.Position = [x0 y0 axw axh];
        end
    end
end

% the white and red isolines are defined in the LaTeX caption, so no key here

exportgraphics(fig, outfile, 'Resolution', RES);
close(fig);
end
