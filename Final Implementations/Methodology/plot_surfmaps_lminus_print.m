%PLOT_SURFMAPS_LMINUS_PRINT  Print-ready surface maps for Section 4.3.1.
%
%   Re-draws Figures 4.3, 4.4 and 4.5 of the thesis from the sweep already
%   solved by SETUP_SIZING_TEST_LMINUS, stored in data_lminus_LFOV20.mat
%   (L_FoV_obj = 2 cm, the value quoted in the captions).
%
%   The maths is untouched: same arrays, same contour levels, same colour
%   limits and the same two isolines (white S = 0.02 m, red s_min = 0.4 mm).
%   Only the layout changes, so the figures stay readable once LaTeX scales
%   them to \textwidth:
%     - the figure is built at its final printed size (15.75 cm = \textwidth
%       for a4paper with hscale=0.75), so 9 pt here is 9 pt on the page;
%     - one shared horizontal colourbar per column instead of one per tile,
%       which gives the panels roughly twice the width they had;
%     - f_# labels moved to row headers, axis labels only on the outer
%       panels, fewer ticks;
%     - the filled contours are drawn without edge lines, which removes the
%       dashed-line moire of the originals.
%
%   Output file names match the ones report.tex already includes, so nothing
%   in the LaTeX source has to change.

clear; clc;

HERE    = fileparts(mfilename('fullpath'));
DATA    = fullfile(HERE, 'data_lminus_LFOV20.mat');
OUT_DIR = fullfile(HERE, 'print figures');
RES     = 300;                      % dpi
TEXTW   = 15.75;                    % \textwidth [cm]
FS      = 8;                        % base font size [pt]; matches the LaTeX
                                    % caption size (footnotesize on a 10pt class)

if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end
D = load(DATA);

l = D.l_data; f = D.f_data; fnum = D.f_num_data;
S     = D.S_data;
smin  = D.smin_data * 1e3;          % mm
Mag   = D.M_data;
df    = D.df_data;
mdist = D.m_data;
dspk  = D.speckle_size;
isoS    = 0.02;                     % white isoline [m]
isoSmin = D.smin_max;               % red isoline [m]

q = @(A, p) quantile(A(:), p);

%% ---- Figure 4.3: sensitivity and s_min -------------------------------
panels = struct( ...
    'data',  {S, smin}, ...
    'title', {'Sensitivity $S$', 'Resolution $s_{min}$'}, ...
    'cbar',  {'S [m]', 's_{min} [mm]'}, ...
    'clim',  {[0 q(S, 0.75)], [0 q(smin, 0.75)]}, ...
    'nlev',  {30, 30});
surfmap_sheet(panels, l, f, fnum, S, D.smin_data, isoS, isoSmin, TEXTW, 1.85, FS, ...
    fullfile(OUT_DIR, 'Sensitivity and dav SBOS surf map.jpg'), RES);

%% ---- Figure 4.4: magnification and distances -------------------------
panels = struct( ...
    'data',  {Mag, df, mdist}, ...
    'title', {'Magnification $M_{ag}$', 'Focus distance $d_f$', 'Object distance $m$'}, ...
    'cbar',  {'M_{ag} [-]', 'd_f [m]', 'm [m]'}, ...
    'clim',  {[0 q(Mag, 0.75)], [0 q(df, 0.75)], [0 q(mdist, 0.75)]}, ...
    'nlev',  {50, 20, 20});
surfmap_sheet(panels, l, f, fnum, S, D.smin_data, isoS, isoSmin, TEXTW, 1.85, FS, ...
    fullfile(OUT_DIR, 'M,m,df surfmaps SBOS.jpg'), RES);

%% ---- Figure 4.5: speckle size ----------------------------------------
panels = struct( ...
    'data',  {dspk}, ...
    'title', {'Speckle size $\Delta s$'}, ...
    'cbar',  {'\Delta s [px]'}, ...
    'clim',  {[0 q(dspk, 0.9)]}, ...
    'nlev',  {20});
surfmap_sheet(panels, l, f, fnum, S, D.smin_data, isoS, isoSmin, TEXTW, 1.85, FS, ...
    fullfile(OUT_DIR, 'Speckle size surf map SBOS.jpg'), RES);

fprintf('written to %s\n', OUT_DIR);

