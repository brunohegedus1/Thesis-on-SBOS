%PLOT_SURFMAPS_LPLUS_PRINT  Print-ready surface maps for Section 4.3.2.
%
%   Re-draws Figures 4.6 and 4.7 of the thesis from the sweep already solved
%   by SETUP_SIZING_TEST_LPLUS, stored in data_lplus_LFOV20.mat (positive
%   defocus, L_FoV_obj = 2 cm, the value quoted in the captions).
%
%   Same deal as PLOT_SURFMAPS_LMINUS_PRINT: the maths is untouched (same
%   arrays, contour levels, colour limits and isolines), only the layout
%   changes so the text survives being scaled to \textwidth. SURFMAP_SHEET
%   does the drawing, so both sections stay consistent.
%
%   Note the column order here is M_ag, m, d_f, which is the order the
%   original figure and the Figure 4.7 caption use; the negative defocus
%   figure uses M_ag, d_f, m.

clear; clc;

HERE    = fileparts(mfilename('fullpath'));
DATA    = fullfile(HERE, 'data_lplus_LFOV20.mat');
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
mdist = D.m_data;
df    = D.df_data;
isoS    = 0.02;                     % white isoline [m]
isoSmin = D.smin_max;               % red isoline [m], 0.0004 m in this run

q = @(A, p) quantile(A(:), p);

%% ---- Figure 4.6: sensitivity and s_min -------------------------------
panels = struct( ...
    'data',  {S, smin}, ...
    'title', {'Sensitivity $S$', 'Resolution $s_{min}$'}, ...
    'cbar',  {'S [m]', 's_{min} [mm]'}, ...
    'clim',  {[0 q(S, 0.75)], [0 q(smin, 0.75)]}, ...
    'nlev',  {30, 30});
surfmap_sheet(panels, l, f, fnum, S, D.smin_data, isoS, isoSmin, TEXTW, 1.85, FS, ...
    fullfile(OUT_DIR, 'Sensitivity and smin BOS surf map.jpg'), RES);

%% ---- Figure 4.7: magnification and distances -------------------------
panels = struct( ...
    'data',  {Mag, mdist, df}, ...
    'title', {'Magnification $M_{ag}$', 'Object distance $m$', 'Focus distance $d_f$'}, ...
    'cbar',  {'M_{ag} [-]', 'm [m]', 'd_f [m]'}, ...
    'clim',  {[0 q(Mag, 0.75)], [0 q(mdist, 0.75)], [0 q(df, 0.75)]}, ...
    'nlev',  {50, 20, 20});
surfmap_sheet(panels, l, f, fnum, S, D.smin_data, isoS, isoSmin, TEXTW, 1.85, FS, ...
    fullfile(OUT_DIR, 'M,m,df surfmaps BOS.jpg'), RES);

fprintf('written to %s\n', OUT_DIR);
