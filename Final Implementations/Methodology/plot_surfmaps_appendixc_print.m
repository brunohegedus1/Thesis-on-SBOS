%PLOT_SURFMAPS_APPENDIXC_PRINT  Print-ready surface maps for Appendix C.
%
%   Re-draws the six supplementary maps from the sweeps already solved by
%   SETUP_SIZING_TEST_LMINUS and SETUP_SIZING_TEST_LPLUS:
%
%     C.1  pixel density and in-focus FoV, negative defocus, 2 cm
%     C.2  speckle size, positive defocus, 2 cm
%     C.3  M_ag, d_f, m, negative defocus, 7 cm
%     C.4  M_ag, m, d_f, positive defocus, 7 cm
%     C.5  speckle size, negative defocus, 7 cm   (subfigure a)
%     C.6  speckle size, positive defocus, 7 cm   (subfigure b)
%
%   Same arrays, contour levels, colour limits and isolines as the originals;
%   only the layout changes. SURFMAP_SHEET does the drawing, so the appendix
%   matches Chapter 4.
%   The red isoline is drawn at s_min = 0.0004 m throughout, matching Chapter 4
%   including the 7 cm maps of Figure 4.8, so every surface map in the thesis
%   carries the same two reference isolines. The 7 cm .mat files still hold
%   smin_max = 0.00035 m from the sweep that produced them; it is not used.
%
%   C.5 and C.6 are the two halves of one float and are included at
%   0.9\textwidth each, so they are built narrower and shorter than the rest.

clear; clc;

HERE    = fileparts(mfilename('fullpath'));
OUT_DIR = fullfile(HERE, 'print figures');
RES     = 300;                      % dpi
TEXTW   = 15.75;                    % \textwidth [cm]
FS      = 8;                        % base font size [pt]; matches the LaTeX
                                    % caption size (footnotesize on a 10pt class)

if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end
q = @(A, p) quantile(A(:), p);
out = @(name) fullfile(OUT_DIR, name);

Lm20 = load(fullfile(HERE, 'data_lminus_LFOV20.mat'));
Lp20 = load(fullfile(HERE, 'data_lplus_LFOV20.mat'));
Lm70 = load(fullfile(HERE, 'data_lminus_LFOV70.mat'));
Lp70 = load(fullfile(HERE, 'data_lplus_LFOV70.mat'));

%% ---- C.1: pixel density and in-focus field of view (negative, 2 cm) ---
ppx = Lm20.Px_pitch * 1e-3;         % 1/mm
panels = struct( ...
    'data',  {ppx, Lm20.LFOV_data}, ...
    'title', {'Pixel density $\rho_{px}$', 'In-focus FoV $L_{FoV}$'}, ...
    'cbar',  {'\rho_{px} [1/mm]', 'L_{FoV} [m]'}, ...
    'clim',  {[0 q(ppx, 0.75)], [0 Lm20.L_FOV_obj]}, ...
    'nlev',  {50, 20});
sheet(panels, Lm20, TEXTW, 3.2, FS, out('Px pitch and LFOV SBOS surf map.jpg'), RES);

%% ---- C.2: speckle size (positive defocus, 2 cm) -----------------------
panels = struct( ...
    'data',  {Lp20.speckle_size}, ...
    'title', {'Speckle size $\Delta s$'}, ...
    'cbar',  {'\Delta s [px]'}, ...
    'clim',  {[0 q(Lp20.speckle_size, 0.75)]}, ...
    'nlev',  {20});
sheet(panels, Lp20, TEXTW, 3.1, FS, out('Speckle size surf map BOS.jpg'), RES);

%% ---- C.3: M_ag, d_f, m (negative defocus, 7 cm) -----------------------
panels = struct( ...
    'data',  {Lm70.M_data, Lm70.df_data, Lm70.m_data}, ...
    'title', {'Magnification $M_{ag}$', 'Focus distance $d_f$', 'Object distance $m$'}, ...
    'cbar',  {'M_{ag} [-]', 'd_f [m]', 'm [m]'}, ...
    'clim',  {[0 q(Lm70.M_data, 0.75)], [0 q(Lm70.df_data, 0.75)], [0 q(Lm70.m_data, 0.75)]}, ...
    'nlev',  {50, 20, 20});
sheet(panels, Lm70, TEXTW, 3.0, FS, out('M,m,df surfmaps SBOS 70mm.jpg'), RES);

%% ---- C.4: M_ag, m, d_f (positive defocus, 7 cm) -----------------------
panels = struct( ...
    'data',  {Lp70.M_data, Lp70.m_data, Lp70.df_data}, ...
    'title', {'Magnification $M_{ag}$', 'Object distance $m$', 'Focus distance $d_f$'}, ...
    'cbar',  {'M_{ag} [-]', 'm [m]', 'd_f [m]'}, ...
    'clim',  {[0 q(Lp70.M_data, 0.75)], [0 q(Lp70.m_data, 0.75)], [0 q(Lp70.df_data, 0.75)]}, ...
    'nlev',  {50, 20, 20});
sheet(panels, Lp70, TEXTW, 3.0, FS, out('M,m,df surfmaps BOS 70mm.jpg'), RES);

%% ---- C.5 and C.6: speckle size at 7 cm, two halves of one float -------
SUBW = 0.9 * TEXTW;                 % both are included at 0.9\textwidth
spk = {Lm70, 'Speckle size surf map SBOS 70mm.jpg'; ...
       Lp70, 'Speckle size surf map BOS 70mm.jpg'};
for i = 1:size(spk, 1)
    D = spk{i,1}; name = spk{i,2};
    panels = struct( ...
        'data',  {D.speckle_size}, ...
        'title', {'Speckle size $\Delta s$'}, ...
        'cbar',  {'\Delta s [px]'}, ...
        'clim',  {[0 q(D.speckle_size, 0.75)]}, ...
        'nlev',  {20});
    sheet(panels, D, SUBW, 1.85, FS, out(name), RES);
end

fprintf('written to %s\n', OUT_DIR);

%% ---- thin wrapper: pulls the grid and isolines out of a dataset -------
function sheet(panels, D, W, axh, FS, outfile, RES)
surfmap_sheet(panels, D.l_data, D.f_data, D.f_num_data, D.S_data, D.smin_data, ...
    0.02, 0.0004, W, axh, FS, outfile, RES);
end
