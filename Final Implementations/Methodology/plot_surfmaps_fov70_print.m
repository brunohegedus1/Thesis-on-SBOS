%PLOT_SURFMAPS_FOV70_PRINT  Print-ready surface maps for Figure 4.8.
%
%   Re-draws the two panels of Figure 4.8 (Section 4.3.3, effect of field of
%   view) from the sweeps already solved by SETUP_SIZING_TEST_LMINUS and
%   SETUP_SIZING_TEST_LPLUS with L_FoV_obj = 7 cm, stored in
%   data_lminus_LFOV70.mat and data_lplus_LFOV70.mat. Both carry
%   The red isoline is drawn at s_min = 0.0004 m, the same value as the 2 cm
%   maps of Chapter 4, so the two field of view cases can be read against each
%   other directly. The .mat files carry smin_max = 0.00035 m from the sweep
%   that produced them; that value is no longer used for the isoline.
%
%   Same treatment as the 2 cm maps, with one difference: Figure 4.8 stacks
%   two full-width panels inside a single float, so each sheet is built
%   shorter (PANEL_H below) to keep the pair plus its captions on one page.
%   Everything else -- arrays, contour levels, colour limits, isolines -- is
%   unchanged.

clear; clc;

HERE    = fileparts(mfilename('fullpath'));
OUT_DIR = fullfile(HERE, 'print figures');
RES     = 300;                      % dpi
TEXTW   = 15.75;                    % \textwidth [cm]
FS      = 8;                        % base font size [pt], as the captions
PANEL_H = 1.85;                     % axes height [cm]; see note above
ISO_SMIN = 0.0004;                  % red isoline [m], same value as the 2 cm maps

if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end

cases = { ...
    'data_lminus_LFOV70.mat', 'Sensitivity and smin SBOS surf map 70cm.jpg'; ...
    'data_lplus_LFOV70.mat',  'Sensitivity and smin BOS surf map 70mm.jpg'};

q = @(A, p) quantile(A(:), p);

for i = 1:size(cases, 1)
    D = load(fullfile(HERE, cases{i,1}));
    smin = D.smin_data * 1e3;       % mm
    panels = struct( ...
        'data',  {D.S_data, smin}, ...
        'title', {'Sensitivity $S$', 'Resolution $s_{min}$'}, ...
        'cbar',  {'S [m]', 's_{min} [mm]'}, ...
        'clim',  {[0 q(D.S_data, 0.75)], [0 q(smin, 0.75)]}, ...
        'nlev',  {30, 30});
    surfmap_sheet(panels, D.l_data, D.f_data, D.f_num_data, D.S_data, D.smin_data, ...
        0.02, ISO_SMIN, TEXTW, PANEL_H, FS, fullfile(OUT_DIR, cases{i,2}), RES);
end

fprintf('written to %s\n', OUT_DIR);
