# Final Implementations

The MATLAB scripts that produce the results in the thesis, one folder per
experiment. Each script was copied from its working folder; the originals are
untouched, so every script here also exists where it came from.

The data each script reads sits next to it, so a script can be run from its
folder as it stands and no path needs editing. See 'Data' at the end of this
file.

## Methodology

Design space exploration of Chapter 4 and the surface maps of Appendix C.
Taken from `SBOS methodology/it7`, the last iteration.

| Script | Role |
|---|---|
| `setup_sizing_test_lminus.m` | Sweep with the object in front of the in-focus plane |
| `setup_sizing_test_lplus.m` | Sweep with the object behind the in-focus plane |
| `plot_surfmaps_lminus_print.m` | Surface maps for the negative-defocus sweep |
| `plot_surfmaps_lplus_print.m` | Surface maps for the positive-defocus sweep |
| `plot_surfmaps_fov70_print.m` | Surface maps for the 70 mm field of view |
| `plot_surfmaps_appendixc_print.m` | The maps collected in Appendix C |
| `surfmap_sheet.m` | Helper that lays out one printed sheet of maps |

## Angled glass experiment

| Script | Role |
|---|---|
| `plot_sensitivity.m` | Sensitivity against defocus distance, both campaigns. Figure 6.1 |

Set `campaign` to `'setup-matrix'` with `groupBy = 'fstop'` for panel (a), and
to `'slider-sweep'` with `groupBy = 'defocus'` for panel (b).

## Speckle Pattern investigation

Entry points are `speckle_size.m` for one map and `speckle_size_batch.m` for a
folder. The three uncertainty scripts and the two plotting scripts read the
batch output.

| Script | Role |
|---|---|
| `speckle_size.m` | Mean speckle size of one intensity map |
| `speckle_size_batch.m` | The same for every export in a folder |
| `speckle_estimators.m` | Autocorrelation and width routines shared by both |
| `read_intensity_map.m` | Reads a DaVis intensity export |
| `sort_results_file.m` | Sorts the batch results by point number |
| `speckle_size_sigma_band.m` | Uncertainty of the flattening width |
| `speckle_size_sensitivity.m` | Response of the measured size to that width |
| `speckle_size_budget.m` | Post-processing uncertainty, point by point |
| `plot_speckle_methods.m` | Comparison of the size estimators |
| `plot_speckle_trends.m` | Size against f-number, magnification and diffuser distance |

## Compressible jet experiment

The chain runs `read_bos_txt` to `bos_pipeline` to `bos_optical_path` to
`abel_invert_columns`, which `run_bos_density.m` drives end to end.

| Script | Role |
|---|---|
| `read_bos_txt.m` | Reads a DaVis displacement export onto a grid |
| `bos_data_dir.m` | Single place that holds the path to the exports |
| `bos_background_offset.m` | Background displacement offset |
| `bos_pipeline.m` | Steps 1 to 8: background, centreline, straighten, crop, fold, Abel |
| `bos_optical_path.m` | Deflection to projected optical path |
| `abel_invert_columns.m` | Abel inversion with Tikhonov regularization |
| `extract_centreline.m` | Centreline extraction helper |
| `run_bos_density.m` | Raw export to density and density gradient |
| `run_bos_gradients_unfolded.m` | Planar gradients on the unfolded grid, centreline gradient |
| `Ls_shock_cell_length.m` | Shock cell length against NPR, Emden and Hartmann-Lazarus fits |
| `plot_cc_displacement_fields.m` | Displacement fields of the six datasets |
| `plot_centreline_panda.m` | Centreline density comparison |
| `plot_displacement_histogram.m` | Displacement histograms |

## Aerospike experiment

| Script | Role |
|---|---|
| `aerospike_density_gradient.m` | Density gradient fields from the displacement exports |
| `aerospike_speckle_size.m` | Speckle size of the aerospike campaign |
| `plot_aerospike_data.m` | Plots of the aerospike measurements |
| `measure_distance.m` | Interactive distance measurement, used for the pixel density |

`aerospike_speckle_size.m` needs `speckle_estimators.m` and
`read_intensity_map.m`. Its `opt.libFolder` points at the sibling folder
`Speckle Pattern investigation`, so the two experiments share one
implementation. That line is the only difference from the original.

## Data

This collection sits in the thesis repository, where `.gitignore` keeps the
large exports out: GitHub takes no file over 100 MB, and these run to 160 MB.
The scripts and the data under 1 MB are tracked; the exports listed below are
not, so a fresh clone needs them copied back.

Every folder holds the data its scripts read, 7.4 GB in total. The duplicate
copies left in the working folders were deleted after every file was checked
against its copy with md5, so this collection now holds the only copy of the
exports. Back it up as it stands.

| Folder | Data | Size | Read by |
|---|---|---|---|
| Methodology | `data_lminus_LFOV20.mat`, `data_lplus_LFOV20.mat`, `data_lminus_LFOV70.mat`, `data_lplus_LFOV70.mat` | 0.5 MB | the four `plot_surfmaps_*_print.m` scripts |
| Angled glass experiment | `data/angled_glass_experimental.csv`, `data/angled_glass_errors.csv` | 2 KB | `plot_sensitivity.m` |
| Speckle Pattern investigation | `speckle_size_batch.csv`, `test_matrix.csv`, `speckle_size_sensitivity.csv`, `speckle_size_sigma_band.csv`, `speckle_size_budget.csv` | 12 KB | the three uncertainty scripts and `plot_speckle_trends.m` |
| Speckle Pattern investigation | `I4.csv` to `I31.csv` and `I1_2.txt` to `I3_2.txt`, the intensity maps of the 31 points, with the parsed cache `I3_2.mat` | 4.2 GB | `speckle_size_batch.m` and `speckle_size.m` |
| Speckle Pattern investigation | `Exp1/intensity_*.txt` and `intensity_*.mat` | 1.3 GB | nothing in this collection. Kept with the campaign it belongs to |
| Compressible jet experiment | `CC results/BOS_3cm_12x120001*.csv`, six runs | 259 MB | `run_bos_density.m`, `run_bos_gradients_unfolded.m`, `plot_cc_displacement_fields.m`, `plot_displacement_histogram.m` |
| Compressible jet experiment | `BOS_data/BOS_3cm_{3..8}bar_12x12.txt` and the `_M0` set | 519 MB | the same scripts with `DATA_SOURCE = 'txt'`, through `bos_data_dir.m` |
| Compressible jet experiment | `Emden_data.csv`, `Pandas data.csv`, `centreline.mat` | 60 KB | `Ls_shock_cell_length.m`, `plot_centreline_panda.m` |
| Aerospike experiment | `Aerospike Data/SBOS1..11.csv` | 333 MB | `aerospike_density_gradient.m`, `plot_aerospike_data.m` |
| Aerospike experiment | `Aerospike Speckle Data/REF*.csv` and `cache/REF*.mat` | 868 MB | `aerospike_speckle_size.m`. The cache comes along so the 150 MB files are not parsed again |

The copies resolve these paths against their own folder, which is the only
change made to them besides the `opt.libFolder` line noted above:

- `bos_data_dir.m` returns the local `BOS_data` folder.
- `CC_RESULTS_DIR` and `CC_DIR` point at the local `CC results` folder.
- `extract_centreline.m` reads `run_bos_density.m` from its own folder.
- `speckle_size.m` and `speckle_size_batch.m` read the intensity maps from
  their own folder.

`Ls_shock_cell_length.m` keeps `EMDEN_CSV` empty, as in the run that
produced Figure 6.11, so the digitised Emden points are not drawn. Point it
at the `Emden_data.csv` in the folder to add them.

## What was left out, and why

| Left out | Reason |
|---|---|
| `SBOS methodology/it5`, `it6` | Earlier iterations of the same sweeps. `it7` supersedes them |
| `SBOS methodology/it5/setup_sizing_test_M.m` | Magnification sweep with no `it7` counterpart; the maps in the thesis come from the `it7` scripts |
| `nfbos_van hinsberg/` | Earlier copies of `read_bos_txt`, `bos_pipeline`, `abel_invert_columns` and `run_bos_pipeline`, plus `nfbos_density_field.m`, which belongs to the older density route that `run_bos_density.m` warns against |
| `BOS jet Post processing/density_gradient.m` | Ad-hoc comparison of two text files, superseded by `run_bos_density.m` |
| `BOS jet Post processing/run_bos_pipeline.m` | Diagnostic driver whose chain is a subset of `run_bos_density.m` |
| `moc jet/` | Method of characteristics solver and its validation. The method appears in Chapter 2 only as theory, and no figure in the thesis comes from it |
| `MIRAGE/` | Ray tracing simulation. MIRAGE is named only in the introduction and produces no result in the thesis |

Any of these can be added back: the script files are still in their folders.
Their data is not -- the raw exports now live here only.
