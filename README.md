# Project Shared Plotting

This repository contains the plotting scripts and the required data structure for reproducing the figures.

## 1. Project structure

The project folder should be organized as follows:

```text
project_shared_plotting/
├── data/
├── codes/
├── figures/
├── logs/
├── environment_minimal_used.yml
├── README.md
└── data_manifest.md
```

### `data/`

This folder stores all input data used by the plotting scripts.

```text
data/
├── behavior/
│   └── Ave_beha_normal.csv
│
├── LMM/
│   ├── corridor/
│   └── turn/
│
├── Lasso/
│   ├── corridor/
│   └── turn/
│
├── AAL3_data/
│   └── roi_cmap.json
│
├── spectrum/
│   └── *_spectrum_data.json
│
├── roi_abbreviations.json
├── feature_names.json
├── PLV_results_8_28_corridor_keyresults.json
├── PLV_results_8_28_turn_keyresults.json
├── PLV_results_8_28_corridor_end_split_keyresults.json
└── PLV_results_8_28_turn_end_split_keyresults.json
```

The exact files required by each figure are listed in `data_manifest.md`.

### `codes/`

This folder contains all plotting scripts and helper functions.

```text
codes/
├── project_paths.py
├── run_all_figures.py
│
├── Figure1D.py
├── Figure4A.py
├── Figure4B.py
├── Figure4C.py
├── Figure4D.py
├── Figure5A.py
├── Figure5A_legend.py
├── Figure5B.py
├── Figure6A.py
├── Figure6B.py
│
├── Sup_fig1.py
├── Sup_fig2a.py
├── Sup_fig2b.py
├── Sup_fig2c.py
├── Sup_fig2d.py
├── Sup_fig3.py
├── Sup_fig4.py
├── Sup_fig5.py
├── Sup_fig6.py
├── Sup_fig7a.py
├── Sup_fig7b.py
├── Sup_fig7c.py
├── Sup_fig7d.py
├── Sup_fig8a.py
├── Sup_fig8b.py
│
└── functions/
    ├── __init__.py
    ├── my_plot_setup.py
    └── figure_helpers.py
```

`project_paths.py` defines the project paths. The plotting scripts use this file to locate the input data and output folders.

### `figures/`

This folder stores the generated figures. Each major figure has its own subfolder.

```text
figures/
├── Figure1/
├── Figure4/
├── Figure5/
├── Figure6/
├── SupFigure1/
├── SupFigure2/
├── SupFigure3/
├── SupFigure4/
├── SupFigure5/
├── SupFigure6/
├── SupFigure7/
└── SupFigure8/
```

For example:

```text
Figure4A.py, Figure4B.py, Figure4C.py, Figure4D.py
```

will save their outputs to:

```text
figures/Figure4/
```

### `logs/`

This folder stores the log files generated when running all plotting scripts.

---

## 2. Environment setup

Activate the Conda environment used for plotting:

```cmd
conda activate DiversityMaping
```

If the environment needs to be recreated on another computer, use:

```cmd
conda env create -f environment_minimal_used.yml
```

If the environment already exists and needs to be updated:

```cmd
conda env update -f environment_minimal_used.yml --prune
```

---

## 3. Run all figures

Go to the project root folder:

```cmd
cd /d "D:\下载\project_shared_plotting_updated\project_shared_plotting"
```

Run all plotting scripts:

```cmd
python ".\codes\run_all_figures.py"
```

The output figures will be saved under:

```text
figures/
```

The running log will be saved under:

```text
logs/
```

---

## 4. Run selected figures

Run one single figure script:

```cmd
python ".\codes\Figure1D.py"
```

Run another example:

```cmd
python ".\codes\Figure4A.py"
```

Run one supplementary figure:

```cmd
python ".\codes\Sup_fig1.py"
```

---

## 5. Run only main figures or supplementary figures

Run only main figures:

```cmd
python ".\codes\run_all_figures.py" --only main
```

Run only supplementary figures:

```cmd
python ".\codes\run_all_figures.py" --only sup
```

Stop immediately when the first error occurs:

```cmd
python ".\codes\run_all_figures.py" --stop-on-error
```

---

## 6. Notes

The scripts use a non-interactive Matplotlib backend. Figures will be saved directly to disk and will not be displayed interactively.

The plotting style, colors, layouts, and statistical logic are defined inside the individual plotting scripts and helper functions.

Input paths and output paths are handled by:

```text
codes/project_paths.py
```

If the project folder is moved to a different location, the scripts should still work as long as the internal folder structure remains unchanged.

## Recent analysis and formal figures (2026-10-11)

Current PAC, Power and Granger source additions, dependencies, reporting conventions and reproduction limits are documented in [the recent snapshot guide](docs/RECENT_ANALYSIS_AND_FIGURES_20261011.md). Existing figure scripts and data-manifest content are retained. No protected input data or credentials are added by this update.

Formal rendered PNGs are managed separately in [the figure gallery repository](https://github.com/sonia1maslova/seeg-navigation-figure-gallery), with [the online gallery](https://sonia1maslova.github.io/seeg-navigation-figure-gallery/).
