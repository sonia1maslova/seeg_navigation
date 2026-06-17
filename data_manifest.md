# Data manifest

Place data files under `data/` using this structure.

```text
data/
├── behavior/
│   └── Ave_beha_normal.csv
├── LMM/
│   ├── corridor/
│   └── turn/
├── Lasso/
│   ├── corridor/
│   └── turn/
├── AAL3_data/
│   └── roi_cmap.json
├── spectrum/
│   └── *_spectrum_data.json
├── roi_abbreviations.json
├── feature_names.json
├── feature_names_6imfs.json
└── PLV_results*_keyresults.json
```

The helper script `Codes/collect_data_to_shared_structure.py` can copy the expected files from the old scattered folders.
