"""
granger_posthoc_config.py
=========================

Central configuration and path helpers for Granger post-processing.

Main output rule
----------------
All post-processing results are saved under:

    output_base_path / parameter_folder / process_type / analysis_output_folder

The filenames inside those folders are intentionally short and do NOT include
the full parameter combination. The parameter combination is represented by
the parameter_folder and by params_manifest.json.

Input rule
----------
Original Granger result files are read from an absolute root:

    input_root_path / process_type / allpairs

The exact input CSV filename is reconstructed from:
    analysis_kind + process_type + brain_area1 + brain_area2 + gc_params
"""

from pathlib import Path
import json
import math
import re


# =========================================================
# Main input/output roots
# =========================================================
from granger_runtime_config import load_settings
_settings = load_settings()
data_path = _settings['data_path']
input_root_path = _settings['input_root_path']
output_base_path = _settings['output_base_path']
analysis_kind = _settings['analysis_kind']
process_types = _settings['process_types']
brain_pairs = _settings['brain_pairs']
gc_params = _settings['gc_params']
alpha = _settings['posthoc']['alpha']
pair_agg_method = _settings['posthoc']['pair_agg_method']
contacts_csv = data_path / 'contacts_22rois.csv'
selected_data_cuts_by_process = None
DATA_CUTS_BY_ANALYSIS = _settings['data_cuts_by_analysis']

# Metric/state definitions shared by step1/step2/step3/edges
# =========================================================
METRIC_ORDER = [
    ('time_domain',   'time_domain'),
    ('spectral_band', 'slow_0p5_4Hz'),
    ('spectral_band', 'theta_alpha_4_12'),
    ('spectral_band', 'beta_12_30'),
    ('spectral_band', 'low_gamma_30_70'),
    ('spectral_band', 'high_gamma_70_150'),
]

METRIC_SHORT_LABELS = {
    ('time_domain',   'time_domain'):       'T',
    ('spectral_band', 'slow_0p5_4Hz'):         'D',
    ('spectral_band', 'theta_alpha_4_12'):  'TA',
    ('spectral_band', 'beta_12_30'):        'B',
    ('spectral_band', 'low_gamma_30_70'):   'LG',
    ('spectral_band', 'high_gamma_70_150'): 'HG',
}

METRIC_SHORT_NAMES = {
    ('time_domain',   'time_domain'):       'time',
    ('spectral_band', 'slow_0p5_4Hz'):         'delta',
    ('spectral_band', 'theta_alpha_4_12'):  'thetaalpha',
    ('spectral_band', 'beta_12_30'):        'beta',
    ('spectral_band', 'low_gamma_30_70'):   'lgamma',
    ('spectral_band', 'high_gamma_70_150'): 'hgamma',
}

STATE_CODE = {
    'missing':   0,
    'dir1_only': 1,
    'dir2_only': 2,
    'both':      3,
    'none':      4,
}


# =========================================================
# String helpers
# =========================================================
def make_safe_filename(s):
    """Compact safe name for output files/folders."""
    s = str(s)
    s = ''.join(ch if ord(ch) >= 32 else '_' for ch in s)
    s = re.sub(r'[\\/:*?"<>|]+', '_', s)
    s = re.sub(r'\s+', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_')


def make_matlab_file_token(s):
    """
    Safe token for exact MATLAB-generated filenames.

    Space becomes "_", but the "__" separators between tokens are preserved by
    make_gc_prefix() because that function joins already-cleaned tokens.
    """
    s = str(s)
    s = ''.join(ch if ord(ch) >= 32 else '_' for ch in s)
    s = re.sub(r'[\\/:*?"<>|]+', '_', s)
    s = re.sub(r'\s+', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_')


def analysis_slug(kind=None):
    kind = analysis_kind if kind is None else str(kind)
    slug = kind.replace('+', '_plus_')
    slug = re.sub(r'[^A-Za-z0-9]+', '_', slug)
    slug = re.sub(r'_+', '_', slug).strip('_')
    return slug


def compact_area_name(s):
    s = str(s)
    replacements = {
        'Temporal_Mid': 'TempMid',
        'Temporal_Sup': 'TempSup',
        'Hippocampus': 'Hipp',
        'ParaHippocampal': 'ParaHipp',
        'SupraMarginal': 'SupraMarg',
        'Precentral': 'Precentral',
        'OFCpost': 'OFCpost',
        '_Anterior': '_Ant',
        '_Posterior': '_Post',
    }
    for old, new in replacements.items():
        s = s.replace(old, new)
    s = s.replace(' L', 'L').replace(' R', 'R')
    s = s.replace(' ', '')
    s = s.replace('__', '_')
    return make_safe_filename(s)


def pair_key(brain_area1, brain_area2):
    return make_safe_filename(compact_area_name(brain_area1) + '_to_' + compact_area_name(brain_area2))


def compact_prefix(process_type, brain_area1, brain_area2, kind=None):
    """
    Short output prefix for files inside the parameter folder.

    It intentionally does NOT include process_type, analysis_kind, or GC params,
    because those are already represented by the folder structure.
    """
    return pair_key(brain_area1, brain_area2)


# =========================================================
# Input filename helpers
# =========================================================
def analysis_file_label(kind=None):
    """Return the exact analysis label segment used in raw all-pair filenames."""
    kind = analysis_kind if kind is None else str(kind)
    if kind.startswith('_'):
        return kind
    if kind in ['order_aligned', 'order_aligned_and_terminal_aligned', 'end_only']:
        return '_' + kind
    return '_' + kind.strip('_')


def surrogate_label(params=None):
    params = gc_params if params is None else params
    if params.get('surrogate_label'):
        return str(params['surrogate_label'])
    if params.get('do_shuffle', False):
        return 'trial_shuffle' + str(params.get('nshuffles', ''))
    return 'noshuffle'


def trialcap_label(params=None):
    params = gc_params if params is None else params
    max_trials = params.get('max_trials_per_pair', float('inf'))
    try:
        if math.isinf(float(max_trials)):
            return 'alltrials'
    except Exception:
        pass
    return 'maxtrials' + str(int(max_trials))


def morder_label(params=None, short=False):
    params = gc_params if params is None else params

    if str(params.get('model_order_mode', '')).lower() == 'fixed':
        return 'fixed' + str(params.get('fixed_morder', ''))

    if short:
        return 'momax' + str(params.get('momax', ''))

    return 'bic_momax' + str(params.get('momax', ''))


def make_gc_prefix(process_type, brain_area1, brain_area2, kind=None, params=None):
    """
    Build the exact MATLAB all-pair CSV prefix from parameters.

    Current short MATLAB format:
      corridor__Temporal_Mid_L_Anterior__to__Hippocampus_L_Anterior_order_aligned__allpairs__shuffle1000__momax20__shuffle
    """
    params = gc_params if params is None else params

    process_token = make_matlab_file_token(process_type)
    a1_token = make_matlab_file_token(brain_area1)
    a2_token = make_matlab_file_token(brain_area2)

    filename_format = str(params.get('filename_format', 'matlab_short')).lower()

    if filename_format == 'matlab_short':
        return (
            f'{process_token}__{a1_token}__to__{a2_token}{analysis_file_label(kind)}'
            f'__allpairs__{surrogate_label(params)}__{morder_label(params, short=True)}'
            f'__{params.get("parallel_level", "")}'
        )

    if filename_format == 'long':
        return (
            f'{process_token}__{a1_token}__to__{a2_token}{analysis_file_label(kind)}'
            f'__allpairs__{surrogate_label(params)}__{trialcap_label(params)}'
            f'__{morder_label(params, short=False)}__nfreq{params.get("nfreqs", "")}'
            f'__{params.get("parallel_level", "")}'
        )

    raise ValueError('Unknown gc filename_format: ' + str(filename_format))


def allpairs_dir(process_type):
    return input_root_path / str(process_type) / 'allpairs'


def pair_csv_path(process_type, brain_area1, brain_area2, kind=None, params=None):
    prefix = make_gc_prefix(process_type, brain_area1, brain_area2, kind, params)
    return allpairs_dir(process_type) / (prefix + '_all_channel_pairs.csv')


def full_mat_path(process_type, brain_area1, brain_area2, kind=None, params=None):
    prefix = make_gc_prefix(process_type, brain_area1, brain_area2, kind, params)
    return allpairs_dir(process_type) / (prefix + '_full.mat')


# =========================================================
# Unified output helpers
# =========================================================
def parameter_folder_name(kind=None, params=None):
    """
    Parameter-combination folder name. Filenames under this folder should stay short.
    """
    kind = analysis_kind if kind is None else str(kind)
    params = gc_params if params is None else params

    parts = [
        analysis_slug(kind),
        surrogate_label(params),
        morder_label(params, short=True),
        str(params.get('parallel_level', '')),
    ]

    if str(params.get('filename_format', 'matlab_short')).lower() == 'long':
        parts.insert(2, trialcap_label(params))
        parts.insert(-1, 'nfreq' + str(params.get('nfreqs', '')))

    return make_safe_filename('__'.join([p for p in parts if str(p) != '']))


def parameter_root(kind=None, params=None):
    return output_base_path / parameter_folder_name(kind, params)


def process_output_root(process_type, kind=None, params=None):
    return parameter_root(kind, params) / str(process_type)


def heatmap_dir(process_type, kind=None, params=None):
    return process_output_root(process_type, kind, params) / 'sig_heatmap'


def heatmap_data_dir(process_type, kind=None, params=None):
    return heatmap_dir(process_type, kind, params) / 'data'


def heatmap_fig_dir(process_type, kind=None, params=None):
    return heatmap_dir(process_type, kind, params) / 'fig'


def dice_dir(process_type, kind=None, params=None):
    return heatmap_dir(process_type, kind, params) / 'dice'


def edge_dir(process_type, kind=None, params=None):
    return process_output_root(process_type, kind, params) / 'sig_edges'


def surfice_fig_dir(process_type, kind=None, params=None):
    return process_output_root(process_type, kind, params) / 'surfice_figures'


def order_distribution_dir(process_type, kind=None, params=None):
    return process_output_root(process_type, kind, params) / 'order_distribution'


def null_hist_dir(process_type, pair_agg_method_this=None, kind=None, params=None):
    if pair_agg_method_this is None:
        pair_agg_method_this = pair_agg_method
    return process_output_root(process_type, kind, params) / ('null_histograms_pairagg_' + str(pair_agg_method_this))


# Backward-compatible folder-name helpers. These now return short folder names only.
def heatmap_folder_name(kind=None):
    return 'sig_heatmap'


def edge_folder_name(kind=None):
    return 'sig_edges'


def surfice_figure_folder_name(kind=None):
    return 'surfice_figures'


def order_distribution_folder_name(kind=None):
    return 'order_distribution'


def null_hist_folder_name(pair_agg_method_this=None, kind=None):
    if pair_agg_method_this is None:
        pair_agg_method_this = pair_agg_method
    return 'null_histograms_pairagg_' + str(pair_agg_method_this)


def params_manifest_dict():
    return {
        'data_path': str(data_path),
        'input_root_path': str(input_root_path),
        'output_base_path': str(output_base_path),
        'parameter_folder_name': parameter_folder_name(),
        'analysis_kind': str(analysis_kind),
        'process_types': list(process_types),
        'brain_pairs': list(brain_pairs),
        'selected_data_cuts_by_process': selected_data_cuts_by_process,
        'data_cuts_by_analysis': DATA_CUTS_BY_ANALYSIS,
        'gc_params': {k: str(v) if isinstance(v, float) and math.isinf(v) else v for k, v in gc_params.items()},
        'alpha': alpha,
        'pair_agg_method': pair_agg_method,
    }


def write_params_manifest():
    root = parameter_root()
    root.mkdir(parents=True, exist_ok=True)
    out_file = root / 'params_manifest.json'
    with open(out_file, 'w', encoding='utf-8') as f:
        json.dump(params_manifest_dict(), f, ensure_ascii=False, indent=2)
    return out_file


# =========================================================
# Data-cut helpers
# =========================================================
def get_selected_data_cuts(process_type, kind=None):
    kind = analysis_kind if kind is None else str(kind)
    if selected_data_cuts_by_process is not None:
        val = selected_data_cuts_by_process.get(process_type, None)
        if val is None:
            return None
        return list(val)
    return list(DATA_CUTS_BY_ANALYSIS.get(kind, {}).get(process_type, []))


def natural_data_cut_order(process_type, existing_cuts, kind=None):
    existing_cuts = [str(x) for x in existing_cuts]
    preferred = get_selected_data_cuts(process_type, kind)
    if preferred is None or len(preferred) == 0:
        preferred = []
    ordered = [x for x in preferred if x in existing_cuts]
    rest = sorted(x for x in existing_cuts if x not in ordered)
    return ordered + rest


def filter_dataframe_by_data_cut(df, process_type, col='data_cut', kind=None):
    cuts = get_selected_data_cuts(process_type, kind)
    if cuts is None or len(cuts) == 0 or col not in df.columns:
        return df.copy()
    return df[df[col].astype(str).isin(cuts)].copy()
