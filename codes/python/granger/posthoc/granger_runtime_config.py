"""Shared S1 settings for posthoc scripts; no analysis runs on import."""
import json
import os
from pathlib import Path

from granger_layout import CONFIG_PATH, PROJECT_ROOT


def load_settings():
    config_path = Path(os.environ.get('GRANGER_CONFIG', CONFIG_PATH))
    raw = json.loads(config_path.read_text(encoding='utf-8'))

    def path(key, env):
        value = Path(os.environ.get(env, raw[key]))
        return value if value.is_absolute() else PROJECT_ROOT / value

    params = {k: float('inf') if v == 'Inf' else v for k, v in raw['gc_params'].items()}
    if params['surrogate_method'] != 'trial_shuffle':
        raise ValueError('Granger surrogate_method must remain trial_shuffle')
    if raw['end_split'] not in (0, 1):
        raise ValueError('end_split must be 0 or 1')
    params.update(surrogate_label=f"trial_shuffle{params['nshuffles']}", filename_format='long')
    cuts = {
        'order_aligned': {'corridor': [f'corr{i}' for i in range(1, 8)],
                          'junction': [f'junction{i}' for i in range(1, 7)]},
        'order_aligned_and_terminal_aligned': {
            'corridor': [f'corr{i}' for i in range(1, 7)] + ['corr_end'],
            'junction': [f'junction{i}' for i in range(1, 6)] + ['junction_end']},
    }
    return {
        **raw,
        'config_path': str(config_path),
        'data_path': path('data_root', 'GRANGER_DATA_ROOT'),
        'input_root_path': path('results_root', 'GRANGER_RESULTS_ROOT'),
        'output_base_path': path('posthoc_root', 'GRANGER_POSTHOC_ROOT'),
        'analysis_kind': 'order_aligned_and_terminal_aligned' if raw['end_split'] else 'order_aligned',
        'brain_pairs': [{'brain_area1': a, 'brain_area2': b}
                        for a in raw['brain_area1s'] for b in raw['brain_area2s']],
        'gc_params': params,
        'data_cuts_by_analysis': cuts,
    }
