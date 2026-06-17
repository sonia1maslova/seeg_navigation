import matplotlib
matplotlib.use("Agg", force=True)
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
plt.ioff()
from functions.my_plot_setup import set_plot
from pathlib import Path
import json

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#------------------------------------------
data_path = DATA_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

#------------------------------------------
process_type = 'turn'
lasso_data_path = data_path / 'LMM' / process_type
encoding_type = 'order+end'

# _selected 仅用于确定绘制哪些 ROI，完整数据用于绘图
signif_lasso_results_df = pd.read_csv(lasso_data_path / ('LMM_' + encoding_type + '_selected.csv'))
lasso_results_df = pd.read_csv(lasso_data_path / ('LMM_' + encoding_type + '.csv'))

#------------------------------------------
if process_type == 'corridor':
    n_stages = 7
    process_label = 'Corridor'
elif process_type == 'turn':
    n_stages = 6
    process_label = 'Junction'

if encoding_type == 'order+end':
    stages = ['_'.join((process_type, str(i))) for i in range(1, n_stages)] + [process_type + '_end']
    stages_label = [str(i) for i in range(1, n_stages)] + ['last']
elif encoding_type == 'order_only':
    stages = ['_'.join((process_type, str(i))) for i in range(1, n_stages + 1)]
    stages_label = [str(i) for i in range(1, n_stages + 1)]
stages = np.array(stages)

#------------------------------------------
rois = signif_lasso_results_df['roi'].unique().tolist()

#------------------------------------------
roi_cmap_df = pd.read_json(AAL3_DIR / 'roi_cmap.json')

with open(DATA_DIR / 'roi_abbreviations.json', 'r', encoding='utf-8') as f:
    roi_abbreviations = json.load(f)

#---------------------------------------------------------
if_save = 1
savefig_type = 'png'
figsize = (3.2, 2.4)

bin_width = 0.65
fontsize = 16

xlim = [0.5, n_stages + 0.5]
ylim = (-0.5, 1.6)

ylabel = 'LME coefficient (β)'
xlabel = process_label + ' number'

star_y = 1.38

#---------------------------------------------------------
plt.close('all')
for roi in rois:
    color_this_roi = roi_cmap_df[roi.replace(' ', '_')].values
    data_this_roi_df = lasso_results_df[lasso_results_df['roi'] == roi]  # 完整数据

    beta_values = []
    corr_p_corrected_values = []
    p_values = []

    for i in range(n_stages):
        row = data_this_roi_df[data_this_roi_df['segment'] == stages[i]]
        if not row.empty:
            beta_values.append(row['fixed_effect'].values.item())
            corr_p_corrected_values.append(row['corr_p_corrected'].values.item())
            p_values.append(row['fixed_effect_p'].values.item())
        else:
            beta_values.append(np.nan)
            corr_p_corrected_values.append(np.nan)
            p_values.append(np.nan)

    #---------------------------------------------------------
    fig, ax = set_plot(figsize=figsize, fontsize=fontsize)
    x = list(range(1, n_stages + 1))

    colors = [
        color_this_roi if (not np.isnan(b) and b > 0 and corr_p < 0.05 and p < 0.05) else [0.7, 0.7, 0.7]
        for b, corr_p, p in zip(beta_values, corr_p_corrected_values, p_values)
    ]
    plot_values = [v if not np.isnan(v) else 0 for v in beta_values]

    ax.bar(x, plot_values, width=bin_width, color=colors, alpha=1)

    for xi, beta, corr_p, p in zip(x, beta_values, corr_p_corrected_values, p_values):
        if np.isnan(beta) or beta <= 0:
            continue
        if corr_p < 0.05 and p < 0.001:
            stars = '***'
        elif corr_p < 0.05 and p < 0.01:
            stars = '**'
        elif corr_p < 0.05 and p < 0.05:
            stars = '*'
        else:
            stars = ''

        if stars:
            ax.text(xi, star_y, stars, ha='center', va='bottom',
                    fontsize=16, fontweight='bold', color='black')

    plt.xlim(xlim)
    plt.ylim(ylim)
    ax.set_xticks(x)
    ax.set_xticklabels(stages_label, rotation=0, ha='center')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plt.title(roi_abbreviations.get(roi, roi), pad=15, fontsize=fontsize)
    plt.tight_layout(pad=0.2)

    if if_save == 1:
        plt.savefig(
            out_path / ('_'.join((roi, process_type, encoding_type, 'beta')) + '.' + savefig_type),
            format=savefig_type, dpi=600, bbox_inches='tight'
        )
    # plt.show()  # disabled: figures are saved directly