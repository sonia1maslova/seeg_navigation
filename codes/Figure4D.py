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
import matplotlib.colors as mcolors
import json

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#------------------------------------------
data_path = DATA_DIR
process_type = 'turn'
lasso_data_path = data_path / 'LMM' /process_type / 'Ante_Post_splited'

out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

#------------------------------------------
encoding_type = 'order+end'
lasso_results_df = pd.read_csv(lasso_data_path/ ('LMM_'+ encoding_type + '.csv'))

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
    stages = ['_'.join((process_type,  str(i))) for i in range(1, n_stages + 1)]
    stages_label = [str(i) for i in range(1, n_stages + 1)]
stages = np.array(stages)
#------------------------------------------
the_roi = 'Temporal_Mid L'
mask = (lasso_results_df['roi'].str.contains(the_roi))
lasso_results_df = lasso_results_df[mask]
rois = lasso_results_df['roi'].unique()
#----------------------------------------------------
with open(DATA_DIR / 'roi_abbreviations.json', 'r', encoding='utf-8') as f:
    roi_abbreviations = json.load(f)
#---------------------------------------------------------
if_save = 1
savefig_type = 'png'
figsize = (3.2, 2.4)

cmaps = ['#f85f73', '#ffb6b9']

linewidth = 1.5
linetype = '-'

xlim = [0.5, n_stages + 0.5]
ylim = (-0.2, 0.6)

axis_linewidth = 1.7
bin_width = 0.65

fontsize = 16
ylabel = 'Correlation'
xlabel = process_label + ' number'

xticks = list(range(1, 8))
yticks = [0.3, 0.6, 0.9, 1.2]

xticks_label = stages_label

#--------------------------------------------------
plt.close('all')
for idx_roi, roi in enumerate(rois):
    color_this_roi = cmaps[idx_roi]
    data_this_roi_df = lasso_results_df[lasso_results_df['roi'] == roi]
    r_values = []
    p_values_corrected = []
    LMM_beta_p_values = []
    for i in range(n_stages):
        current_df = data_this_roi_df[data_this_roi_df['segment'] == stages[i]]
        r = current_df['corr'].values.item()
        p_corrected = current_df['corr_p_corrected'].values.item()
        LMM_beta_p = current_df['fixed_effect_p'].values.item()
        r_values.append(r)
        p_values_corrected.append(p_corrected)
        LMM_beta_p_values.append(LMM_beta_p)

    #---------------------------------------------------------
    fig, ax = set_plot(figsize=figsize, fontsize=fontsize)
    x = list(range(1, n_stages+1))
    y = r_values

    # r > 0 且 corr_p < 0.05 才显示 ROI 颜色，否则灰色
    colors = [color_this_roi if (r > 0 and p < 0.05) else [0.7, 0.7, 0.7]
              for r, p in zip(r_values, p_values_corrected)]

    ax.bar(x, y, width=bin_width, color=colors, alpha=1)

    for xi, yi, r, corr_p in zip(x, y, r_values, p_values_corrected):
        # 只有 r > 0 且显著才标注星号
        if r > 0 and corr_p < 0.05:
            if corr_p < 0.001:
                stars = '***'
            elif corr_p < 0.01:
                stars = '**'
            else:
                stars = '*'
            ax.text(xi, 0.48, stars, ha='center', va='bottom',
                    fontsize=16, fontweight='bold', color='black')

    plt.xlim(xlim)
    plt.ylim(ylim)
    xticks = x
    plt.xticks(xticks)
    ax.set_xticklabels(xticks_label, rotation=0, ha='center')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plt.title(roi_abbreviations[roi], pad=15, fontsize=fontsize)
    plt.tight_layout(pad=0.2)

    if if_save == 1:
        plt.savefig(out_path / ('_'.join((roi, process_type, encoding_type)) + '.' + savefig_type), format=savefig_type, dpi=600,
                    bbox_inches='tight')
    # plt.show()  # disabled: figures are saved directly