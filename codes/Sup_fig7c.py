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
lasso_data_path = data_path / 'LMM' /process_type
encoding_type = 'order_only'
lasso_results_df = pd.read_csv(lasso_data_path/ ('LMM_'+ encoding_type + '.csv'))
current_rois = lasso_results_df['roi'].unique().tolist()

#------------------------------------------
with open(DATA_DIR / 'LMM' / process_type / ('sig_dict_' + process_type + '.json'), 'r', encoding='utf-8') as f:
    sig_dict = json.load(f)
all_rois_sig = list(sig_dict.keys())

rois = list(set(all_rois_sig) & set(current_rois))
#----------------------------------------------------
with open(DATA_DIR / 'roi_abbreviations.json', 'r', encoding='utf-8') as f:
    roi_abbreviations = json.load(f)

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
roi_cmap_df = pd.read_json(AAL3_DIR / 'roi_cmap.json')

#---------------------------------------------------------
if_save = 1
savefig_type = 'png'
figsize = (3.2, 2.4)

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
for roi in rois:
    color_this_roi = roi_cmap_df[roi.replace(' ', '_')].values
    data_this_roi_df = lasso_results_df[lasso_results_df['roi'] == roi]
    r_values = []
    p_values_corrected = []
    for i in range(n_stages):
        current_df = data_this_roi_df[data_this_roi_df['segment'] == stages[i]]
        r = current_df['corr'].values.item()
        p_corrected = current_df['corr_p_corrected'].values.item()
        r_values.append(r)
        p_values_corrected.append(p_corrected)

    #---------------------------------------------------------
    fig, ax = set_plot(figsize=figsize, fontsize=fontsize)
    x = list(range(1, n_stages+1))
    y = r_values
    colors = [color_this_roi if (p < 0.05 and r > 0) else [0.7, 0.7, 0.7]
              for r, p in zip(r_values, p_values_corrected)]

    ax.bar(x, y, width=bin_width, color=colors, alpha=1)

    for i, (xi, yi, r, p) in enumerate(zip(x, y, r_values, p_values_corrected)):
        if r > 0:
            if p < 0.001:
                stars = '***'
            elif p < 0.01:
                stars = '**'
            elif p < 0.05:
                stars = '*'
            else:
                stars = ''

            if stars:
                star_y = 0.48
                ax.text(xi, star_y, stars, ha='center', va='bottom',
                        fontsize=16, fontweight='bold', color='black')

    plt.xlim(xlim)
    plt.ylim(ylim)
    xticks = x
    plt.xticks(xticks)
    ax.set_xticklabels(xticks_label, rotation=0, ha='center')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plt.title(roi_abbreviations[roi], pad=15, fontsize = fontsize)
    plt.tight_layout(pad=0.2)

    if if_save == 1:
        plt.savefig(out_path / ('_'.join((roi, process_type, encoding_type)) + '.' + savefig_type), format=savefig_type, dpi=600,
                    bbox_inches='tight')
    # plt.show()  # disabled: figures are saved directly



