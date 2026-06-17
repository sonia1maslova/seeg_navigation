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
from statsmodels.stats.multitest import fdrcorrection

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#------------------------------------------
data_path = DATA_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

#------------------------------------------
process_type = 'corridor'
encoding_type = 'order+end'
#------------------------------------------
data_df = pd.read_json(DATA_DIR / 'Lasso' / 'Temp_Mid_L_corridor_6IMF.json')
roi = 'Temporal_Mid L'

#----------------------------------------------------
with open(DATA_DIR / 'roi_abbreviations.json', 'r', encoding='utf-8') as f:
    roi_abbreviations = json.load(f)

#------------------------------------------
process_label = 'Corridor'
n_stages = 7

stages = list(range(1, 8))
stages = np.array(stages)
stages_label = [str(i) for i in range(1, n_stages)] + ['last']

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

color_this_roi = roi_cmap_df[roi.replace(' ', '_')].values
r_values = []
p_values = []
for i in range(n_stages):
    current_df = data_df[data_df['segi'] == stages[i]]
    r = current_df['Corr'].values.item()
    p = current_df['p-value'].values.item()
    r_values.append(r)
    p_values.append(p)

#---------------------------------------------------------
reject, p_values_corrected = fdrcorrection(p_values, alpha=0.05, method='indep')

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
    plt.savefig(out_path / ('_'.join((roi, process_type, encoding_type, '6imfs')) + '.' + savefig_type), format=savefig_type, dpi=600,
                bbox_inches='tight')
# plt.show()  # disabled: figures are saved directly



