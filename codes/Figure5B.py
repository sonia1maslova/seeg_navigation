import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.colors as mcolors
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import seaborn as sns
from functions.figure_helpers import dice_similarity
import pandas as pd
from pathlib import Path
import json
from functions.my_plot_setup import set_plot


from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#------------------------------------------
data_path = DATA_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)
#------------------------------------
with open(data_path / 'feature_names.json', 'r', encoding='utf-8') as f:
    feature_names = json.load(f)

#-----------------------------------------
results_df = pd.read_json(data_path / 'Lasso' / 'corridor' / 'Ante_Post_splited' / 'lasso_results_order+end.json')

roi = 'Temporal_Mid L_Anterior'

results_TML_Ante_df = results_df[results_df['loc_name'] == roi]
results_TML_Ante_df.reset_index(drop=True, inplace=True)
print(results_TML_Ante_df.columns)
segi_labels = results_TML_Ante_df['segi'].tolist()
segi_labels = [d.replace('corridor_', 'C') for d in segi_labels]
segi_labels[-1] = 'end'

coeff_data = results_TML_Ante_df['coefficients']
coeff_data_list = []
for row in range(len(coeff_data)):
    coeff_data_list.append(coeff_data[row])

coeff_data = np.array(coeff_data_list)
coeff_data[coeff_data!= 0] = 1

dice_matrix = np.full((len(coeff_data), len(coeff_data)), np.nan)
for idx_1 in range(len(coeff_data)):
    for idx_2 in range(len(coeff_data)):
        coeff_1 = coeff_data[idx_1]
        coeff_2 = coeff_data[idx_2]
        dice_matrix[idx_1, idx_2] = dice_similarity(coeff_1, coeff_2)

np.fill_diagonal(dice_matrix, 0)
dice_df = pd.DataFrame(dice_matrix, index=segi_labels, columns=segi_labels)

mask = np.triu(np.ones_like(dice_df, dtype=bool))

#---------------------------------------------------------
savefig_type = 'pdf'
figsize = (7/2.54, 6/2.54)
fontsize = 12

cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ['#FFFFFF', '#C25A5A'])
vmin = 0
vmax = 0.5

fig, ax = set_plot(figsize=figsize, fontsize=fontsize)
cax = fig.add_axes([0.7, 0.3, 0.025, 0.4])  # [left, bottom, width, height]
sns.heatmap(dice_df.apply(pd.to_numeric, errors='coerce'),
            annot=True, cmap=cmap, yticklabels=segi_labels, xticklabels=segi_labels,
            fmt=".2f", annot_kws={"size": 8}, mask=mask,
            ax=ax, cbar_ax=cax, square=True, vmin=vmin, vmax=vmax)

ax.set_xticklabels(ax.get_xticklabels())
ax.tick_params(axis='both', which='both',
               direction='out', length=0, width=0,
               bottom=True, top=False, left=True, right=False)
for spine in cax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.6)
    spine.set_color('black')

cax.tick_params(
    axis='y',                    # 'x', 'y', 'both'
    which='major',               # 'major', 'minor', 'both'
    direction='out',             # 'in', 'out', 'inout'
    length=2,                    # 刻度长度
    width=0.5,                     # 刻度宽度
    color='black',               # 刻度颜色
    labelsize=9,                # 标签字体大小
    labelcolor='black',          # 标签颜色
    pad=3,                       # 标签与刻度距离
    rotation=0                   # 标签旋转角度
)
cax.set_yticks([0, 0.2, 0.4])
cax.set_yticklabels(['0', '0.2', '0.4'])
cax.set_ylabel('Dice Coefficient ', fontsize=10, rotation=270, labelpad=10)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# plt.show()  # disabled: figures are saved directly
plt.savefig(out_path / ('Figure_5b.' + savefig_type) , format=savefig_type, dpi=600, bbox_inches='tight' )

