import matplotlib
matplotlib.use("Agg", force=True)


import json
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib
matplotlib.use("Agg", force=True)
import pandas as pd
from pathlib import Path
from functions.my_plot_setup import set_plot


from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#--------------------------------------------------------------------------------------
data_path = DATA_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)
#---------------------------------------------------------------------------------------
process_type = 'turn'
encoding_type = ''

with open(data_path/ ('PLV_results_8_28_' + process_type + encoding_type + '_keyresults.json'), 'r', encoding='utf-8') as f:
    PLV_data = json.load(f)
PLV_data = pd.DataFrame(PLV_data)
print(PLV_data.keys())


data_cuts = PLV_data['data_cut'].unique()
stages = ['1', '2', '3', '4', '5', '6']
filter_bands = ["[1, 4]", "[4, 12]", "[12, 30]", "[30, 70]", "[70, 150]"]



roi1s = ['Temporal_Mid L_Anterior', 'Temporal_Mid L_Posterior',
       'SupraMarginal L', 'Temporal_Sup L','Precentral L',
          'Frontal_Inf_Orb_2 L']
roi1s_labels = ['Ant Temp Mid L','Post Temp Mid L',
                'SupraMarginal L', 'Temp Sup L','Precentral L',
                'Fro_Inf_Orb_2 L']
roi1s_selected = [0, 1, 2, 3, 5]
roi2s = ['Hippocampus L_Anterior', 'Hippocampus L_Posterior',
         'ParaHippocampal L', 'OFCpost L', ]
roi2s_labels = ['Ant Hipp L','Post Hipp L','ParaHipp L','OFCpost L']
roi2s_selected = [0, 1, 2, 3]

#------------------------------------------------
savefig_type = 'png'
figsize = (5+ 4*len(roi2s_selected)/2.54, (3+ 4*len(roi1s_selected))/2.54)
fontsize = 12

if process_type == 'turn':
    xlim = (-0.5, 5.5)
elif process_type == 'corridor':
    xlim = (-0.5, 6.5)
ylim = (0,6,5)

xlabel = 'Junction number'

#-------------------------------------------------------------
plt.close('all')
fig, axs = set_plot(len(roi1s_selected), len(roi2s_selected),figsize=figsize, fontsize=fontsize)

#-------------------------------------------------------------
for i, idx_roi1 in enumerate(roi1s_selected):
    idx_row = i
    roi1 = roi1s[idx_roi1]
    for j, idx_roi2 in enumerate(roi2s_selected):
        idx_column = j
        roi2 = roi2s[idx_roi2]
        sub_title = roi2s_labels[j]
        if len(roi1s_selected) == 1:
            current_ax = axs[idx_column]
        elif len(roi2s_selected) == 1:
            current_ax = axs[idx_row]
        else:
            current_ax = axs[idx_row, idx_column]


        current_data = PLV_data[(PLV_data['brain_area1'] == roi1) & (PLV_data['brain_area2'] == roi2)].copy()
        current_data.loc[:, 'filter_band'] = current_data['filter_band'].astype(str)

        line_data = []
        for band in filter_bands:
            values = []
            for cut in data_cuts:
                value = -np.log10(current_data[(current_data['filter_band'] == band) &
                                          (current_data['data_cut'] == cut)]['pValue2'])
                values.append(value.values[0] if not value.empty else 0)
            line_data.append(values)

        if int(np.max((np.array(line_data)))) > 6.5:
            ylim = (0, int(np.max((np.array(line_data)))) + 1)
        else:
            ylim = (0, 6.5)

        significant_value = -np.log10(0.05)
        significant_values = [significant_value] * len(data_cuts)

        for i, band in enumerate(filter_bands):
            current_ax.plot(stages, line_data[i], marker='o', label=band)

        current_ax.plot(stages, significant_values, linestyle='dashed', color='black', label='p=0.05')

        current_ax.set_xlim(xlim)
        current_ax.set_ylim(ylim)

        if idx_row == len(roi1s_selected) - 1:
            current_ax.set_xlabel(xlabel, fontsize=fontsize + 2, fontweight='bold')
        else:
            current_ax.set_xticks([])

        if idx_row == 0:
            current_ax.set_title(f'{sub_title}', fontsize=fontsize + 2, pad=13, fontweight='bold')

        if idx_column == 0:
            current_ax.text(-2.5, 3, roi1s_labels[idx_roi1],
                                         rotation = 90,
                                          ha='center',
                                          va='center',
                                          fontweight = 'bold'
                                          )
            current_ax.set_ylabel('-log10(p)')

axs[0, len(roi2s_selected) -1 ].legend(fontsize=8, ncol=1, bbox_to_anchor=(1.5, 1))

plt.tight_layout(pad=1, rect=(0.0, 0.0, 0.99, 1))
# plt.show()  # disabled: figures are saved directly
plt.savefig(out_path / ('Sup_Figure_11.' + savefig_type) , format=savefig_type, dpi=600, bbox_inches='tight' )


