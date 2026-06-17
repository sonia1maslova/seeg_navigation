import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.ticker import MaxNLocator
from pathlib import Path
import pandas as pd
from functions.my_plot_setup import set_plot
import numpy as np
from matplotlib.ticker import MaxNLocator


from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#----------------------------------------------
data_path = SPECTRUM_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

#----------------------------------------------
roi = 'Hippocampus_L_Anterior'
spectrum_data = pd.read_json(data_path / (roi + '_spectrum_data.json'))
print(spectrum_data.columns)

spectrum_data_copy = spectrum_data.copy()
spectrum_data_good = spectrum_data_copy[spectrum_data_copy['performance_group'] == 'good']
spectrum_data_bad = spectrum_data_copy[spectrum_data_copy['performance_group'] == 'bad']
spectrum_data_normal = spectrum_data_copy[spectrum_data_copy['performance_group'] == 'normal']

Avg_power = np.mean(np.vstack(spectrum_data_copy['power'].values), axis=0)
Avg_power_good = np.mean(np.vstack(spectrum_data_good['power'].values), axis=0)
Avg_power_bad = np.mean(np.vstack(spectrum_data_bad['power'].values), axis=0)
Avg_power_normal = np.mean(np.vstack(spectrum_data_normal['power'].values), axis=0)

all_frequency = spectrum_data_good.loc[0, 'frequency']
frequency_band = list(range(1, 31))
idx_freq = []
for idx, freq in enumerate(all_frequency):
    if freq in frequency_band:
        idx_freq.append(idx)


#-------------------------------------------
savefig_type = 'png'
figsize = (6, 4)
fontsize = 16

individual_linewidth = 0.3
Ave_linewidth = 2.5
alpha = 0.05

cmap =  np.array([[247,70,89], [35,140,240], [30,166,70],
                        [206,170,208], [251,140,98], [54, 52, 242]])/ 255

xlim = [1, 30]

xlabel = 'Frequency'
ylabel = 'Power Spectral Density'
#-------------------------------------------
fig, ax = set_plot(figsize=figsize, fontsize=fontsize, axis_linewidth=1.2)

for idx_row, row in spectrum_data_copy.iterrows():
    x = frequency_band
    y = [row['power'][idx] for idx in idx_freq]

    cl = cmap[1]
    plt.plot(x, y, color=cl, alpha=alpha, linewidth=individual_linewidth)

y_Ave_power = Avg_power[idx_freq]
x = frequency_band
plt.plot(x, y_Ave_power, color = cmap[1], alpha = 1, linewidth = Ave_linewidth,
            label= '')

plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=6))
plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=5))
plt.legend(fontsize=fontsize-4)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim(xlim)
plt.tight_layout(pad=1)
# plt.plot(x, y_Ave_power_normal, color = cmap[2], alpha = 1, linewidth = Ave_linewidth )
# plt.show()  # disabled: figures are saved directly
plt.savefig(out_path / ('Sup_Figure_3.' + savefig_type) , format=savefig_type, dpi=600, bbox_inches='tight' )



