import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib
matplotlib.use("Agg", force=True)
from matplotlib.ticker import MaxNLocator
from pathlib import Path
import seaborn as sns
from functions.my_plot_setup import set_plot

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#-------------------------------------------------
data_path = BEHAVIOR_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

#-------------------------------------------------
savefig_type = 'png'
figsize = (5, 4)
ylabel = 'Pointing Error'
xlabel = ''
fontsize = 18

#-------------------------------------------------
turn_groups = pd.read_csv(data_path / 'Ave_beha_normal.csv', index_col=False)
fig, ax = set_plot(figsize=figsize, fontsize=fontsize)
sns.violinplot(data=turn_groups, inner="box", palette=["white", "white", "white"],
               linewidth=1.5)
for violin in plt.gca().collections:
    violin.set_edgecolor("black")

plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=4))

plt.xlabel(xlabel,  weight = 'normal')
plt.ylabel(ylabel,  weight = 'normal')

plt.tight_layout(pad =1.5)
# plt.show()  # disabled: figures are saved directly

plt.savefig(out_path /  ('Figure_1c.' + savefig_type) , format=savefig_type, dpi=600, bbox_inches='tight')

