import matplotlib
matplotlib.use("Agg", force=True)


import matplotlib
import matplotlib.patches as patches
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
plt.ioff()
from functions.my_plot_setup import set_plot
from pathlib import Path

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#------------------------------------------
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

custom_colors=['#9DD0C7', '#8AB1D2', '#9180AC', '#D9BDD8', '#E58579']
labels = ['IMF1', 'IMF2', 'IMF3', 'IMF4' , 'IMF5']
#---------------------------------------------------------
savefig_type = 'png'
figsize = (3, 3)
fontsize = 12

#---------------------------------------------------------
fig, ax = set_plot(figsize=figsize, fontsize=fontsize)

#---------------------------------------------------------
legend_elements = []

for i in range(len(custom_colors)):
    cl = custom_colors[i]
    label = labels[i]
    legend_elements.append(patches.Patch(color=cl, label=label))
ax.axis('off')
plt.legend(handles=legend_elements, loc='center', frameon=False)
plt.tight_layout()
# plt.show()  # disabled: figures are saved directly
plt.savefig(out_path/ ('Figure_5a_legend.' + savefig_type) , format=savefig_type, dpi=600, bbox_inches='tight')

