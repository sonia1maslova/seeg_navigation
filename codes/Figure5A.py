import matplotlib
matplotlib.use("Agg", force=True)

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.colors as mcolors
import re
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
plt.ioff()
from pathlib import Path
import json

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
def plot_features(coeff_df, stages, title, prefix_color_map):

    pattern = re.compile(r'(.*_)?Higuchi_fractal_dimension')
    coeff_df.columns = [pattern.sub(r'\1HFD', col) for col in coeff_df.columns]

    G = nx.Graph()

    fig, ax = plt.subplots(figsize=(29.7/2.54, 25/2.54))

    fig.suptitle(title, fontsize=12)
    fig.subplots_adjust(top=0.97)

    feature_frequency = coeff_df.apply(lambda x: x > 0).sum(axis=0).to_dict()

    for i, stage in enumerate(stages):

        ax.text(i, len(coeff_df.columns), stage, fontsize=10, ha='center', va='bottom' )

        stage_features = coeff_df.iloc[i,:]

        sorted_features = stage_features.sort_values(ascending=True)

        for j in range(len(sorted_features)):
            feature = sorted_features.index[j]
            value = sorted_features.values[j]

            if value > 0:
                G.add_node((i, j), value=value, feature=feature)

                if i > 0:
                    for k in range(len(coeff_df.columns)):
                        if G.has_node((i-1, k)) and G.nodes[(i-1, k)]['feature'] == feature:
                            G.add_edge((i-1, k), (i, j), weight=feature_frequency[feature])

    pos = {node: (node[0], node[1]) for node in G.nodes()}
    edges = G.edges(data=True)
    weights = [edge[2]['weight'] for edge in edges]
    max_weight = max(weights) if weights else 1  # \u907f\u514d\u9664\u4ee5\u96f6
    alphas = [weight / max_weight for weight in weights]  # \u8ba1\u7b97\u900f\u660e\u5ea6
    for edge, alpha in zip(edges, alphas):
        nx.draw_networkx_edges(G, pos, edgelist=[edge], ax=ax, alpha=alpha, edge_color='#D0D2D4', width=edge[2]['weight']/1.5)

    for i, stage in enumerate(stages):
        stage_features = coeff_df.iloc[i]
        sorted_features = stage_features.sort_values(ascending=True)
        for j in range(len(sorted_features)):
            feature = sorted_features.index[j]
            value = sorted_features.values[j]
            if value > 0:
                prefix = feature.split('_')[0]
                color = prefix_color_map[prefix]
                ax.scatter(i, j, s=value*70+5, c=color, alpha=1)
                ax.text(i, j, feature, fontsize=7, ha='right', va='center')

    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages)
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    return fig

#==========================================================
data_path = DATA_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)

#------------------------------------
process_type= 'corridor'

#------------------------------------
with open(data_path / 'feature_names.json', 'r', encoding='utf-8') as f:
    feature_names = json.load(f)

#-----------------------------------------
results_df = pd.read_json(data_path / 'Lasso' / process_type / 'Ante_Post_splited' / 'lasso_results_order+end.json')

roi = 'Temporal_Mid L_Anterior'


results_TML_Ante_df = results_df[results_df['loc_name'] == roi]
results_TML_Ante_df.reset_index(drop=True, inplace=True)
print(results_TML_Ante_df.columns)
segi_labels = results_TML_Ante_df['segi'].tolist()
segi_labels = [d.replace('corridor_', 'Corridor') for d in segi_labels]
segi_labels[-1] = 'Last corridor'

coeff_data = results_TML_Ante_df['coefficients']
coeff_data_list = []
for row in range(len(coeff_data)):
    coeff_data_list.append(coeff_data[row])

coeff_data = np.array(coeff_data_list)
coeff_data = np.abs(coeff_data)
coeff_df = pd.DataFrame(coeff_data, index=segi_labels, columns=feature_names)
stages = segi_labels

#---------------------------------------------
savefig_type = 'pdf'

imf_prefix = ['imf1', 'imf2', 'imf3', 'imf4', 'imf5', ]
# custom_colors=['#9DD0C7', '#D9BDD8', '#9180AC', '#8AB1D2', '#E58579']
custom_colors=['#9DD0C7', '#8AB1D2', '#9180AC', '#D9BDD8', '#E58579']
# custom_colors=['#FBB9BA', '#D4E5F4', '#C2E3D0', '#F5EBF4', '#F6E4D0']
prefix_color_map = dict(zip(imf_prefix, custom_colors))

fig = plot_features(coeff_df, stages, '', prefix_color_map)
# plt.show()  # disabled: figures are saved directly
plt.savefig(out_path/  ('Figure_5a.' + savefig_type) , format=savefig_type, dpi=600, bbox_inches='tight')

