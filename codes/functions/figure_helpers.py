"""Minimal helper functions used by the shared plotting scripts."""

from __future__ import annotations

import re
import numpy as np


def dice_similarity(x, y):
    """Calculate Dice coefficient between two binary arrays."""
    x = np.asarray(x)
    y = np.asarray(y)
    intersection = np.sum(np.logical_and(x, y))
    x_sum = np.sum(x)
    y_sum = np.sum(y)
    return (2 * intersection) / (x_sum + y_sum) if (x_sum + y_sum) != 0 else 0


def plot_features_multiple_rois(coeff_df, stages, rois, title, prefix_color_map):
    """Plot selected Lasso features across multiple ROIs.

    This function is a trimmed version of the original helper and keeps only
    the logic used by Sup_fig4.py and Sup_fig5.py.
    """
    import networkx as nx
    import matplotlib.pyplot as plt

    n_stage = len(stages)
    figsize = ((5 * n_stage + 2) / 2.54, 25 / 2.54)

    pattern = re.compile(r'(.*_)?Higuchi_fractal_dimension')
    coeff_df.columns = [pattern.sub(r'\1HFD', col) for col in coeff_df.columns]

    x_pos = []
    unique_rois = []
    unique_rois_x_po = []
    x_pos_this_roi = []
    for i, roi in enumerate(rois):
        if i >= 1:
            if rois[i] == rois[i - 1]:
                x_po = x_pos[-1] + 1
            else:
                x_po = x_pos[-1] + 1.3
        else:
            x_po = 0
        x_pos.append(x_po)

        if roi not in unique_rois:
            if x_pos_this_roi:
                unique_rois_x_po.append(sum(x_pos_this_roi) / len(x_pos_this_roi))
            x_pos_this_roi = []
            unique_rois.append(roi)
        x_pos_this_roi.append(x_pos[-1])
        if i == len(rois) - 1:
            unique_rois_x_po.append(sum(x_pos_this_roi) / len(x_pos_this_roi))

    G = nx.Graph()

    fig, ax = plt.subplots(figsize=figsize)

    fig.suptitle(title, fontsize=12)
    fig.subplots_adjust(top=0.95)

    feature_frequency = coeff_df.apply(lambda x: x > 0).sum(axis=0).to_dict()

    for i, stage in enumerate(stages):
        ax.text(x_pos[i], len(coeff_df.columns) + 0.2, stage,
                fontsize=10, ha='center', va='bottom')

        stage_features = coeff_df.iloc[i, :]
        sorted_features = stage_features.sort_values(ascending=True)

        for j in range(len(sorted_features)):
            feature = sorted_features.index[j]
            value = sorted_features.values[j]

            if value > 0:
                G.add_node((i, j), pos=(x_pos[i], j), value=value, feature=feature)
                if (i > 0) & (rois[i] == rois[i - 1]):
                    for k in range(len(coeff_df.columns)):
                        if G.has_node((i - 1, k)) and G.nodes[(i - 1, k)]['feature'] == feature:
                            G.add_edge((i - 1, k), (i, j), weight=feature_frequency[feature])

    pos = {}
    for node in G.nodes():
        i, j = node
        x_coord = x_pos[i] if i < len(x_pos) else i
        pos[node] = (x_coord, j)

    edges = G.edges(data=True)
    weights = [edge[2]['weight'] for edge in edges]
    max_weight = max(weights) if weights else 1
    alphas = [weight / max_weight for weight in weights]
    for edge, alpha in zip(edges, alphas):
        nx.draw_networkx_edges(
            G, pos, edgelist=[edge], ax=ax, alpha=alpha,
            edge_color='#D0D2D4', width=edge[2]['weight'] / 1.5
        )

    for i, stage in enumerate(stages):
        stage_features = coeff_df.iloc[i]
        sorted_features = stage_features.sort_values(ascending=True)
        for j in range(len(sorted_features)):
            feature = sorted_features.index[j]
            value = sorted_features.values[j]
            if value > 0:
                prefix = feature.split('_')[0]
                color = prefix_color_map[prefix]
                ax.scatter(x_pos[i], j, s=value * 70 + 5, c=color, alpha=1)
                ax.text(x_pos[i], j, feature, fontsize=7, ha='right', va='center')

    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages)
    ax.set_yticks([])

    for idx, unique_roi in enumerate(unique_rois):
        ax.text(unique_rois_x_po[idx], len(coeff_df.columns) + 1.5, unique_roi,
                fontsize=12, ha='center', va='bottom', weight='bold')

    for spine in ax.spines.values():
        spine.set_visible(False)

    return fig
