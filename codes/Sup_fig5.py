import matplotlib
matplotlib.use("Agg", force=True)
import pandas as pd
import numpy as np
import matplotlib
from functions.figure_helpers import plot_features_multiple_rois
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
plt.ioff()
from pathlib import Path
import json

from project_paths import DATA_DIR, BEHAVIOR_DIR, SPECTRUM_DIR, AAL3_DIR, get_figure_output_dir
#==========================================================
data_path = DATA_DIR
out_path = get_figure_output_dir(__file__)
out_path.mkdir(exist_ok=True)
#------------------------------------
with open(data_path / 'feature_names.json', 'r', encoding='utf-8') as f:
    feature_names = json.load(f)

#-------------------------------------------
process_type = 'turn'
encoding_type = 'order+end'

#--------------------------------------------
LMM_results_df = pd.read_csv(data_path / 'LMM' /process_type / ('LMM_'+ encoding_type + '.csv'))

mask = ((LMM_results_df['corr_p_corrected'] < 0.05) &
        (LMM_results_df['corr'] > 0)&
        (LMM_results_df['fixed_effect'] > 0)&
        (LMM_results_df['fixed_effect_p'] < 0.05))
signif_lasso_results_df = LMM_results_df[mask]

#-----------------------------------------
lasso_results_with_coefficient_df = pd.read_json(data_path /
                        'Lasso' / process_type/ ('lasso_results_' +encoding_type + '.json'))

lasso_results_with_coefficient_df['segi'].unique()

rois = signif_lasso_results_df['roi'].unique()
rois = rois[rois != 'Temporal_Mid L']

#----------------------------------------------------
with open(DATA_DIR / 'roi_abbreviations.json', 'r', encoding='utf-8') as f:
    roi_abbreviations = json.load(f)
#----------------------------------------
savefig_type = 'png'
custom_colors=['#9DD0C7', '#D9BDD8', '#9180AC', '#8AB1D2', '#E58579'] #  IMF3 IMF5 IMF4 IMF1 IMF2

#-----------------------------------------
coeff_data_list = []
segi_labels = []
rois_label = []
for roi in rois:
    signif_lasso_results_this_roi_df = signif_lasso_results_df[signif_lasso_results_df['roi'] ==roi]
    segments  = signif_lasso_results_this_roi_df['segment'].to_numpy()

    lasso_results_with_coefficient_this_roi_df = (
        lasso_results_with_coefficient_df[lasso_results_with_coefficient_df['loc_name'] == roi])

    mask =lasso_results_with_coefficient_this_roi_df['segi'].isin(segments)
    current_data_df = lasso_results_with_coefficient_this_roi_df[mask]
    current_data_df.reset_index(drop=True, inplace=True)

    #----------------------------------------
    segi_labels_this_roi = [d.replace('turn_', 'Junction') for d in segments]
    print(segi_labels_this_roi)
    segi_labels_this_roi = [d  if 'end' not in d else 'Last junction' for d in segi_labels_this_roi]

    segi_labels = segi_labels + segi_labels_this_roi

    # ---------------------------------------
    rois_label_this_roi = current_data_df['loc_name'].tolist()
    rois_label_this_roi = [roi_abbreviations[d] for d in rois_label_this_roi]
    rois_label = rois_label + rois_label_this_roi

    #---------------------------------------
    coeff_data_this_roi = current_data_df['coefficients']
    for row in range(len(coeff_data_this_roi)):
        coeff_data_list.append(coeff_data_this_roi[row])

#--------------------------------------------
coeff_data = np.array(coeff_data_list)
coeff_data = np.abs(coeff_data)
coeff_df = pd.DataFrame(coeff_data, index=segi_labels, columns=feature_names)
#----------------------------------------------
imf_prefix = ['imf1', 'imf2', 'imf3', 'imf4', 'imf5', ]
custom_colors=['#9DD0C7', '#8AB1D2', '#9180AC', '#D9BDD8', '#E58579']
prefix_color_map = dict(zip(imf_prefix, custom_colors))

#----------------------------------------------
fig = plot_features_multiple_rois(coeff_df, segi_labels, rois_label, '', prefix_color_map)
# plt.show()  # disabled: figures are saved directly
plt.savefig(out_path / ('S5.' + savefig_type), format=savefig_type, dpi=600, bbox_inches='tight')



