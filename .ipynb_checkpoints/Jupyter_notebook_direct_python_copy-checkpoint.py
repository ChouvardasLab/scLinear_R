from sklearnex import patch_sklearn
# Reimplementations of sklearn functions with optimization for speed
patch_sklearn(global_patch=True)
from sklearn.decomposition import TruncatedSVD
# from sklearn.mixture import GaussianMixture
# from sklearn.linear_model import LinearRegression
# from sklearn.model_selection import GridSearchCV
import numpy as np
import pickle
# from sklearn.linear_model import LogisticRegression
# from scipy.stats import norm
# from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
# import stepmix
from sklearn import svm

# with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/Tabula_sapiens_classifier_input/gene_expression_tabsap','rb') as file:
#         gex_matrix = pickle.load(file)
#         gex_matrix = gex_matrix.transpose()

with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/Tabula_sapiens_classifier_input/adt_predictions_for_tabsap','rb') as file:
        adt_predictions = pickle.load(file)

with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/Tabula_sapiens_classifier_input/tabsap_annotations','rb') as file:
        annotations = pickle.load(file)

# with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/Tabula_sapiens_classifier_input/tsvd_projector_from_3tmodel','rb') as file:
#         tsvd_v = pickle.load(file)

# with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/Tabula_sapiens_classifier_input/gex_filtered_for_tsvd','rb') as file:
#         gex_filtered = pickle.load(file)

# Keep csc since it uses only about 1/10th of memory compared to array
with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/Tabula_sapiens_classifier_input/gene_expression_tabsap_normed_filtered','rb') as file:
        gex_normed_filtered = pickle.load(file)

csc_minor_stored_values = np.bincount(gex_normed_filtered.indices, minlength=gex_normed_filtered.shape[0])
keep_rows = np.where(csc_minor_stored_values)[0]
gex_normed_filtered = gex_normed_filtered[keep_rows,]
gex_normed_filtered = gex_normed_filtered.transpose()
# No empty cols / droplets --> probably already filtered before

gex_projector = TruncatedSVD(n_components = 300)
gex_projected = gex_projector.fit_transform(gex_normed_filtered)