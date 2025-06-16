from sklearn.decomposition import TruncatedSVD
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import LinearRegression
import numpy as np
import pickle
from sklearn.linear_model import LogisticRegression
from scipy.stats import norm
# Breakpoint on matplotlib.pyplot import --> otherwise debugger crashes as soon as a plot is generated (works fine if not in debugger)
import matplotlib.pyplot as plt
import stepmix


# with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/ADT_CLR','rb') as file:
#         adt_matrix = pickle.load(file)

# adt_matrix = adt_matrix.transpose()
# adt_matrix = adt_matrix[0:7532,0:2]

# with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/GEX_normed','rb') as file:
#         gex_matrix = pickle.load(file)

# gex_matrix = gex_matrix.transpose()

with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/GEX_raw','rb') as file:
        gex_matrix = pickle.load(file)

with open('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/ADT_raw','rb') as file:
        gex_matrix = pickle.load(file)




# print('a')


# smix = stepmix.StepMix(n_components = 3,n_steps = 3, measurement = 'covariate', structural = 'gaussian_full')
# svd = TruncatedSVD(n_components= 300)
# X_svd = svd.fit_transform(gex_matrix)
# smix.fit(X_svd,adt_matrix)



# for i in range(1,6):
#     gmm = GaussianMixture(n_components= i)
#     gmm.fit(X_svd)
#     print(gmm.bic(X_svd))

# gmm1 = GaussianMixture(n_components=2)
# gmm2 = GaussianMixture(n_components=3)
# gmm1.fit(X_svd)
# gmm2.fit(X_svd)

# gmm1_idents = gmm1.predict(X_svd)
# gmm2_idents = gmm2.predict(X_svd)

# adt_svd = TruncatedSVD(n_components= 31)
# adt_matrix = adt_matrix[0:7532,:]
# Y_svd = adt_svd.fit_transform(adt_matrix)
# Y_svd = adt_matrix

# classification_init = np.ones(Y_svd[:,0].shape[0])
# for i in range(Y_svd[:,0].shape[0]):
#         if Y_svd[i,0] <= np.quantile(Y_svd[:,0],0.33):
#                 classification_init[i] = 0
#         elif Y_svd[i,0] <= np.quantile(Y_svd[:,0],0.66):
#                 classification_init[i] = 1
#         else:
#                 classification_init[i] = 2
        

# classifier = LogisticRegression(multi_class = 'multinomial', solver = 'lbfgs', max_iter = 10000)
# classifier.fit(X_svd,classification_init)

# classes = classifier.predict(X_svd)
# lreg1 = LinearRegression()
# lreg2 = LinearRegression()
# lreg3 = LinearRegression()

# for i in range(100):
#  print(i)
#  lreg1.fit(X_svd[classes == 0,:],Y_svd[classes == 0,0:1])
#  lreg2.fit(X_svd[classes == 1,:],Y_svd[classes == 1,0:1])
#  lreg3.fit(X_svd[classes == 2,:],Y_svd[classes == 2,0:1])
#  prediction_1 = lreg1.predict(X_svd)
#  prediction_2 = lreg2.predict(X_svd)
#  prediction_3 = lreg3.predict(X_svd)
#  p1 = np.zeros((Y_svd.shape[0],1))
#  p2 = np.zeros((Y_svd.shape[0],1))
#  p3 = np.zeros((Y_svd.shape[0],1))

#  sd1 = np.std(prediction_1[classes == 0,0])
#  sd2 = np.std(prediction_2[classes == 1,0])
#  sd3 = np.std(prediction_3[classes == 2,0])
#  p1[:,0] =  norm.pdf(Y_svd[:,0],prediction_1[:,0],sd1)
#  p2[:,0] =  norm.pdf(Y_svd[:,0],prediction_2[:,0],sd2)
#  p3[:,0] =  norm.pdf(Y_svd[:,0],prediction_3[:,0],sd3)
# #  lp1 = np.log(p1)
# #  lp2 = np.log(p2)
# #  lp3 = np.log(p3)
#  pp1 = np.prod(p1,axis=1)
#  pp2 = np.prod(p2,axis=1)
#  pp3 = np.prod(p3,axis=1)
#  arr = np.zeros((len(pp1),3))
#  arr[:,0] = pp1
#  arr[:,1] = pp2
#  arr[:,2] = pp3
 
#  classes_previous = classes
#  classes = np.argmax(arr,axis = 1)
#  print(np.bincount(classes))
#  print(sum(classes_previous != classes),'class labels changed')
# #  classifier.fit(X_svd,classes)
# #  classes = classifier.predict(X_svd)
#  predictions_all = np.zeros((len(pp1),1,3))
#  predictions_all[:,:,0] = prediction_1
#  predictions_all[:,:,1] = prediction_2
#  predictions_all[:,:,2] = prediction_3
#  prediction_best = np.zeros((len(pp1),1))
#  for c in range(len(classes)):
#         prediction_best[c,:] = predictions_all[c,:,classes[c]]
#  if i % 10 == 0:
#         plt.scatter(Y_svd[classes == 0,0],prediction_best[classes == 0,0])
#         plt.scatter(Y_svd[classes == 1,0],prediction_best[classes == 1,0])
#         plt.scatter(Y_svd[classes == 2,0],prediction_best[classes == 2,0])
#         ax = plt.gca()
#         ax.set_aspect('equal', adjustable='box')
#         plt.show()
#  if sum(classes_previous != classes) == 0:
#         print('No class labels changed - algorithm converged')
#         break


# lregs = [lreg1,lreg2,lreg3]
# for i in range(100):
#  print(i)
#  for g in range(len(lregs)):
#         lregs[g].fit(X_svd[classes == g,:],Y_svd[classes == g,0:1])
#  predictions = [lreg.predict(X_svd) for]
#  prediction_1 = lreg1.predict(X_svd)
#  prediction_2 = lreg2.predict(X_svd)
#  prediction_3 = lreg3.predict(X_svd)
#  p1 = np.zeros((Y_svd.shape[0],1))
#  p2 = np.zeros((Y_svd.shape[0],1))
#  p3 = np.zeros((Y_svd.shape[0],1))

#  sd1 = np.std(prediction_1[classes == 0,0])
#  sd2 = np.std(prediction_2[classes == 1,0])
#  sd3 = np.std(prediction_3[classes == 2,0])
#  p1[:,0] = (sum(classes == 0)/len(classes))* norm.pdf(Y_svd[:,0],prediction_1[:,0],sd1)
#  p2[:,0] = (sum(classes == 1)/len(classes))*norm.pdf(Y_svd[:,0],prediction_2[:,0],sd2)
#  p3[:,0] = (sum(classes == 2)/len(classes))*norm.pdf(Y_svd[:,0],prediction_3[:,0],sd3)
# #  lp1 = np.log(p1)
# #  lp2 = np.log(p2)
# #  lp3 = np.log(p3)
#  pp1 = np.prod(p1,axis=1)
#  pp2 = np.prod(p2,axis=1)
#  pp3 = np.prod(p3,axis=1)
#  arr = np.zeros((len(pp1),3))
#  arr[:,0] = pp1
#  arr[:,1] = pp2
#  arr[:,2] = pp3
 
#  classes_previous = classes
#  classes = np.argmax(arr,axis = 1)
#  print(np.bincount(classes))
#  print(sum(classes_previous != classes),'class labels changed')
# #  classifier.fit(X_svd,classes)
# #  classes = classifier.predict(X_svd)
#  predictions_all = np.zeros((len(pp1),1,3))
#  predictions_all[:,:,0] = prediction_1
#  predictions_all[:,:,1] = prediction_2
#  predictions_all[:,:,2] = prediction_3
#  prediction_best = np.zeros((len(pp1),1))
#  for c in range(len(classes)):
#         prediction_best[c,:] = predictions_all[c,:,classes[c]]
#  if i % 10 == 0:
#         plt.scatter(Y_svd[classes == 0,0],prediction_best[classes == 0,0])
#         plt.scatter(Y_svd[classes == 1,0],prediction_best[classes == 1,0])
#         plt.scatter(Y_svd[classes == 2,0],prediction_best[classes == 2,0])
#         ax = plt.gca()
#         ax.set_aspect('equal', adjustable='box')
#         plt.show()
#  if sum(classes_previous != classes) == 0:
#         print('No class labels changed - algorithm converged')
#         break


# lreg1 = LinearRegression()
# lreg2 = LinearRegression()
# lreg3 = LinearRegression()
# lreg1.fit(X_svd[classes == 0,:],Y_svd[classes == 0,:])
# lreg2.fit(X_svd[classes == 1,:],Y_svd[classes == 1,:])
# lreg3.fit(X_svd[classes == 2,:],Y_svd[classes == 2,:])
# prediction_1 = lreg1.predict(X_svd)
# prediction_2 = lreg2.predict(X_svd)
# prediction_3 = lreg3.predict(X_svd)
# predictions_all = np.zeros((prediction_1.shape[0],prediction_1.shape[1],3))
# predictions_all[:,:,0] = prediction_1
# predictions_all[:,:,1] = prediction_2
# predictions_all[:,:,2] = prediction_3
# prediction_best = np.zeros((prediction_1.shape[0],prediction_1.shape[1]))
# for c in range(len(classes)):
#         prediction_best[c,:] = predictions_all[c,:,classes[c]]

# # prediction_best_remapped = np.matmul(prediction_best,adt_svd.components_)
# prediction_best_remapped = prediction_best

# for col in range(prediction_best.shape[1]):
#         plt.scatter(Y_svd[classes == 0,col],prediction_best[classes == 0,col])
#         plt.scatter(Y_svd[classes == 1,col],prediction_best[classes == 1,col])
#         plt.scatter(Y_svd[classes == 2,col],prediction_best[classes == 2,col])
#         plt.show()

# for col in range(prediction_best_remapped.shape[1]):

#         plt.scatter(adt_matrix[classes == 0,col],prediction_best_remapped[classes == 0,col])
#         plt.scatter(adt_matrix[classes == 1,col],prediction_best_remapped[classes == 1,col])
#         plt.scatter(adt_matrix[classes == 2,col],prediction_best_remapped[classes == 2,col])
#         plt.show()

# print('a')