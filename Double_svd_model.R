kidney1_seurat <- Seurat::Read10X_h5('C:/Users/mz24b548/Downloads/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Kidney_Cancer1_BC1_AB1_count_sample_filtered_feature_bc_matrix.h5')
devtools::load_all('R/')
library(ggplot2)
plot_comparison <- function(adt1,adt2,xaxis = 'ADT Values 1', yaxis = 'ADT Values 2',title = '',groupings = rep(1,length(adt1))){
  plot_dat <- data.frame(cbind(adt1,adt2,as.factor(groupings)))
  plot_dat[,3] <- as.factor(plot_dat[,3])
  colnames(plot_dat) <- c('ADT1','ADT2','group')
  comp_plot <- ggplot(plot_dat,aes(x=ADT1,y=ADT2,color = group))+
    geom_point()+
    geom_abline(intercept=0,slope=1,color = "red")+
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab(xaxis) + ylab(yaxis) +
    theme_minimal() +
    coord_equal() +
    #If x axis = true values and y = prediction --> predictions have much bigger outlier range --> set limits based on true range
    ylim(c(quantile(adt1,0.02),max(quantile(adt1,0.98),0))) +
    scale_color_discrete(name = "Cell types", labels = levels(as.factor(groupings))) +
    ggtitle(title) + theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "black"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = "black"),
      plot.background = element_rect(fill = "white")

    )
  return(comp_plot)
}


kidney1_gex <- kidney1_seurat$`Gene Expression`
kidney1_adt <- kidney1_seurat$`Antibody Capture`

combi_seurat <- Seurat::CreateSeuratObject(counts = kidney1_gex, assay = 'RNA')
adt_seurat <- Seurat::CreateAssayObject(counts = kidney1_adt)

combi_seurat[["ADT"]] <- adt_seurat
prepped <- prepare_data(combi_seurat)

prepped_gex <- Seurat::GetAssayData(prepped[['RNA']], layer = 'counts')
prepped_adt <- Seurat::GetAssayData(prepped[['ADT']], layer = 'counts')

y <- t(as.matrix(prepped_adt))


# Which ADT counts are 0 --> clr doesn't add pseudocount but instead leaves them as 0 which introduces potentially even a bigger error
# Use no log normalization for gene expression BUT potentially for the TSVD projected values --> they should almost never be exactly
# 0 since they are a linear combination of many different genes
# Can't simply drop 0 values after CLR since we need the value for tSVD IF we go with ADT tsvd / pca
# Calculate geometric mean for each cell, calculate CLR for all non-zero counts as if no pseudocount added, replace CLR of zero-count with CLR of values as if it had been 1
# gmeans_per_cell <- vector(mode = "list", length = ncol(prepped_adt))
# for(i in 1:ncol(prepped_adt)){
#   gmeans_per_cell[i] <- exp(mean(log(prepped_adt[prepped_adt[,i] != 0,i])))
# }
# prepped_adt <- Matrix::Matrix((t(compositions::clr(t(as.matrix(prepped_adt))))),sparse = TRUE)

# to_replace <- which(as.vector(prepped_adt_no_clr) == 0)
# for(entry in to_replace){
#   col <- floor((entry-1)/nrow(prepped_adt)) + 1
#   row <- ((entry-1) %% nrow(prepped_adt)) + 1
#   prepped_adt[row,col] <- log(1/gmeans_per_cell[[col]])
# }


# mu = colMeans(t(prepped_adt))
# adt_pca <- prcomp(t(prepped_adt))
# # Reconstruction if needed
# unscaled = adt_pca$x %*% t(adt_pca$rotation)
# reconstructed = scale(unscaled, center = -mu, scale = FALSE)
#
# adt_pca_centered <- adt_pca$x
# adt_pca_df <- data.frame(as.factor(prepped$cell_type),adt_pca_centered)
# adt_pca_df[,1] <- as.factor(adt_pca_df[,1])
# colnames(adt_pca_df)[1] <- 'cell_type'
# ggplot(adt_pca_df,aes(x = `PC1`, color = cell_type)) + geom_histogram(position = 'identity', fill = NA) +
#   scale_color_discrete(name = "", labels = levels(as.factor(prepped$cell_type)))
#
# adt_clr_centered <- as.matrix(t(prepped_adt))
# adt_clr_df <- data.frame(as.factor(prepped$cell_type),adt_clr_centered)
# adt_clr_df[,1] <- as.factor(adt_clr_df[,1])
# colnames(adt_clr_df)[1] <- 'cell_type'
#
# adt_kec_pca <- prcomp(t(as.matrix(prepped_adt)[,prepped$cell_type == 'Kidney epithelial cell']))
# regular_normed <- t(as.matrix(gexp_normalize(as.matrix(prepped_gex)[,prepped$cell_type == 'Kidney epithelial cell'])))

# Use un-clred ADT data
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = prepped_gex))
clusters <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, clusters=clusters)
gex_normed_no_log <- scuttle::normalizeCounts(sce,center.size.factors = FALSE, log = FALSE)

all_predis_all_cells <- matrix(NA, nrow = 32,ncol = 1)

for(cell_t in unique(prepped$cell_type)){
dreg_y <- t(as.matrix(prepped_adt)[,prepped$cell_type == cell_t])
#Pseudocount --> 0.99 instead of one so it can later be more easily distinguished from true 1 count
dreg_y[dreg_y == 0] <- 0.99
# dreg_y[dreg_y == 0.99] <- 0
dreg_y_scaled <- dreg_y / rowSums(dreg_y)
train_indices <- 1:floor(nrow(dreg_y_scaled)*2/3)
test_indices <- (floor(nrow(dreg_y_scaled)*2/3)+1):nrow(dreg_y_scaled)
gex_tsvd <- sparsesvd::sparsesvd(Matrix::Matrix(t(as.matrix(gex_normed_no_log)[,prepped$cell_type == cell_t]),sparse = TRUE),min(300,length(test_indices)))
gex_projected <- t(as.matrix(gex_normed_no_log)[,prepped$cell_type == cell_t]) %*% gex_tsvd$v
c_ilr <- Matrix::Matrix((t(compositions::ilr(as.matrix(dreg_y_scaled)))),sparse = TRUE)

predi_all_this_cell <- matrix(nrow = nrow(c_ilr),ncol = length(test_indices)  )

for(i in 1:nrow(c_ilr)){
  print(i)
  lmod <- lm(as.matrix(c_ilr)[i,train_indices] ~ gex_projected[train_indices,])
  predi_all_this_cell[i,] <- lmod$coefficients[1]+lmod$coefficients[2:length(lmod$coefficients)] %*% t(gex_projected[test_indices,])
}
remapped <- Matrix::Matrix((t(compositions::ilrInv(t(as.matrix(predi_all_this_cell))))),sparse = TRUE)
cell_names <- rownames(dreg_y_scaled)[test_indices]
colnames(remapped) <- cell_names
all_predis_all_cells <- cbind(all_predis_all_cells,remapped)
}
all_predis_all_cells <- as.matrix(all_predis_all_cells)[,-1]
which_cells <- colnames(all_predis_all_cells)
reference <- t(as.matrix(prepped_adt)[,which_cells])
reference[reference == 0] <- 0.99
reference <- reference / rowSums(reference)
for(i in 1:ncol(reference)){
plot_comparison(log(reference[,i]),log(all_predis_all_cells[i,]), xaxis = 'Log(true fraction)',
                yaxis = 'Log(Prediction)', groupings = prepped$cell_type[which_cells],title = colnames(reference)[i])
ggsave(paste0('Prediction_evaluation_plots/',colnames(reference)[i],'_prediction.png'))
  }


predi_all <- matrix(nrow = nrow(adt_kec_pca$x),ncol = ncol(adt_kec_pca$x))




for(i in 1:ncol(adt_kec_pca$x)){
  print(i)
  mod <- cv.glmnet(regular_normed,adt_kec_pca$x[,i], alpha = 0)
  best_lambda <- mod$lambda.min
  best_model <- glmnet(regular_normed, adt_kec_pca$x[,i], alpha = 0, lambda = best_lambda)
  predi_all[,i] <- predict.glmnet(best_model,regular_normed)
}

predi_all_remapped <- predi_all %*% t(adt_kec_pca$rotation)
mu <- colMeans(t(as.matrix(prepped_adt)[,prepped$cell_type == 'Kidney epithelial cell']))
predi_all_remapped = scale(predi_all_remapped, center = -mu, scale = FALSE)
comparison <- t(as.matrix(prepped_adt)[,prepped$cell_type == 'Kidney epithelial cell'])

mod <- cv.glmnet(regular_normed,adt_kec_pca$x[,1], alpha = 0)
best_lambda <- mod$lambda.min
best_model <- glmnet(regular_normed, adt_kec_pca$x[,1], alpha = 0, lambda = best_lambda)
predi <- predict.glmnet(best_model,regular_normed)




adt_df <- data.frame(prepped_adt)
adt_df <- data.frame(cbind(as.factor(prepped$cell_type),t(adt_df)))
# adt_df <- data.frame(cbind(prepped$seurat_clusters,t(adt_df)))
adt_df[,1] <- as.factor(adt_df[,1])
colnames(adt_df)[1] <- 'cell_type'

library(ggplot2)
ggplot(adt_df,aes(x = `CD3.CD3E`, y = `CD4.CD4`, color = cell_type)) + geom_point() +
  scale_fill_discrete(name = "", labels = levels(as.factor(prepped$cell_type)))

kidney1_adt_mat <- as.matrix(kidney1_adt)

#Checked --> representation of original is very accurate (remapping = pretty much perfect fit. Makes sense since almost FULL svd)
adt_svd <- sparsesvd::sparsesvd(Matrix::t(prepped_adt),nrow(prepped_adt)-1)
adt_projected <- Matrix::t(prepped_adt) %*% adt_svd$v
adt_p_df <- data.frame(cbind(as.factor(prepped$cell_type),adt_projected))
adt_p_df[,1] <- as.factor(adt_p_df[,1])
ggplot(adt_p_df,aes(x = `X2`, y = `X3`, color = X1)) + geom_point() +
  scale_fill_discrete(name = "", labels = levels(as.factor(prepped$cell_type)))

ggplot(adt_p_df,aes(x = `X2`, color = X1)) + geom_histogram(position = 'identity', fill = NA) +
  scale_color_discrete(name = "", labels = levels(as.factor(prepped$cell_type)))

ggplot(adt_df,aes(x = `CD3.CD3E`, color = cell_type)) + geom_histogram(position = 'identity', fill = NA) +
  scale_color_discrete(name = "", labels = levels(as.factor(prepped$cell_type)))


sce <- SingleCellExperiment::SingleCellExperiment(list(counts = prepped_gex))
clusters <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, clusters=clusters)
gex_normed_no_log <- scuttle::normalizeCounts(sce,center.size.factors = FALSE, log = FALSE)
#Filtering still seems to have missed a few multiplexes or otherwise way to high rna counts cells

regular_normed <- gexp_normalize(prepped_gex)

#try_all --> performs about the same as at the very beginning (ok clustering, bad actual prediction)
# library(glmnet)
# model <- cv.glmnet(t(regular_normed), adt_projected[,1], alpha = 0)
# best_lambda <- model$lambda.min
# best_model <- glmnet(t(regular_normed), adt_projected[,1], alpha = 0, lambda = best_lambda)

# Gex log or no log --> not a massive difference i'd say --> slightly worse for no log
gex_kex <- as.matrix(gex_normed_no_log)[,prepped$cell_type == 'Kidney epithelial cell']
# gex_kec <- as.matrix(regular_normed)[,prepped$cell_type == 'Kidney epithelial cell']
adt_kec <- as.matrix(adt_projected)[prepped$cell_type == 'Kidney epithelial cell',]
# adt_kec <- t(as.matrix(prepped_adt)[,prepped$cell_type == 'Kidney epithelial cell'])
# reference  <- t(as.matrix(prepped_adt)[,prepped$cell_type == 'Kidney epithelial cell'])

library(glmnet)

# row = cells, cols = components / adt
predi_all <- matrix(nrow = nrow(adt_kec),ncol = ncol(adt_kec))

for(i in 1:ncol(adt_kec)){
print(i)
model <- cv.glmnet(t(gex_kec), adt_kec[,i], alpha = 0)
best_lambda <- model$lambda.min
best_model <- glmnet(t(gex_kec), adt_kec[,i], alpha = 0, lambda = best_lambda)
predi_all[,i] <- predict.glmnet(best_model,t(gex_kec))
}

predi_remapped <- predi_all %*% t(adt_svd$v)


plot_comparison(adt_kec[,1],predict.glmnet(best_model,t(gex_kec)))



#Remove some ADTs and turn into 'composition' to check if dependence on these ADTs now disappears / becomes lower
# Here removing CD3, CD8a, CD14, CD127_IL7R, CD10_MME
# kidney1_adt <- kidney1_adt_mat[-c(1,3,5,32,9),]
# kidney1_adt_comp <- apply(kidney1_adt,2,function(x){x/sum(x)})
# kidney1_adt <- kidney1_adt_comp

kidney1_adt <- as.matrix(t(compositions::clr(t(as.matrix(kidney1_adt)))))
# kidney1_adt <- t(kidney1_adt)
# kidney1_adt_projector <- svd(kidney1_adt,nu=nrow(kidney1_adt),nv = ncol(kidney1_adt))
# kidney1_adt_projection <- kidney1_adt %*% kidney1_adt_projector$v
# kidney1_adt <- t(kidney1_adt_projection[,-32])

#Just to confirm that the operation is invertible
#kidney1_remapped <- kidney1_adt_projection %*% t(kidney1_adt_projector$v)

indx <- sample(1:ncol(kidney1_seurat$`Gene Expression`), size = ncol(kidney1_seurat$`Gene Expression`), replace = FALSE)


# should have format #features x #cells
# Here the ADT names are Protein_GeneCodingForProtein OR Mouse_

kidney1_gex_train <- kidney1_gex[,indx[1:floor(0.75*length(indx))]]
kidney1_gex_test <- kidney1_gex[,indx[(floor(0.75*length(indx))+1):length(indx)]]


kidney1_adt_train <- kidney1_adt[,indx[1:floor(0.75*length(indx))]]
kidney1_adt_test <- kidney1_adt[,indx[(floor(0.75*length(indx))+1):length(indx)]]


# kidney1_predictor <- fit_predictor(kidney1_gex_train,kidney1_adt_train,kidney1_gex_test,normalize_adt = FALSE)
# kidney1_prediction <- predict(kidney1_predictor,kidney1_gex_test)
#
#
# # for(name in rownames(kidney1_adt_test)){
# #   plot_comparison(kidney1_adt_test[name,],prediction_remapped[name,])
# #   ggplot2::ggsave(paste0('Prediction_evaluation_plots/',name,'.png'))
# # }
# prediction_remapped <- kidney1_prediction %*% t(kidney1_adt_projector$v[,-32])
# test_remapped <- t(kidney1_adt_test) %*% t(kidney1_adt_projector$v[,-32])
# errors <- prediction_remapped - test_remapped
# err_sq <- errors^2
#
# p_adt <- subset(prediction_remapped,features = which(rownames(prediction_remapped) %in% rownames(test_remapped)) )
# t_adt <- subset(test_remapped,features = which(rownames(test_remapped) %in% rownames(prediction_remapped)) )
#
# # t_adt <- t_adt[match(rownames(p_adt),rownames(t_adt)),match(colnames(p_adt),colnames(t_adt))]
#
# # Not sure how to interpret the means of the single model metrics but for now just replicating the behaviour of python code
# rmse <- mean(sqrt(colSums(err_sq)/nrow(err_sq)))
# # Pearson calculated in a really strange way in python --> why would we average across cells and not across models (if at all?)
# # Here reporting the correlation coefficients of each ADT prediction separately
# pearson <- diag(cor(t_adt,p_adt,method = "pearson"))
# spearman <- diag(cor(t_adt,p_adt,method = "spearman"))
# # Original behaviour:
# mean_pearson <- mean(unlist(lapply(1:nrow(t_adt),function(i) cor(t_adt[i,],p_adt[i,],method = 'pearson'))))
# mean_spearman <- mean(unlist(lapply(1:nrow(t_adt),function(i) cor(t_adt[i,],p_adt[i,],method = 'spearman'))))
# plot(test_remapped[,1],errors[,1])
#

kidney1_gex_train_normed <- gexp_normalize(kidney1_gex_train)
kidney1_gex_train_normed <- as.matrix(kidney1_gex_train_normed)
# kidney1_gex_train_normed <- kidney1_gex_train_normed[kidney1_predictor$genes_considered,]
# kidney1_gex_train_projected <- t(kidney1_gex_train_normed) %*% kidney1_predictor$tsvd_v
#
# kidney1_gex_test_normed <- gexp_normalize(kidney1_gex_test)
# kidney1_gex_test_normed <- as.matrix(kidney1_gex_test_normed)
# kidney1_gex_test_normed <- kidney1_gex_test_normed[kidney1_predictor$genes_considered,]
# kidney1_gex_test_projected <- t(kidney1_gex_test_normed) %*% kidney1_predictor$tsvd_v
#
# X1 <- kidney1_gex_train_projected
# Y <- kidney1_adt_train
# # # mod1 <- stepFlexmix(kidney1_adt_train[1,] ~ kidney1_gex_train_projected,control = list(verbose = 0), k = 1:5,nrep = 5)
# library(flexmix)
# fm_models <- flexmix(~X1, k = 2, model = list(
#   FLXMRglm(Y[1,]~.),
#   FLXMRglm(Y[2,]~.),
#   FLXMRglm(Y[3,]~.),
#   FLXMRglm(Y[4,]~.),
#   FLXMRglm(Y[5,]~.),
#   FLXMRglm(Y[6,]~.),
#   FLXMRglm(Y[7,]~.)
# ))
#
#
# fm_prediction <- flexmix::predict(fm_models)
# fm_predi_dim1 <- matrix(nrow = nrow(X1),ncol = 1)
# clusters <- flexmix::clusters(fm_models)
#
# comp <- 7
# for(i in 1:nrow(X1)){
#   fm_predi_dim1[i,1] <- fm_prediction[[clusters[i]]][i,comp]
# }
#
#
# coeff_mat <- matrix(nrow = length(kidney1_predictor$lm_coefficients),ncol = length(kidney1_predictor$lm_coefficients[[1]]))
# for(i in 1:length(kidney1_predictor$lm_coefficients)){
#   coeff_mat[i,] <- kidney1_predictor$lm_coefficients[[i]]
# }
#
# clusters_f <- as.factor(clusters)
# rfc = randomForest::randomForest(X1,clusters_f) #poor fit even for train itself
# rft <- randomForest::randomForest(X1,Y[1,]) #good fit for train but extremely slow --> check if also works for test
#
# rf_prediction <- predict(rft,X1)
#
#

#
#
#


drop <- which(apply(kidney1_gex_train_normed,1,function(x){all(x == 0)}))
kidney1_gex_train_normed_redux <- kidney1_gex_train_normed[-drop,]
tsvdin <- Matrix::Matrix(kidney1_gex_train_normed_redux,sparse  = TRUE)
trained_tsvd <- sparsesvd::sparsesvd(Matrix::t(tsvdin),rank = 300)
kidney1_gex_train_projected <- t(kidney1_gex_train_normed_redux) %*% trained_tsvd$v

# Best at 5? --> probably overfitting
mod1 <- stepFlexmix(kidney1_adt_train[1,] ~ kidney1_gex_train_projected,control = list(verbose = 0), k = 1:5,nrep = 3)

fm_models <- flexmix(~kidney1_gex_train_projected, k = 5, model = list(
  FLXMRglm(kidney1_adt_train[1,]~.),
  FLXMRglm(kidney1_adt_train[2,]~.),
  FLXMRglm(kidney1_adt_train[3,]~.),
  FLXMRglm(kidney1_adt_train[4,]~.),
  FLXMRglm(kidney1_adt_train[5,]~.),
  FLXMRglm(kidney1_adt_train[6,]~.),
  FLXMRglm(kidney1_adt_train[7,]~.)
))

fo <- sample(rep(seq(10), length = 2000))
mod1 <- stepFlexmix(~t(kidney1_gex_train_normed_redux[,1:2000]),model = list(
  FLXMRglmnet(kidney1_adt_train[1,1:2000]~.,foldid = fo),
  FLXMRglmnet(kidney1_adt_train[2,1:2000]~.,foldid = fo),
  FLXMRglmnet(kidney1_adt_train[3,1:2000]~.,foldid = fo),
  FLXMRglmnet(kidney1_adt_train[4,1:2000]~.,foldid = fo),
  FLXMRglmnet(kidney1_adt_train[5,1:2000]~.,foldid = fo),
  FLXMRglmnet(kidney1_adt_train[6,1:2000]~.,foldid = fo),
  FLXMRglmnet(kidney1_adt_train[7,1:2000]~.,foldid = fo)
),control = list(verbose = 0), k = 1:5,nrep = 3)
