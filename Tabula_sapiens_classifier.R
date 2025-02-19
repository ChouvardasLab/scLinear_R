# tsap_dat <- readRDS('Tabula_sapiens_data/tabula_sapiens_blood.rds')
# #Cell_types directly under tsap_dat$cell_type. They also have other annotation info stored but this seems to be
# #What they report as THE cell type in their own data exploration tool
# rna_counts <- Seurat::GetAssayData(tsap_dat@assays[["RNA"]],layer = 'counts')
# #Columns = cells, Rows = features
# feature_names <- tsap_dat@assays$RNA@meta.features$feature_name
#
#
# reticulate::use_condaenv("C:/Users/mz24b548/AppData/Local/miniconda3/envs/scLinearDependencies")
predictor <- readRDS('three_tissue_model_full')
# print(paste0(
#   length(intersect(predictor$genes_considered,feature_names))
#   ,'/', length(predictor$genes_considered),' genes used in training also measured/matched in this dataset'
#   )
# )
# rownames(rna_counts) <- feature_names
#
#
devtools::load_all('R/')
# subsample <- colnames(rna_counts)[1:2000]
#
# #Selecting rows for dgC much faster than columns it seems --> transpose --> select rows --> transpose
# adt_prediction <- predict(predictor,Matrix::t(Matrix::t(rna_counts)[subsample,]))
#
# adt_df <- data.frame(cbind(adt_prediction[,1],as.factor(tsap_dat$cell_type)))
# adt_df[,2] <- as.factor(adt_df[,2])
# colnames(adt_df) <- c('prediction','cell_type')
# ggplot(adt_df,aes(x = `prediction`, color = cell_type)) + geom_histogram(position = 'identity', fill = NA) +
#   scale_color_discrete(name = "", labels = levels(as.factor(tsap_dat$cell_type)))
tsap_dat <- readRDS('Tabula_sapiens_data/tabula_sapiens_blood.rds')
rna_counts <- Seurat::GetAssayData(tsap_dat@assays[["RNA"]],layer = 'counts')
feature_names <- tsap_dat@assays$RNA@meta.features$feature_name
rownames(rna_counts) <- feature_names
predictions <- readRDS('tabula_sapiens_blood_adt_predictions')
indx <- sample(1:nrow(predictions), size = nrow(predictions), replace = FALSE)
train_indices <- indx[1:floor(0.75*length(indx))]
test_indices <- indx[(floor(0.75*length(indx))+1):length(indx)]
svm <- e1071::svm(predictions[train_indices,],tsap_dat$cell_type[train_indices])
test_predictions <- predict(svm,predictions[test_indices,])
correct <- test_predictions == tsap_dat$cell_type[test_indices]
wrong <- test_predictions != tsap_dat$cell_type[test_indices]
print(
  paste0(100*round(sum(correct)/(sum(correct)+sum(wrong)),3),'% classified correctly'
  )
)


classifications <- as.data.frame(table(test_predictions,tsap_dat$cell_type[test_indices]))
colnames(classifications)[1:2] <- c('Predicted class','True annotation')
order <- c(
  "hematopoietic precursor cell"
  ,"hematopoietic stem cell"
  ,"common myeloid progenitor"
  ,"CD4-positive, alpha-beta T cell"
  ,"naive thymus-derived CD4-positive, alpha-beta T cell"
  ,"CD8-positive, alpha-beta T cell"
  ,"mature NK T cell"
  ,"regulatory T cell"
  ,"monocyte"
  ,"classical monocyte"
  ,"intermediate monocyte"
  ,"non-classical monocyte"
  ,"erythrocyte"
  ,"platelet"
  ,"macrophage"
  ,"B cell"
  ,"natural killer cell"
  ,"basophil"
  ,"neutrophil"
  ,"myeloid dendritic cell"
  ,"plasmacytoid dendritic cell"
  ,"plasma cell"
)
classifications$`Predicted class` <- factor(classifications$`Predicted class`,
levels = order)
classifications$`True annotation` <- factor(classifications$`True annotation`,
                                            levels = order)


rna_counts_normed <-

gex_filtered <- filter_input_genes(rna_counts,predictor)
gex_projected <- Matrix::t(gex_filtered) %*% predictor$tsvd_v
gex_train <- Matrix::t(Matrix::t(gex_projected)[,train_indices])
gex_projected_svd_train <- reticulate::r_to_py(gex_train)

gex_tabsap <- reticulate::r_to_py(rna_counts)
#Easier than passing all row-/ & columnnames to py and doing filtering there
gex_filtered <- reticulate::r_to_py(gex_filtered)
adt_predictions_for_tabsap <- reticulate::r_to_py(predictions)
tsvd_projector <- reticulate::r_to_py(predictor$tsvd_v)
tabsap_annotations <- reticulate::r_to_py(tsap_dat$cell_type)

reticulate::py_save_object(gex_tabsap,'Tabula_sapiens_classifier_input/gene_expression_tabsap')
reticulate::py_save_object(adt_predictions_for_tabsap,'Tabula_sapiens_classifier_input/adt_predictions_for_tabsap')
reticulate::py_save_object(tsvd_projector,'Tabula_sapiens_classifier_input/tsvd_projector_from_3tmodel')
reticulate::py_save_object(tabsap_annotations,'Tabula_sapiens_classifier_input/tabsap_annotations')
reticulate::py_save_object(gex_filtered,'Tabula_sapiens_classifier_input/gex_filtered_for_tsvd')

svm_projected_gex <- e1071::svm(gex_train,tsap_dat$cell_type[train_indices])

# library(dplyr)
# total_occurences_in_prediction <- classifications %>% group_by(`Predicted class`) %>% summarise(total_occurences = sum(Freq))
# total_occurences_in_annotation <- classifications %>% group_by(`True annotation`) %>% summarise(total_occurences = sum(Freq))
# drop_prediction_levels <- total_occurences_in_prediction[total_occurences_in_prediction$total_occurences == 0,'Predicted class']
# drop_annotation_levels <- total_occurences_in_annotation[total_occurences_in_annotation$total_occurences == 0,'True annotation']
# classifications <- classifications[classifications$`Predicted class` n ]
# redux <- classifications[!(classifications$`Predicted class` %in% as.list(drop_prediction_levels)[[1]])
#                 &!(classifications$`True annotation` %in% as.list(drop_annotation_levels)[[1]])
#                 ,]



ggplot(classifications,aes(x = `True annotation`,y = `Predicted class`)) +
  theme_bw() +
  geom_tile(colour = "black",fill = 'white') +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5)) +
  geom_text(aes(label = Freq), color = "black", fontface = "bold")

# Good chance this will crash / never finish
# gene_svm <- e1071::svm(Matrix::t(rna_counts)[train_indices,],tsap_dat$cell_type[train_indices])





# misclassifications <- as.data.frame(table(test_predictions[wrong],tsap_dat$cell_type[test_indices][wrong]))
# colnames(misclassifications)[1:2] <- c('Prediction','True Annotation')
# misclassifications <- misclassifications[misclassifications[,3] != 0,]

# Right of the bat --> 74901 correct, 10332 wrong --> split into train/test first for proper testing



