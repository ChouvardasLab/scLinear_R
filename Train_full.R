#Trained on model 'as is' --> using the somewhat 'wrong' Seurat CLR and no further TSVD on ADT for ILR transform
devtools::load_all('R/')
set.seed(50)
for(file in dir('Three_tissue_input_data/')){
   data_source <- stringr::str_match(file,
                     '.*multiplex_([a-zA-Z]+_[a-zA-Z1-2]+)_.*')[,2]
   data <- Seurat::Read10X_h5(paste0('Three_tissue_input_data/',file))
   gex <- data$`Gene Expression`
   adt <- data$`Antibody Capture`
   combi_seurat <- Seurat::CreateSeuratObject(counts = gex, assay = 'RNA')
   adt_seurat <- Seurat:: CreateAssayObject(counts = adt)
   combi_seurat[["ADT"]] <- adt_seurat
   assign(data_source,prepare_data(combi_seurat))
   rm(adt,adt_seurat,combi_seurat,data,gex)
}

merged <- merge(Breast_Cancer,y=c(Kidney_Cancer1,Kidney_Cancer2,Lung_Cancer), add.cell.ids = c('Breast','Kidney1','Kidney2','Lung'),project = 'full')
rm(Lung_Cancer,Kidney_Cancer1,Kidney_Cancer2,Breast_Cancer)
merged <- SeuratObject::JoinLayers(merged)
indx <- sample(1:ncol(merged), size = ncol(merged), replace = FALSE)
train_cols <- colnames(merged)[indx[1:floor(0.75*ncol(merged))]]
test_cols <- colnames(merged)[indx[(floor(0.75*ncol(merged))+1):ncol(merged)]]

# In R it seems like rowaccess is faster for dgC and colaccess faster for dgE even though it should generally be the other way around
# For some reason, reordering is not respected as in
# gex_train <- Seurat::GetAssayData(all_train@assays[["RNA"]],layer = 'counts')
# gex_test <- Seurat::GetAssayData(all_test@assays[["RNA"]],layer = 'counts')
# adt_train <- Seurat::GetAssayData(all_train@assays[["ADT"]],layer = 'counts')
# adt_test <- Seurat::GetAssayData(all_test@assays[["ADT"]],layer = 'counts')
gex = Seurat::GetAssayData(merged@assays[["RNA"]],layer = 'counts')
adt = Seurat::GetAssayData(merged@assays[["ADT"]],layer = 'counts')
gex_train = Matrix::t(Matrix::t(gex)[train_cols,])
adt_train = Matrix::t(Matrix::t(adt)[train_cols,])
gex_test = Matrix::t(Matrix::t(gex)[test_cols,])
adt_test = Matrix::t(Matrix::t(adt)[test_cols,])


train_model <- fit_predictor(gex_train,adt_train,gexp_test = gex_test)
ev <- evaluate_predictor(train_model,gex_test,adt_test)
test_predictions <- predict(train_model,gex_test)
adt_test[adt_test == 0] <- 0.99
reference <- Matrix::Matrix(((compositions::clr(t(as.matrix(adt_test))))),sparse = TRUE)
errors <- test_predictions - reference
# # Unnormed counts are integers so all 0.99s are artificial --> no chance of accidentally modifying an original non-zero value
# adt_test[adt_test == 0.99] <- 0

# Train on ALL data (no train test split since we won't evaluate it again)
gex_all <- Seurat::GetAssayData(merged@assays[["RNA"]],layer = 'counts')
adt_all <- Seurat::GetAssayData(merged@assays[["ADT"]],layer = 'counts')
full_model <- fit_predictor(gex_all,adt_all)

library(ggplot2)
error_waterfall <- function(model_errors, xlabel = 'ADT', ylabel = 'Prediction error (CLR-transformed)',
                            title = 'Errors of different ADT count predictions for a single cell'){
df <- data.frame(x = names(model_errors), y = model_errors)
display_numbers <- round(model_errors,2)
colnames(df) <- c('ADT','error')
waterfalls::waterfall(df, rect_text_labels = display_numbers) + theme_minimal() +
  xlab(xlabel) + ylab(ylabel) +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

error_waterfall(errors[1,])
error_waterfall(test_predictions[1,],ylabel = 'Predicted count (CLR-transformed)',title = 'ADT count predictions for a single cell')

gmeans <- apply(adt_test,2,function(x){exp(mean(log(x)))})




true_vals <- reference[2,]
preds <- test_predictions[2,]
fixed_effects <- abind::abind(train_model$lm_coefficients,along = 2)[1,]
df <- data.frame(x = names(true_vals), true_val = true_vals, prediction = preds, fixed = fixed_effects)
df$dynamic <- df$prediction - df$fixed
combi <- data.frame(matrix(nrow = 3 * nrow(df), ncol = 2))
combi[seq(1,nrow(combi),3),1] <- paste0(df$x,'_true')
combi[seq(2,nrow(combi),3),1] <- paste0(df$x,'_fixed')
combi[seq(3,nrow(combi),3),1] <- paste0(df$x,'_offset')
combi[seq(1,nrow(combi),3),2] <- df$true_val
combi[seq(2,nrow(combi),3),2] <- df$fixed
combi[seq(3,nrow(combi),3),2] <- df$dynamic
colnames(combi) <- c('label','value')
display_numbers <- round(combi$value,2)
waterfalls::waterfall(combi, rect_text_labels = display_numbers, fill_by_sign = FALSE,fill_colours = rep(c(1,2,3),nrow(combi)/3)) + theme_minimal() +
  xlab('Component') + ylab('Value') +
  ggtitle('Error decomposition') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

comparison <- t(apply(reference,1,function(x){x-test_predictions[1,]}))


#positive = #00BFC4
#negative = #F8766D
error_bar <- function(model_errors, xlabel = 'ADT', ylabel = 'Prediction error (CLR-transformed)',
                            title = 'Errors of different ADT count predictions for a single cell'){
  df <- data.frame(adt = as.factor(names(model_errors)), val = model_errors, sign = as.factor(sign(model_errors)))
  # Reordering to preserve original order --> easier comparison to waterfall which preserves order by default
  df$adt <- factor(df$adt, levels = df$adt)
  display_numbers <- round(model_errors,2)
  ggplot(df,aes(x = adt,y = val, fill = sign))+ geom_bar(stat = 'identity',colour = 'black', width = 0.7) +
    theme_minimal() + scale_fill_manual(values = c('1'='#00BFC4','0' = 'grey','-1'='#F8766D')) +
    xlab(xlabel) + ylab(ylabel) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'None') +
    geom_hline(yintercept = 0)
}

error_waterfall(errors[1,])
error_bar(errors[1,])


source('Plot_comparison.R')

adt_df <- data.frame(reference)
adt_df <- data.frame(cbind(as.factor(merged$cell_type[rownames(reference)]),adt_df))
# adt_df <- data.frame(cbind(prepped$seurat_clusters,t(adt_df)))
adt_df[,1] <- as.factor(adt_df[,1])
colnames(adt_df)[1] <- 'cell_type'
ggplot(adt_df,aes(x = `CD3.CD3E`)) + geom_histogram(fill ='grey',color = 'black') +
  xlab('CD3 count (CLR-transformed)') +
  ylab('Number of cells') + ggtitle('ADT level distribution') + ylim(c(0,1500))
ggsave('ADT_dist_total.jpg')
ggplot(adt_df,aes(x = `CD3.CD3E`, fill = cell_type, alpha = 0.5)) + geom_histogram(position = 'identity',color = 'black') +
  scale_color_discrete(name = "", labels = levels(as.factor(merged$cell_type[rownames(reference)]))) + xlab('CD3 count (CLR-transformed)') +
  ylab('Number of cells') + ggtitle('ADT level distribution') + ylim(c(0,1500))
ggsave('ADT_dist_by_celltype.jpg')

plot_comparison(reference[,1],test_predictions[,1], xaxis = 'True ADT value (CLR transformed)', yaxis = 'Prediction', title = 'CD3 predictions')

# From here train split by celltype
test_predictions_split <- matrix(nrow = nrow(reference),ncol = ncol(reference))
rownames(test_predictions_split) <- rownames(reference)
colnames(test_predictions_split) <- colnames(reference)
gex_projected <-  t(filter_input_genes(gexp_normalize(gex),train_model)) %*% train_model$tsvd_v
for(ct in unique(merged$cell_type)){
# Which rows of the combined predictions matrix should be filled in with the results from this loop
cs <- names(merged$cell_type)[merged$cell_type == ct]
gex_train_ct <- Matrix::t(gex_projected)[,intersect(cs,colnames(gex_train))]
training_set <- Matrix::t(apply(gex_train_ct, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
adt_train_ct <- Matrix::t(Matrix::t(adt_train)[intersect(cs,colnames(adt_train)),])
adt_train_ct[adt_train_ct == 0] <- 0.99
adt_train_ct <- Matrix::Matrix((t(compositions::clr(t(as.matrix(adt_train_ct))))),sparse = TRUE)
gex_test_ct <- Matrix::t(gex_projected)[,intersect(cs,colnames(gex_test))]
test_set <- Matrix::t(apply(gex_test_ct, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
for(i in 1:nrow(adt_train_ct)){
lmod <- lm(adt_train_ct[i,] ~ t(training_set))
ct_predi <- lmod$coefficients[1] + t(test_set) %*% lmod$coefficients[2:301]
test_predictions_split[intersect(cs,colnames(gex_test)),i] <- ct_predi
}
}

for(i in colnames(test_predictions_split)){
  plot_comparison(reference[,i],test_predictions_split[,i],xaxis = paste0('True ',i,' (CLR transformed)'),yaxis = paste0('Predicted ',i,' (CLR transformed)'),
                  title = 'Prediction vs. true ADT values', groupings = merged$cell_type[rownames(reference)])
    ggsave(paste0('./Prediction_evaluation_plots/',i,'_split.jpg'))
  plot_comparison(reference[,i],test_predictions[,i],xaxis = paste0('True ',i,' (CLR transformed)'),yaxis = paste0('Predicted ',i,' (CLR transformed)'),
                    title = 'Prediction vs. true ADT values')
  ggsave(paste0('./Prediction_evaluation_plots/',i,'.jpg'))
}

i <- 'CD10-MME'
p <- plot_comparison(reference[,i],test_predictions[,i],xaxis = paste0('True ',i,' (CLR transformed)'),yaxis = paste0('Predicted ',i,' (CLR transformed)'),
                title = 'Prediction vs. true ADT values')
p <- p+ xlim(2.5,6) + ylim(2.5,5.5)
ggsave(paste0('./Prediction_evaluation_plots/',i,'_zoomed.jpg'))


#Somehow show the lm offsets and 'effect' of tsvd components in the plot as well - maybe this will somehow give a clue regarding the source of error linkage
plot_comparison(as.vector(reference[1:1000,]),as.vector(test_predictions[1:1000,]),groupings = rep(colnames(reference),1000))


train_model$lm_coefficients

