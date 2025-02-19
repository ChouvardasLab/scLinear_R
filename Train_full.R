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
merged <- SeuratObject::JoinLayers(merged)
indx <- sample(1:ncol(merged), size = ncol(merged), replace = FALSE)
all_train <- merged[,indx[1:floor(0.75*ncol(merged))]]
all_test <- merged[,indx[(floor(0.75*ncol(merged))+1):ncol(merged)]]

gex_train <- Seurat::GetAssayData(all_train@assays[["RNA"]],layer = 'counts')
gex_test <- Seurat::GetAssayData(all_test@assays[["RNA"]],layer = 'counts')
adt_train <- Seurat::GetAssayData(all_train@assays[["ADT"]],layer = 'counts')
adt_test <- Seurat::GetAssayData(all_test@assays[["ADT"]],layer = 'counts')

train_model <- fit_predictor(gex_train,adt_train,gexp_test = gex_test)
ev <- evaluate_predictor(train_model,gex_test,adt_test)
test_predictions <- predict(train_model,gex_test)
# adt_test[adt_test == 0] <- 0.99
# reference <- Matrix::Matrix(((compositions::clr(t(as.matrix(adt_test))))),sparse = TRUE)
# # Unnormed counts are integers so all 0.99s are artificial --> no chance of accidentally modifying an original non-zero value
# adt_test[adt_test == 0.99] <- 0

rm(Breast_Cancer,Kidney_Cancer1,Kidney_Cancer2,Lung_Cancer)

# Train on ALL data (no train test split since we won't evaluate it again)
gex_all <- Seurat::GetAssayData(merged@assays[["RNA"]],layer = 'counts')
adt_all <- Seurat::GetAssayData(merged@assays[["ADT"]],layer = 'counts')
full_model <- fit_predictor(gex_all,adt_all)

