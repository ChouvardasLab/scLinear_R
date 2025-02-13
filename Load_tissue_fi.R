fi_kidney <- readRDS('feature_importance_kidney1')
fi_breast <- readRDS('feature_importance_breast')
fi_lung <- readRDS('feature_importance_lung')

topX <- function(feature_importance,adt_name,X = 50){
  return(feature_importance[names(head(sort(abs(feature_importance[,adt_name]),decreasing = TRUE),X)),adt_name])
}

DrawVenn <- function(fi1,fi2,fi3,n_top,adt_name,names){
  top_fi1 <- topX(fi1,adt_name,n_top)
  top_fi2 <- topX(fi2,adt_name,n_top)
  top_fi3 <- topX(fi3,adt_name,n_top)
  VennDiagram::venn.diagram(list(names(top_fi1),names(top_fi2),names(top_fi3)),category.names = names,paste0('Top_',n_top,'genes_overlap_',adt_name,'.png'),disable.logging = TRUE)
}

show <- 100
for(ADT in colnames(fi_kidney)){
check <- topX(fi_kidney,ADT,show)
check <- sort(check,decreasing = TRUE)
refgenes <- names(check)
refgenes_importance_breast <- fi_breast[refgenes,ADT]
refgenes_importance_lung <- fi_lung[refgenes,ADT]

#norm by most important
check <- check/max(abs(check))
refgenes_importance_breast <- refgenes_importance_breast/max(abs(refgenes_importance_breast))
refgenes_importance_lung <- refgenes_importance_lung/max(abs(refgenes_importance_lung))

plotdat <- matrix(nrow = 3*show, ncol = 3)
plotdat[,1] <- c(rep('Kidney',show),rep('Breast',show),rep('Lung',show))
plotdat[,2] <- rep(refgenes,3)
plotdat[,3] <- c(check,refgenes_importance_breast,refgenes_importance_lung)
colnames(plotdat) <- c('Tissue','Gene','Feature_importance')
ggplot(plotdat,aes(x = forcats::fct_inorder(Tissue), y = forcats::fct_inorder(Gene), fill = as.numeric(Feature_importance))) + geom_tile() + scale_fill_gradient() +
xlab('Cancer Type') + ylab('Gene') + ggtitle(ADT)
ggsave(paste0('FI_comparison_plots/',ADT,'_feature_importance.png'))
}

topX_up <- function(feature_importance,adt_name,X){
  return(head(sort(feature_importance[,adt_name],decreasing =  TRUE),X))
}

topX_down <- function(feature_importance,adt_name,X){
  return(head(sort(feature_importance[,adt_name]),X))
}

count <- 50
for(prot in colnames(fi_kidney)){
  kidney_upregulators <- names(topX_up(fi_kidney,prot,count))
  kidney_downregulators <- names(topX_down(fi_kidney,prot,count))
  effect_lung_upregulators <- fi_lung[kidney_upregulators,prot]
  effect_lung_downregulators <- fi_lung[kidney_downregulators,prot]
  effect_breast_upregulators <- fi_breast[kidney_upregulators,prot]
  effect_breast_downregulators <- fi_breast[kidney_downregulators,prot]
  lung_effects <- data.frame(rbind(
    cbind(effect_lung_upregulators,effect_in_kidney = rep(1,length(kidney_upregulators)))
    ,cbind(effect_lung_downregulators,effect_in_kidney = rep(-1,length(kidney_downregulators)))
    )
  )
  breast_effects <- data.frame(rbind(
    cbind(effect_breast_upregulators,effect_in_kidney = rep(1,length(kidney_upregulators)))
    ,cbind(effect_breast_downregulators,effect_in_kidney = rep(-1,length(kidney_downregulators)))
    )
  )

  colnames(lung_effects)[1] <- 'effect_lung'
  colnames(breast_effects)[1] <- 'effect_breast'

  lung_effects[,2] <- as.factor(lung_effects[,2])
  breast_effects[,2] <- as.factor(breast_effects[,2])

  plung <- ggplot(lung_effects,aes(x = effect_in_kidney, y = effect_lung)) +
    geom_boxplot(outlier.colour = 'black',outlier.shape=16, outlier.size=2, notch=FALSE) +
    xlab('Effect in kidney') + ylab('Feature importance') + ggtitle('Effect in lung') +
    scale_x_discrete(labels=c("Negative effect","Positive effect")) + ylim(-1,1)

  pbreast <- ggplot(breast_effects,aes(x = effect_in_kidney, y = effect_breast)) +
    geom_boxplot(outlier.colour = 'black',outlier.shape=16, outlier.size=2, notch=FALSE) +
    xlab('Effect in kidney') + ylab('Feature importance') + ggtitle('Effect in breast') +
    scale_x_discrete(labels=c("Negative effect","Positive effect")) + ylim(-1,1)



  ggsave(paste0('FI_comparison_plots/',prot,'_effect_tissue_comparison.png'),
         gridExtra::arrangeGrob(plung,pbreast,top = paste0('Effect of top ',count,' positive and negative ',prot,' predictors in other tissues'),ncol = 2)
         )
  }




