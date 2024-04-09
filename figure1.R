##########################################################################################
## The temporal progression of lung immune remodeling during breast cancer metastasis  ###
## Companion code: Figure 1, Chris McGinnis, PhD, Stanford, October 2023  ################
##########################################################################################
## objects
load('seu_imm_fin_2.Robj')
load('immune_markers.Robj')
load('jsd_imm.Robj')
load('anno_freq_summary_imm.Robj')
load('prop_test_imm.Robj')

## functions
load("jsd_comp.Robj")
load('stderror.Robj')

######################################################################################################################################
## Figure 1: Longitudinal scRNA-seq cell atlas of PyMT mouse lung immune cells captures dynamics of the metastatic microenvironment ##
######################################################################################################################################
## Fig. 1A: Experimental schematic
## Fig. 1B: Immune cell type annotation UMAP, marker gene dot plots 
DimPlot(seu_imm_fin_2, group.by='celltype_fig1', cols=c('skyblue3','indianred1','black','darkred','red','turquoise4','turquoise','dodgerblue','lightsalmon','darkorchid4','royalblue3','navy','palegreen3','maroon','turquoise3','grey','pink3'), raster=F, pt.size = 0.01) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

# markers_imm <- FindAllMarkers(seu_imm_fin_2, only.pos = T, logfc.threshold = log(2.5))
# immune_markers <- NULL
# for (i in unique(markers_imm$cluster)) { immune_markers <- c(immune_markers, markers_imm$gene[which(markers_imm$cluster == i)[1:3]]) }
g <- DotPlot(seu_imm_fin_2, features=unique(immune_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(levels(g$data$id)))

## Fig. 1C: Timepoint JSD heatmap 
# jsd_imm <- jsd_comp(seu = seu_imm_fin_2, condition = 'age_rep')
Heatmap(jsd_imm, col = brewer.pal(10,'RdGy'), show_column_names = F, row_names_gp = gpar(fontsize=8), show_heatmap_legend = F)

## Fig. 1D: Cell type proportion bar charts
# clade_wt <- c('wt_rep2','wt_rep3')
# clade_early <- c(paste0('6wk_rep',1:3), paste0('7wk_rep',c(1,3)), paste0('8wk_rep',1:3), paste0('9wk_rep',2:3))
# clade_mid <- c('9wk_rep1', paste0('10wk_rep',c(1,2,4,6)), '11wk_rep1', '14wk_rep3')
# clade_late <- c('10wk_rep5', paste0('11wk_rep',2:3), paste0('12wk_rep',1:3), paste0('13wk_rep',2:3), paste0('14wk_rep',1:2))
# 
# anno_freq <- table(seu_imm_fin_2@meta.data$age_rep, seu_imm_fin_2@meta.data$celltype_fig1)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_imm <- as.data.frame(matrix(0L, nrow=17*4, ncol=4))
# colnames(anno_freq_summary_imm) <- c('celltype','clade','mean','se')
# anno_freq_summary_imm$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_imm$clade <- rep(unique(anno_freq$clade), 17)
# for (i in 1:nrow(anno_freq_summary_imm)) {
#   ct <- anno_freq_summary_imm$celltype[i]
#   cl <- anno_freq_summary_imm$clade[i]
#   anno_freq_summary_imm$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_imm$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# anno_freq_summary_imm$clade <- factor(anno_freq_summary_imm$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_imm$celltype <- factor(anno_freq_summary_imm$celltype, levels=levels(seu_imm_fin_2@meta.data$celltype_fig1))
ggplot(anno_freq_summary_imm, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics
# sce_imm <- as.SingleCellExperiment(seu_imm_fin_2)
# i.imm <- seq_len(ncol(sce_imm))
# boot.imm1 <- sample(i.imm, replace=TRUE)
# boot.imm2 <- sample(i.imm, replace=TRUE)
# sce_imm_rep2 <- sce_imm[,boot.imm1]
# sce_imm_rep3 <- sce_imm[,boot.imm2]
# 
# sample <- c(sce_imm$age_rep,
#             paste0(sce_imm_rep2$age_rep,'_boot2'),
#             paste0(sce_imm_rep3$age_rep,'_boot3'))
# cluster <- c(sce_imm$celltype_fig1,
#              sce_imm_rep2$celltype_fig1,
#              sce_imm_rep3$celltype_fig1)
# group <- c(sce_imm$clade,
#            sce_imm_rep2$clade,
#            sce_imm_rep3$clade)
# allcounts <- cbind(counts(sce_imm),
#                    counts(sce_imm_rep2),
#                    counts(sce_imm_rep3))
# 
# sce_imm_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_imm_boot$sample <- sample
# sce_imm_boot$group <- group
# sce_imm_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_imm_boot$group)
# ind_early <- grep('early',sce_imm_boot$group)
# ind_mid <- grep('mid',sce_imm_boot$group)
# ind_late <- grep('late',sce_imm_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_imm_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_imm_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_imm_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_imm_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_imm_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_imm_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_imm_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_imm_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_imm_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_imm_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_imm_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_imm_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_imm_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_imm_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_imm_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_imm_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_imm_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_imm_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_imm <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_imm) <- c('comparison','celltype','pval')
# prop_test_imm$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_imm$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_imm$pval <- pval_vec
# prop_test_imm$pval_bin <- 'ns'
# prop_test_imm$pval_bin[which(prop_test_imm$pval <= 0.05 & prop_test_imm$pval > 1e-5)] <- '*'
# prop_test_imm$pval_bin[which(prop_test_imm$pval <= 1e-5 & prop_test_imm$pval > 1e-10)] <- '**'
# prop_test_imm$pval_bin[which(prop_test_imm$pval <= 1e-10)] <- '***'
# prop_test_imm$comparison <- factor(prop_test_imm$comparison, levels = c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'))
# prop_test_imm$celltype <- factor(prop_test_imm$celltype, levels=levels(seu_imm_fin_2@meta.data$celltype_fig1))
# prop_test_imm$pval_bin <- factor(prop_test_imm$pval_bin, levels=c('ns','*','**','***'))
prop_test_imm <- acast(prop_test_imm, comparison~celltype, value.var="pval_bin")
Heatmap(prop_test_imm, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[2:4],0.8)))


