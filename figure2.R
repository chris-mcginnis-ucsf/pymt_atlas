#######################################################################
## The temporal progression of immune remodeling during metastasis  ###
## Companion code: Figure 2, Chris McGinnis, PhD, Stanford, 05/2023  ##
#######################################################################
## objects
load('seu_int.Robj')
load('seu_alv.Robj')
load('anno_freq_summary_int.Robj')
load('anno_freq_summary_alv.Robj')
load('prop_test_int.Robj')
load('prop_test_alv.Robj')
load('crip1_zscores.Robj')
load('umap_alv.Robj')

## functions
load("stderror.Robj")
load("density_plot.Robj")

############################################################################################################################
## Figure 2. AM inflammation and IM wound healing transcriptional signatures are linked to pre-metastatic niche formation ##
############################################################################################################################
## Fig. 2A: IM subtype annotations, feature plots, and subtype proportion bar charts
# GEX space
DimPlot(seu_int, group.by = 'subtype', cols = c('darkred','goldenrod','red','grey'), pt.size = 1) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Mrc1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Cd74', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Crip1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Vim', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

# Subtype proportions
# anno_freq <- table(seu_int@meta.data$age_rep, seu_int@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_int <- as.data.frame(matrix(0L, nrow=4*4, ncol=4))
# colnames(anno_freq_summary_int) <- c('celltype','clade','mean','se')
# anno_freq_summary_int$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_int$clade <- rep(unique(anno_freq$clade), 4)
# for (i in 1:nrow(anno_freq_summary_int)) {
#   ct <- anno_freq_summary_int$celltype[i]
#   cl <- anno_freq_summary_int$clade[i]
#   anno_freq_summary_int$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_int$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# anno_freq_summary_int$clade <- factor(anno_freq_summary_int$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_int$celltype <- factor(anno_freq_summary_int$celltype, levels=c('antigen','mrc1','crip1','prolif'))
ggplot(anno_freq_summary_int, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics
# sce_int <- as.SingleCellExperiment(seu_int)
# i.int <- seq_len(ncol(sce_int))
# boot.int1 <- sample(i.int, replace=TRUE)
# boot.int2 <- sample(i.int, replace=TRUE)
# sce_int_rep2 <- sce_int[,boot.int1]
# sce_int_rep3 <- sce_int[,boot.int2]
# 
# sample <- c(sce_int$age_rep,
#             paste0(sce_int_rep2$age_rep,'_boot2'),
#             paste0(sce_int_rep3$age_rep,'_boot3'))
# cluster <- c(sce_int$subtype,
#              sce_int_rep2$subtype,
#              sce_int_rep3$subtype)
# group <- c(sce_int$clade,
#            sce_int_rep2$clade,
#            sce_int_rep3$clade)
# allcounts <- cbind(counts(sce_int),
#                    counts(sce_int_rep2),
#                    counts(sce_int_rep3))
# 
# sce_int_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_int_boot$sample <- sample
# sce_int_boot$group <- group
# sce_int_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_int_boot$group)
# ind_early <- grep('early',sce_int_boot$group)
# ind_mid <- grep('mid',sce_int_boot$group)
# ind_late <- grep('late',sce_int_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_int_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_int_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_int_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_int_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_int_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_int_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_int_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_int_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_int_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_int_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_int_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_int_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_int_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_int_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_int_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_int_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_int_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_int_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_int <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_int) <- c('comparison','celltype','pval')
# prop_test_int$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_int$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_int$pval <- pval_vec
# prop_test_int$pval_bin <- 'ns'
# prop_test_int$pval_bin[which(prop_test_int$pval <= 0.05)] <- '*'

## Fig. 2B: Crip1-high Cav1+ signature z-score heatmap 
# crip1_genes <- c('Tppp3','Crip1','Ecm1','Vim','S100a4','Cspg4','Cav1','Aqp1','Spp1','S100a10','Lgals1','Ckb','Pmepa1','Tagln2','Lmna','Fn1','Nfic','Plec','Aopep','Olfml3', 'Tubb5','Pdlim1','Rgcc','Myof','Cyb5r3','Klf4','Ahnak','Anxa1','Lrp1','Capn2','Pltp','Mgst3','Lair1','Pf4','Emp1')
# crip1_genes <- crip1_genes[which(crip1_genes %in% seu_int@assays$RNA@var.features)]
# crip1_zscores <- as.data.frame(matrix(0L, nrow=4, ncol=length(crip1_genes)))
# rownames(crip1_zscores) <- c('antigen','mrc1','crip1','prolif')
# colnames(crip1_zscores) <- crip1_genes
# crip1_genes_mean <- rowMeans(seu_int@assays$RNA@data[crip1_genes, ])
# crip1_genes_sd <- apply(seu_int@assays$RNA@data[crip1_genes, ], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(crip1_zscores)) {
#   ind <- which(seu_int@meta.data$subtype == i)
#   temp <-  rowMeans(seu_int@assays$RNA@data[crip1_genes,ind])
#   crip1_zscores[i, ] <- (temp-crip1_genes_mean)/crip1_genes_sd
# }
Heatmap(crip1_zscores[1:3,], cluster_rows = F, show_heatmap_legend = F)

## Fig. 2C: AM subtype annotations and Cd14+ inflammatory signature feature plots
DimPlot(seu_alv, group.by = 'subtype', cols = c('maroon','goldenrod','black','tan4','red','lightsalmon','pink2','grey'), pt.size = 0.1)+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Tnf', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Cd14', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Cxcl2', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Nlrp3', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. 2D: AM metastatic stage density UMAPs and subtype proportion barcharts 
# Density UMAPs
# umap_alv <- as.data.frame(seu_alv@reductions$umap@cell.embeddings)
# umap_alv[,"clade"] <- seu_alv@meta.data$clade
density_plot(data = umap_alv, clade = 'wt', highlightcol = 'black', alpha = 0.5, bins=5)
density_plot(data = umap_alv, clade = 'early', highlightcol = 'black', alpha = 0.5, bins=5)
density_plot(data = umap_alv, clade = 'mid', highlightcol = 'black', alpha = 0.5, bins=5)
density_plot(data = umap_alv, clade = 'late', highlightcol = 'black', alpha = 0.5, bins=5)

# Subtype proportions
# anno_freq <- table(seu_alv@meta.data$age_rep, seu_alv@meta.data$subtype_new)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_alv <- as.data.frame(matrix(0L, nrow=4*8, ncol=4))
# colnames(anno_freq_summary_alv) <- c('celltype','clade','mean','se')
# anno_freq_summary_alv$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_alv$clade <- rep(unique(anno_freq$clade), 8)
# for (i in 1:nrow(anno_freq_summary_alv)) {
#   ct <- anno_freq_summary_alv$celltype[i]
#   cl <- anno_freq_summary_alv$clade[i]
#   anno_freq_summary_alv$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_alv$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_alv$clade <- factor(anno_freq_summary_alv$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_alv$celltype <- factor(anno_freq_summary_alv$celltype, levels=c('homeo','inf','lipid','allergy','antigen','mt','IFN','prolif'))
ggplot(anno_freq_summary_alv, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_alv <- as.SingleCellExperiment(seu_alv)
# i.alv <- seq_len(ncol(sce_alv))
# boot.alv1 <- sample(i.alv, replace=TRUE)
# boot.alv2 <- sample(i.alv, replace=TRUE)
# sce_alv_rep2 <- sce_alv[,boot.alv1]
# sce_alv_rep3 <- sce_alv[,boot.alv2]
# 
# sample <- c(sce_alv$age_rep,
#             paste0(sce_alv_rep2$age_rep,'_boot2'),
#             paste0(sce_alv_rep3$age_rep,'_boot3'))
# cluster <- c(sce_alv$subtype_new,
#              sce_alv_rep2$subtype_new,
#              sce_alv_rep3$subtype_new)
# group <- c(sce_alv$clade,
#            sce_alv_rep2$clade,
#            sce_alv_rep3$clade)
# allcounts <- cbind(counts(sce_alv),
#                    counts(sce_alv_rep2),
#                    counts(sce_alv_rep3))
# 
# sce_alv_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_alv_boot$sample <- sample
# sce_alv_boot$group <- group
# sce_alv_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_alv_boot$group)
# ind_early <- grep('early',sce_alv_boot$group)
# ind_mid <- grep('mid',sce_alv_boot$group)
# ind_late <- grep('late',sce_alv_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_alv_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_alv_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_alv_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_alv_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_alv_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_alv_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_alv_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_alv_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_alv_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_alv_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_alv_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_alv_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_alv_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_alv_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_alv_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_alv_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_alv_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_alv_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_alv <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_alv) <- c('comparison','celltype','pval')
# prop_test_alv$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_alv$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_alv$pval <- pval_vec
# prop_test_alv$pval_bin <- 'ns'
# prop_test_alv$pval_bin[which(prop_test_alv$pval <= 0.05)] <- '*'