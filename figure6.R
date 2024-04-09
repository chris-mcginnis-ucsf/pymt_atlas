#######################################################################
## The temporal progression of immune remodeling during metastasis  ###
## Companion code: Figure 6, Chris McGinnis, PhD, Stanford, 03/2024  ##
#######################################################################
## objects
load('seu_nk.Robj')
load('umap_nk.Robj')
load('anno_freq_summary_nk.Robj')

load('seu_t.Robj')
load('seu_b.Robj')
load('seu_nk_kim2020.Robj')
load('anno_freq_summary_nk_kim.Robj')
load('anno_freq_summary_t.Robj')
load('anno_freq_summary_b.Robj')
load('prop_test_t.Robj')
load('prop_test_b.Robj')
load('nk_facs_data_summary.Robj')
load('nk_facs_counts.Robj')
load('nk_facs_pval_summary.Robj')
load('t_cell_zscores.Robj')

## functions
load("stderror.Robj")
load("density_plot.Robj")

#########################################################################################################################################
## Figure 6: Lymphocyte subtype characterization reveals details of the inflammatory and immunosuppressive metastatic microenvironment ##
#########################################################################################################################################
## Fig. 6A: NK subtype annotations, metastatic stage density UMAPs, and subtype proportion barcharts
# GEX space
DimPlot(seu_nk, group.by = 'subtype', cols=c('forestgreen','goldenrod','palegreen3','black'), pt.size = 1) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

# Density UMAPs
# umap_nk <- as.data.frame(seu_nk@reductions$umap@cell.embeddings)
# umap_nk[,"clade"] <- seu_nk@meta.data$clade
density_plot(data = umap_nk, clade = 'wt', highlightcol = 'black', alpha = 0.5, bins=5)
density_plot(data = umap_nk, clade = 'early', highlightcol = 'black', alpha = 0.5, bins=5)
density_plot(data = umap_nk, clade = 'mid', highlightcol = 'black', alpha = 0.5, bins=5)
density_plot(data = umap_nk, clade = 'late', highlightcol = 'black', alpha = 0.5, bins=5)

# Subtype proportions
# anno_freq <- table(seu_nk@meta.data$age_rep, seu_nk@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_nk <- as.data.frame(matrix(0L, nrow=4*4, ncol=4))
# colnames(anno_freq_summary_nk) <- c('celltype','clade','mean','se')
# anno_freq_summary_nk$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_nk$clade <- rep(unique(anno_freq$clade), 4)
# for (i in 1:nrow(anno_freq_summary_nk)) {
#   ct <- anno_freq_summary_nk$celltype[i]
#   cl <- anno_freq_summary_nk$clade[i]
#   anno_freq_summary_nk$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_nk$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_nk$clade <- factor(anno_freq_summary_nk$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_nk$celltype <- factor(anno_freq_summary_nk$celltype, levels=c('mod','cyto','cyto_ccl','prolif'))
ggplot(anno_freq_summary_nk, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_nk <- as.SingleCellExperiment(seu_nk)
# i.nk <- seq_len(ncol(sce_nk))
# boot.nk1 <- sample(i.nk, replace=TRUE)
# boot.nk2 <- sample(i.nk, replace=TRUE)
# sce_nk_rep2 <- sce_nk[,boot.nk1]
# sce_nk_rep3 <- sce_nk[,boot.nk2]
# 
# sample <- c(sce_nk$age_rep,
#             paste0(sce_nk_rep2$age_rep,'_boot2'),
#             paste0(sce_nk_rep3$age_rep,'_boot3'))
# cluster <- c(sce_nk$subtype,
#              sce_nk_rep2$subtype,
#              sce_nk_rep3$subtype)
# group <- c(sce_nk$clade,
#            sce_nk_rep2$clade,
#            sce_nk_rep3$clade)
# allcounts <- cbind(counts(sce_nk),
#                    counts(sce_nk_rep2),
#                    counts(sce_nk_rep3))
# 
# sce_nk_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_nk_boot$sample <- sample
# sce_nk_boot$group <- as.character(group)
# sce_nk_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_nk_boot$group)
# ind_early <- grep('early',sce_nk_boot$group)
# ind_mid <- grep('mid',sce_nk_boot$group)
# ind_late <- grep('late',sce_nk_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_nk_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_nk_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_nk_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_nk_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_nk_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_nk_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_nk_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_nk_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_nk_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_nk_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_nk_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_nk_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_nk_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_nk_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_nk_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_nk_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_nk_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_nk_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_nk <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_nk) <- c('comparison','celltype','pval')
# prop_test_nk$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_nk$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_nk$pval <- pval_vec
# prop_test_nk$pval_bin <- 'ns'
# prop_test_nk$pval_bin[which(prop_test_nk$pval <= 0.05)] <- '*'

## Fig. 6B: Flow cytometry analysis of NK cell subtype proportions
# Proportion bar charts
ggplot(nk_facs_data_summary, aes(x=celltype, y=proportion, fill=clade)) + 
  geom_boxplot(position = position_dodge(0.8), color='black') + 
  geom_jitter(position = position_dodge(0.8)) + 
  theme_classic() + 
  scale_fill_manual(values=alpha(c('lightcoral','red','maroon'),0.75)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  ylim(c(0,80))

# Statistics
# nk_facs_pval_summary <- as.data.frame(matrix(0L, nrow=9, ncol=3))
# colnames(nk_facs_pval_summary) <- c('time','subtype','pval')
# nk_facs_pval_summary$time <- rep(c('early_mid','early_late','mid_late'),3)
# nk_facs_pval_summary$subtype <- rep(c('mod','int','cyto'), each=3)
# for (i in 1:nrow(nk_facs_pval_summary)) {
#   temp <- unlist(strsplit(nk_facs_pval_summary$time[i], split='_'))
#   x <- temp[1]
#   y <- temp[2]
#   ct <- nk_facs_pval_summary$subtype[i]
#   pval <- t.test(nk_facs_data_summary$proportion[which(nk_facs_data_summary$celltype == ct & nk_facs_data_summary$clade == x)],
#                  nk_facs_data_summary$proportion[which(nk_facs_data_summary$celltype == ct & nk_facs_data_summary$clade == y)])
#   nk_facs_pval_summary$pval[i] <- pval$p.value
# }
# nk_facs_pval_summary[,'pval_bin'] <- 'ns'
# nk_facs_pval_summary$pval_bin[which(nk_facs_pval_summary$pval <= 0.05)] <- '*'
# nk_facs_pval_summary$pval_bin <- factor(nk_facs_pval_summary$pval_bin, levels=c('ns','*'))
# nk_facs_pval_summary$subtype <- factor(nk_facs_pval_summary$subtype, levels=c('mod','int','cyto'))

## Fig. 6C: Kim et al human NK cell analysis
# GEX space
DimPlot(seu_nk_kim2020, group.by = 'subtype', cols=c('forestgreen','palegreen3'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_nk_kim2020, group.by = 'source', cols=c('darkred','lightcoral','dodgerblue','navy'), pt.size=0.25)+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_nk_kim2020, 'FCGR3A', max.cutoff = 'q95', pt.size = 0.25) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_nk_kim2020, 'GZMK', max.cutoff = 'q95', pt.size = 0.25) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

# Subtype proportions
# mLN_group <- c('BRONCHO_11',paste0('EBUS_',c(10,12,13,15,19,51)))
# nLN_group <- unique(seu_nk_kim2020@meta.data$sample[grep('LN_',seu_nk_kim2020@meta.data$sample)])
# tLung_group <- unique(seu_nk_kim2020@meta.data$sample[grep('LUNG_T',seu_nk_kim2020@meta.data$sample)])
# nLung_group <- unique(seu_nk_kim2020@meta.data$sample[grep('LUNG_N',seu_nk_kim2020@meta.data$sample)])
# 
# anno_freq <- table(seu_nk_kim2020@meta.data$sample, seu_nk_kim2020@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$source <- 'unknown'
# anno_freq$source[which(anno_freq$Var1 %in% mLN_group)] <- 'mLN'
# anno_freq$source[which(anno_freq$Var1 %in% nLN_group)] <- 'nLN'
# anno_freq$source[which(anno_freq$Var1 %in% tLung_group)] <- 'tLung'
# anno_freq$source[which(anno_freq$Var1 %in% nLung_group)] <- 'nLung'
# 
# anno_freq_summary_nk_kim <- as.data.frame(matrix(0L, nrow=4*2, ncol=4))
# colnames(anno_freq_summary_nk_kim) <- c('source','subtype','mean','se')
# anno_freq_summary_nk_kim$subtype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_nk_kim$source <- rep(unique(anno_freq$source), 2)
# for (i in 1:nrow(anno_freq_summary_nk_kim)) {
#   ct <- anno_freq_summary_nk_kim$subtype[i]
#   cl <- anno_freq_summary_nk_kim$source[i]
#   anno_freq_summary_nk_kim$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$source == cl)])
#   anno_freq_summary_nk_kim$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$source == cl)])
# }
# anno_freq_summary_nk_kim$source <- factor(anno_freq_summary_nk_kim$source, levels=c('nLN','mLN','nLung','tLung'))
# anno_freq_summary_nk_kim$subtype <- factor(anno_freq_summary_nk_kim$subtype, levels=c('mod','cyto'))
ggplot(anno_freq_summary_nk_kim[which(anno_freq_summary_nk_kim$subtype == 'cyto'), ], aes(x=subtype, y=mean, fill=source)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('lightcoral','darkred','dodgerblue','navy'))

## Fig. 6D: T cell subtype annotations and subtype proportion barcharts
# GEX space
DimPlot(seu_t, group.by = 'subtype', cols = c('olivedrab3','darkseagreen4','navy','steelblue3','deepskyblue3','grey','black','forestgreen')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
  
# Subtype proportions
# anno_freq <- table(seu_t@meta.data$age_rep, seu_t@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_t <- as.data.frame(matrix(0L, nrow=4*8, ncol=4))
# colnames(anno_freq_summary_t) <- c('celltype','clade','mean','se')
# anno_freq_summary_t$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_t$clade <- rep(unique(anno_freq$clade), 8)
# for (i in 1:nrow(anno_freq_summary_t)) {
#   ct <- anno_freq_summary_t$celltype[i]
#   cl <- anno_freq_summary_t$clade[i]
#   anno_freq_summary_t$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_t$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_t$clade <- factor(anno_freq_summary_t$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_t$celltype <- factor(anno_freq_summary_t$celltype, levels=c('cd4_naive','cd8_naive','cd8_memory','cd8_effector','treg','cd4_cd29','gdT_ILC','prolif'))
ggplot(anno_freq_summary_t, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_t <- as.SingleCellExperiment(seu_t)
# i.t <- seq_len(ncol(sce_t))
# boot.t1 <- sample(i.t, replace=TRUE)
# boot.t2 <- sample(i.t, replace=TRUE)
# sce_t_rep2 <- sce_t[,boot.t1]
# sce_t_rep3 <- sce_t[,boot.t2]
# 
# sample <- c(sce_t$age_rep,
#             paste0(sce_t_rep2$age_rep,'_boot2'),
#             paste0(sce_t_rep3$age_rep,'_boot3'))
# cluster <- c(sce_t$subtype,
#              sce_t_rep2$subtype,
#              sce_t_rep3$subtype)
# group <- c(sce_t$clade,
#            sce_t_rep2$clade,
#            sce_t_rep3$clade)
# allcounts <- cbind(counts(sce_t),
#                    counts(sce_t_rep2),
#                    counts(sce_t_rep3))
# 
# sce_t_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_t_boot$sample <- sample
# sce_t_boot$group <- as.character(group)
# sce_t_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_t_boot$group)
# ind_early <- grep('early',sce_t_boot$group)
# ind_mid <- grep('mid',sce_t_boot$group)
# ind_late <- grep('late',sce_t_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_t_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_t_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_t_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_t_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_t_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_t_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_t_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_t_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_t_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_t_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_t_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_t_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_t_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_t_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_t_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_t_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_t_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_t_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_t <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_t) <- c('comparison','celltype','pval')
# prop_test_t$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_t$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_t$pval <- pval_vec
# prop_test_t$pval_bin <- 'ns'
# prop_test_t$pval_bin[which(prop_test_t$pval <= 0.05)] <- '*'

## Fig. 7E: Z-score heatmap of T-cell functional markers
# t_cell_genes <- c('Sell','Tcf7','Lef1','Foxp3','Entpd1','Nt5e','Il10','Tbx21','Ifng','Rorc','Il17a','Mki67','Prf1','Gzma','Gzmb','Gzmk','Cd27','Cd28','Tnfrsf9','Tnfrsf4','Icos','Btla','Pdcd1','Lag3','Havcr2','Tigit','Ctla4','Vsir','Tnfrsf18')
# t_cell_zscores <- as.data.frame(matrix(0L, nrow=8, ncol=length(t_cell_genes)))
# rownames(t_cell_zscores) <- c('cd4_naive','cd8_naive','cd8_memory','cd8_effector','treg','cd4_cd29','gdT_ILC','prolif')
# colnames(t_cell_zscores) <- t_cell_genes
# t_cell_genes_mean <- rowMeans(seu_t@assays$RNA@data[t_cell_genes, ])
# t_cell_genes_sd <- apply(seu_t@assays$RNA@data[t_cell_genes, ], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(t_cell_zscores)) {
#   ind <- which(seu_t@meta.data$subtype == i)
#   temp <-  rowMeans(seu_t@assays$RNA@data[t_cell_genes,ind])
#   t_cell_zscores[i, ] <- (temp-t_cell_genes_mean)/t_cell_genes_sd
# }
Heatmap(t_cell_zscores[,-which(colnames(t_cell_zscores) == 'Mki67')], cluster_columns = F)

## Fig. 7F: B cell subtype annotations and subtype proportion barcharts
DimPlot(seu_b, group.by = 'subtype', cols=c('forestgreen','olivedrab3','darkseagreen4')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

# Subtype proportions
# anno_freq <- table(seu_b@meta.data$age_rep, seu_b@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_b <- as.data.frame(matrix(0L, nrow=4*3, ncol=4))
# colnames(anno_freq_summary_b) <- c('celltype','clade','mean','se')
# anno_freq_summary_b$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_b$clade <- rep(unique(anno_freq$clade), 3)
# for (i in 1:nrow(anno_freq_summary_b)) {
#   ct <- anno_freq_summary_b$celltype[i]
#   cl <- anno_freq_summary_b$clade[i]
#   anno_freq_summary_b$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_b$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_b$clade <- factor(anno_freq_summary_b$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_b$celltype <- factor(anno_freq_summary_b$celltype, levels=c('naive','ighm','breg'))
ggplot(anno_freq_summary_b, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_b <- as.SingleCellExperiment(seu_b)
# i.b <- seq_len(ncol(sce_b))
# boot.b1 <- sample(i.b, replace=TRUE)
# boot.b2 <- sample(i.b, replace=TRUE)
# sce_b_rep2 <- sce_b[,boot.b1]
# sce_b_rep3 <- sce_b[,boot.b2]
# 
# sample <- c(sce_b$age_rep,
#             paste0(sce_b_rep2$age_rep,'_boot2'),
#             paste0(sce_b_rep3$age_rep,'_boot3'))
# cluster <- c(sce_b$subtype,
#              sce_b_rep2$subtype,
#              sce_b_rep3$subtype)
# group <- c(sce_b$clade,
#            sce_b_rep2$clade,
#            sce_b_rep3$clade)
# allcounts <- cbind(counts(sce_b),
#                    counts(sce_b_rep2),
#                    counts(sce_b_rep3))
# 
# sce_b_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_b_boot$sample <- sample
# sce_b_boot$group <- as.character(group)
# sce_b_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_b_boot$group)
# ind_early <- grep('early',sce_b_boot$group)
# ind_mid <- grep('mid',sce_b_boot$group)
# ind_late <- grep('late',sce_b_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_b_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_b_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_b_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_b_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_b_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_b_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_b_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_b_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_b_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_b_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_b_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_b_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_b_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_b_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_b_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_b_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_b_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_b_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_b <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_b) <- c('comparison','celltype','pval')
# prop_test_b$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_b$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_b$pval <- pval_vec
# prop_test_b$pval_bin <- 'ns'
# prop_test_b$pval_bin[which(prop_test_b$pval <= 0.05)] <- '*'
