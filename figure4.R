#######################################################################
## The temporal progression of immune remodeling during metastasis  ###
## Companion code: Figure 4, Chris McGinnis, PhD, Stanford, 05/2023  ##
#######################################################################
## objects
load('seu_mono.Robj')
load('seu_neu.Robj')
load('seu_dc.Robj')
load('anno_freq_summary_mono.Robj')
load('anno_freq_summary_neu.Robj')
load('anno_freq_summary_dc.Robj')
load('prop_test_mono.Robj')
load('prop_test_neu.Robj')
load('prop_test_dc.Robj')
load("seu_neu_mature.Robj")
load('umap_neu_mature.Robj')
load('mature_neu_zscores.Robj')

## functions
load("stderror.Robj")

######################################################################################################
## Figure 4: BM-derived myeloid subtype characterization uncovers metastasis-associated DC subtype  ##
## frequency shifts and induction of TLR-NFκB inflammation in pre-metastatic niche neutrophils  ######
######################################################################################################
## Fig. 4A: Monocyte subtype annotations and subtype proportion barcharts
# GEX Space
DimPlot(seu_mono, group.by = 'subtype', cols = c('dodgerblue','steelblue3','navy')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

# Subtype proportions
# anno_freq <- table(seu_mono@meta.data$age_rep, seu_mono@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_mono <- as.data.frame(matrix(0L, nrow=4*3, ncol=4))
# colnames(anno_freq_summary_mono) <- c('celltype','clade','mean','se')
# anno_freq_summary_mono$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_mono$clade <- rep(unique(anno_freq$clade), 3)
# for (i in 1:nrow(anno_freq_summary_mono)) {
#   ct <- anno_freq_summary_mono$celltype[i]
#   cl <- anno_freq_summary_mono$clade[i]
#   anno_freq_summary_mono$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_mono$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_mono$clade <- factor(anno_freq_summary_mono$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_mono$celltype <- factor(anno_freq_summary_mono$celltype, levels=c('cm','intm','ncm'))
ggplot(anno_freq_summary_mono, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_mono <- as.SingleCellExperiment(seu_mono)
# i.mono <- seq_len(ncol(sce_mono))
# boot.mono1 <- sample(i.mono, replace=TRUE)
# boot.mono2 <- sample(i.mono, replace=TRUE)
# sce_mono_rep2 <- sce_mono[,boot.mono1]
# sce_mono_rep3 <- sce_mono[,boot.mono2]
# 
# sample <- c(sce_mono$age_rep,
#             paste0(sce_mono_rep2$age_rep,'_boot2'),
#             paste0(sce_mono_rep3$age_rep,'_boot3'))
# cluster <- c(sce_mono$subtype,
#              sce_mono_rep2$subtype,
#              sce_mono_rep3$subtype)
# group <- c(sce_mono$clade,
#            sce_mono_rep2$clade,
#            sce_mono_rep3$clade)
# allcounts <- cbind(counts(sce_mono),
#                    counts(sce_mono_rep2),
#                    counts(sce_mono_rep3))
# 
# sce_mono_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_mono_boot$sample <- sample
# sce_mono_boot$group <- as.character(group)
# sce_mono_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_mono_boot$group)
# ind_early <- grep('early',sce_mono_boot$group)
# ind_mid <- grep('mid',sce_mono_boot$group)
# ind_late <- grep('late',sce_mono_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_mono_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_mono_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_mono_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_mono_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_mono_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_mono_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_mono_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_mono_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_mono_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_mono_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_mono_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_mono_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_mono_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_mono_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_mono_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_mono_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_mono_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_mono_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_mono <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_mono) <- c('comparison','celltype','pval')
# prop_test_mono$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_mono$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_mono$pval <- pval_vec
# prop_test_mono$pval_bin <- 'ns'
# prop_test_mono$pval_bin[which(prop_test_mono$pval <= 0.05)] <- '*'

## Fig. 4B: Neutrophil subtype annotations, subtype proportion barcharts, and Cd14+ inflammatory signature feature plots
# GEX Space
DimPlot(seu_neu, group.by = 'subtype', cols=c('navy','red','dodgerblue','steelblue3'), pt.size = 0.1) + NoAxes() + theme(plot.title = element_blank()) + NoLegend()

# act_genes <- c('Cd14','Ccl4', 'Ccl3', 'Cxcl2', 'Cxcl3', 'Spp1', 'Il1b', 'Nfkbia', 'Socs3', 'Mif', 'Klf6', 'Atf3', 'Ptgs2', 'Xbp1','Jun','Ccrl2','Saa3','Gad45b','Ninj1','Clec4n','Hcar2','Basp1','Btg1','Il1rn','Ifrd1','Txnip','Ccl9','Ier3','Ier5','Rs1','Thbs1','Cxcl1','Hilpda','Hist1h1c','Srgn','Hspa5','Csf1','Pgts2')
# seu_neu <- AddModuleScore(seu_neu, features = list(act_genes[which(act_genes %in% seu_neu@assays$RNA@var.features)]), ctrl=5, name='mdsc')
FeaturePlot(seu_neu, 'Cxcl2', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()
FeaturePlot(seu_neu, 'Cd14', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()
FeaturePlot(seu_neu, 'Nlrp3', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()
FeaturePlot(seu_neu, 'mdsc1', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()

# Subtype proportions
# anno_freq <- table(seu_neu@meta.data$age_rep, seu_neu@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_neu <- as.data.frame(matrix(0L, nrow=4*4, ncol=4))
# colnames(anno_freq_summary_neu) <- c('celltype','clade','mean','se')
# anno_freq_summary_neu$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_neu$clade <- rep(unique(anno_freq$clade), 4)
# for (i in 1:nrow(anno_freq_summary_neu)) {
#   ct <- anno_freq_summary_neu$celltype[i]
#   cl <- anno_freq_summary_neu$clade[i]
#   anno_freq_summary_neu$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_neu$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_neu$clade <- factor(anno_freq_summary_neu$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_neu$celltype <- factor(anno_freq_summary_neu$celltype, levels=c('mature','transition','immature','inf'))
ggplot(anno_freq_summary_neu, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_neu <- as.SingleCellExperiment(seu_neu)
# i.neu <- seq_len(ncol(sce_neu))
# boot.neu1 <- sample(i.neu, replace=TRUE)
# boot.neu2 <- sample(i.neu, replace=TRUE)
# sce_neu_rep2 <- sce_neu[,boot.neu1]
# sce_neu_rep3 <- sce_neu[,boot.neu2]
# 
# sample <- c(sce_neu$age_rep,
#             paste0(sce_neu_rep2$age_rep,'_boot2'),
#             paste0(sce_neu_rep3$age_rep,'_boot3'))
# cluster <- c(sce_neu$subtype_new,
#              sce_neu_rep2$subtype_new,
#              sce_neu_rep3$subtype_new)
# group <- c(sce_neu$clade,
#            sce_neu_rep2$clade,
#            sce_neu_rep3$clade)
# allcounts <- cbind(counts(sce_neu),
#                    counts(sce_neu_rep2),
#                    counts(sce_neu_rep3))
# 
# sce_neu_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_neu_boot$sample <- sample
# sce_neu_boot$group <- as.character(group)
# sce_neu_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_neu_boot$group)
# ind_early <- grep('early',sce_neu_boot$group)
# ind_mid <- grep('mid',sce_neu_boot$group)
# ind_late <- grep('late',sce_neu_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_neu_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_neu_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_neu_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_neu_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_neu_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_neu_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_neu_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_neu_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_neu_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_neu_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_neu_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_neu_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_neu_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_neu_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_neu_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_neu_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_neu_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_neu_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_neu <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_neu) <- c('comparison','celltype','pval')
# prop_test_neu$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_neu$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_neu$pval <- pval_vec
# prop_test_neu <- prop_test_neu[which(prop_test_neu$celltype %in% c('mature','transitional','immature','inf')), ]
# prop_test_neu$pval_bin <- 'ns'
# prop_test_neu$pval_bin[which(prop_test_neu$pval <= 0.05)] <- '*'

## Fig. 4C: Mature neutrophil metastatic stage UMAPs and clade-specific zscore heatmap
# umap_neu_mature <- as.data.frame(seu_neu_mature@reductions$umap@cell.embeddings)
# umap_neu_mature[,'clade'] <- seu_neu_mature@meta.data$clade
ggplot(umap_neu_mature, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature[which(umap_neu_mature$clade == 'wt'), ], color = 'black') +
  theme_void() 
ggplot(umap_neu_mature, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature[which(umap_neu_mature$clade == 'early'), ], color = 'lightcoral') +
  theme_void() 
ggplot(umap_neu_mature, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature[which(umap_neu_mature$clade == 'mid'), ], color = 'red') +
  theme_void() 
ggplot(umap_neu_mature, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature[which(umap_neu_mature$clade == 'late'), ], color = 'maroon') +
  theme_void() 

# mature_neu_stage_markers <- FindAllMarkers(seu_neu_mature, only.pos = T, logfc.threshold = log(1.75), min.pct = 0.5)
# mature_neu_genes <- c(mature_neu_stage_markers$gene[which(mature_neu_stage_markers$cluster == 'wt')][1:10],
#                       mature_neu_stage_markers$gene[which(mature_neu_stage_markers$cluster == 'early')][1:10],
#                       mature_neu_stage_markers$gene[which(mature_neu_stage_markers$cluster == 'late')])
# mature_neu_genes <- unique(mature_neu_genes)
# mature_neu_zscores <- as.data.frame(matrix(0L, nrow=4, ncol=length(mature_neu_genes)))
# rownames(mature_neu_zscores) <- c('wt','early','mid','late')
# colnames(mature_neu_zscores) <- mature_neu_genes
# mature_neu_genes_mean <- rowMeans(seu_neu_mature@assays$RNA@data[mature_neu_genes, ])
# mature_neu_genes_sd <- apply(seu_neu_mature@assays$RNA@data[mature_neu_genes, ], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(mature_neu_zscores)) {
#   ind <- which(seu_neu_mature@meta.data$clade == i)
#   temp <-  rowMeans(seu_neu_mature@assays$RNA@data[mature_neu_genes,ind])
#   mature_neu_zscores[i, ] <- (temp-mature_neu_genes_mean)/mature_neu_genes_sd
# }
Heatmap(mature_neu_zscores, cluster_rows = F, show_heatmap_legend = F)

## Fig. 4D: DC subtype annotations, metastatic stage density UMAPs, and subtype proportion barcharts
# GEX Space
DimPlot(seu_dc, group.by = 'subtype', cols=c('tan4','deepskyblue3','steelblue3','navy','black','grey')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

# Subtype proportions
# anno_freq <- table(seu_dc@meta.data$age_rep, seu_dc@meta.data$subtype)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_dc <- as.data.frame(matrix(0L, nrow=4*6, ncol=4))
# colnames(anno_freq_summary_dc) <- c('celltype','clade','mean','se')
# anno_freq_summary_dc$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_dc$clade <- rep(unique(anno_freq$clade), 6)
# for (i in 1:nrow(anno_freq_summary_dc)) {
#   ct <- anno_freq_summary_dc$celltype[i]
#   cl <- anno_freq_summary_dc$clade[i]
#   anno_freq_summary_dc$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_dc$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# 
# anno_freq_summary_dc$clade <- factor(anno_freq_summary_dc$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_dc$celltype <- factor(anno_freq_summary_dc$celltype, levels=c('cdc1','cdc2','IFN','ccr7','pdc','prolif'))
ggplot(anno_freq_summary_dc, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey','lightcoral','red','maroon'))

# Statistics
# sce_dc <- as.SingleCellExperiment(seu_dc)
# i.dc <- seq_len(ncol(sce_dc))
# boot.dc1 <- sample(i.dc, replace=TRUE)
# boot.dc2 <- sample(i.dc, replace=TRUE)
# sce_dc_rep2 <- sce_dc[,boot.dc1]
# sce_dc_rep3 <- sce_dc[,boot.dc2]
# 
# sample <- c(sce_dc$age_rep,
#             paste0(sce_dc_rep2$age_rep,'_boot2'),
#             paste0(sce_dc_rep3$age_rep,'_boot3'))
# cluster <- c(sce_dc$subtype,
#              sce_dc_rep2$subtype,
#              sce_dc_rep3$subtype)
# group <- c(sce_dc$clade,
#            sce_dc_rep2$clade,
#            sce_dc_rep3$clade)
# allcounts <- cbind(counts(sce_dc),
#                    counts(sce_dc_rep2),
#                    counts(sce_dc_rep3))
# 
# sce_dc_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_dc_boot$sample <- sample
# sce_dc_boot$group <- as.character(group)
# sce_dc_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_dc_boot$group)
# ind_early <- grep('early',sce_dc_boot$group)
# ind_mid <- grep('mid',sce_dc_boot$group)
# ind_late <- grep('late',sce_dc_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_dc_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_dc_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_dc_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_dc_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_dc_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_dc_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_dc_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_dc_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_dc_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_dc_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_dc_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_dc_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_dc_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_dc_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_dc_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_dc_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_dc_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_dc_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_dc <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_dc) <- c('comparison','celltype','pval')
# prop_test_dc$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_dc$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_dc$pval <- pval_vec
# prop_test_dc$pval_bin <- 'ns'
# prop_test_dc$pval_bin[which(prop_test_dc$pval <= 0.05)] <- '*'