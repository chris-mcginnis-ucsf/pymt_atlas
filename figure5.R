#######################################################################
## The temporal progression of immune remodeling during metastasis  ###
## Companion code: Figure 5, Chris McGinnis, PhD, Stanford, 05/2023  ##
#######################################################################
## objects
load('nmf_myl.Robj')
load('seu_mono.Robj')
load('seu_int.Robj')
load('seu_mam_sub.Robj')
load('anno_freq_summary_mono_inf.Robj')
load('anno_freq_summary_int_inf.Robj')
load('prop_test_mono_inf.Robj')
load('prop_test_int_inf.Robj')
load('gsea_myl_inf.Robj')
load('myl_facs_data_summary.Robj')
load('myl_facs_data_raw.Robj')
load('seu_fresh_mono.Robj')
load('seu_fresh_neu.Robj')
load('seu_4t1_mono.Robj')
load('seu_4t1_neu.Robj')

##############################################################################################################################################################
## Figure 5: TLR-NFκB inflammation signature is present in BM-derived and tissue-resident myeloid cells and correlates with pre-metastatic niche formation  ##
##############################################################################################################################################################
## Fig. 5A: Monocyte and IM inflammatory subtype annotations with NMF30 and 'activated' MDSC signature feature plots
# act_genes <- c('Cd14','Ccl4', 'Ccl3', 'Cxcl2', 'Cxcl3', 'Spp1', 'Il1b', 'Nfkbia', 'Socs3', 'Mif', 'Klf6', 'Atf3', 'Ptgs2', 'Xbp1','Jun','Ccrl2','Saa3','Gad45b','Ninj1','Clec4n','Hcar2','Basp1','Btg1','Il1rn','Ifrd1','Txnip','Ccl9','Ier3','Ier5','Rs1','Thbs1','Cxcl1','Hilpda','Hist1h1c','Srgn','Hspa5','Csf1','Pgts2')
# seu_mono <- AddModuleScore(seu_mono, features = list(act_genes[which(act_genes %in% seu_mono@assays$RNA@var.features)]), ctrl=5, name='mdsc')
# seu_int <- AddModuleScore(seu_int, features = list(act_genes[which(act_genes %in% seu_int@assays$RNA@var.features)]), ctrl=5, name='mdsc')
# seu_mono@meta.data[,'nmf19'] <- nmf_myl$h[19,rownames(seu_mono@meta.data)]
# seu_int@meta.data[,'nmf19'] <- nmf_myl$h[19,rownames(seu_int@meta.data)]
# nmf19_genes <- rownames(nmf_myl$w)[which(nmf_myl$w[,19] >= 0.005)]
# seu_mono <- AddModuleScore(seu_mono, features = list(nmf19_genes), ctrl=5, name='nmf19_mod')
# seu_int <- AddModuleScore(seu_int, features = list(nmf19_genes), ctrl=5, name='nmf19_mod')

FeaturePlot(seu_mono, 'nmf19_mod1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
DimPlot(seu_mono, group.by = 'subtype_inf', cols=c('dodgerblue','red','steelblue3','navy'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'mdsc1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

FeaturePlot(seu_int, 'nmf19_mod1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
DimPlot(seu_int, group.by = 'subtype_inf', cols=c('darkred','goldenrod','red','black','grey'), pt.size = 1)  + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'mdsc1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. 5B: Monocyte and IM inflammatory subtype proportion bar charts
# Mono
# anno_freq <- table(seu_mono@meta.data$age_rep, seu_mono@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_mono_inf <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_mono_inf) <- c('celltype','clade','mean','se')
# anno_freq_summary_mono_inf$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_mono_inf$clade <- rep(unique(anno_freq$clade), 2)
# for (i in 1:nrow(anno_freq_summary_mono_inf)) {
#   ct <- anno_freq_summary_mono_inf$celltype[i]
#   cl <- anno_freq_summary_mono_inf$clade[i]
#   anno_freq_summary_mono_inf$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_mono_inf$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# anno_freq_summary_mono_inf$clade <- factor(anno_freq_summary_mono_inf$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_mono_inf$celltype <- factor(anno_freq_summary_mono_inf$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_mono_inf[which(anno_freq_summary_mono_inf$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

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
# cluster <- c(sce_mono$inf,
#              sce_mono_rep2$inf,
#              sce_mono_rep3$inf)
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
# prop_test_mono_inf <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_mono_inf) <- c('comparison','celltype','pval')
# prop_test_mono_inf$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_mono_inf$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_mono_inf$pval <- pval_vec
# prop_test_mono_inf$pval_bin <- 'ns'
# prop_test_mono_inf$pval_bin[which(prop_test_mono_inf$pval <= 0.05)] <- '*'

# IM
# anno_freq <- table(seu_int@meta.data$age_rep, seu_int@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_int_inf <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_int_inf) <- c('celltype','clade','mean','se')
# anno_freq_summary_int_inf$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_int_inf$clade <- rep(unique(anno_freq$clade), 2)
# for (i in 1:nrow(anno_freq_summary_int_inf)) {
#   ct <- anno_freq_summary_int_inf$celltype[i]
#   cl <- anno_freq_summary_int_inf$clade[i]
#   anno_freq_summary_int_inf$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_int_inf$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# anno_freq_summary_int_inf$clade <- factor(anno_freq_summary_int_inf$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_int_inf$celltype <- factor(anno_freq_summary_int_inf$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_int_inf[which(anno_freq_summary_int_inf$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=clade)) + 
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
# cluster <- c(sce_int$inf,
#              sce_int_rep2$inf,
#              sce_int_rep3$inf)
# group <- c(sce_int$clade,
#            sce_int_rep2$clade,
#            sce_int_rep3$clade)
# allcounts <- cbind(counts(sce_int),
#                    counts(sce_int_rep2),
#                    counts(sce_int_rep3))
# 
# sce_int_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_int_boot$sample <- sample
# sce_int_boot$group <- as.character(group)
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
# prop_test_int_inf <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_int_inf) <- c('comparison','celltype','pval')
# prop_test_int_inf$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_int_inf$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_int_inf$pval <- pval_vec
# prop_test_int_inf$pval_bin <- 'ns'
# prop_test_int_inf$pval_bin[which(prop_test_int_inf$pval <= 0.05)] <- '*'

## Fig. 5C: 'Activated' MDSC vs NMF20 module score scatter plots
ggplot(seu_mono@meta.data, aes(x=nmf19_mod1, y=mdsc1)) + 
  theme_classic() + 
  geom_smooth(method=lm, color='black') + 
  geom_point(data=seu_mono@meta.data, aes(x=nmf19_mod1, y=mdsc1, color=inf), size=0.01) + 
  scale_color_manual(values=c('red','black')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + ylim(c(-0.6, 1))
cor(seu_mono@meta.data$nmf19_mod1, seu_mono@meta.data$mdsc1)

ggplot(seu_int@meta.data, aes(x=nmf19_mod1, y=mdsc1)) + 
  theme_classic() + 
  geom_smooth(method=lm, color='black') + 
  geom_point(data=seu_int@meta.data, aes(x=nmf19_mod1, y=mdsc1, color=inf), size=0.5) + 
  scale_color_manual(values=c('red','black')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + ylim(c(-0.8, 0.8))
cor(seu_int@meta.data$nmf19_mod1, seu_int@meta.data$mdsc1)

## Fig. 5D: NMF39 Hallmark GSEA lollipop plot 
# hall <- msigdb_gsets(species="Mus musculus", category="H")
# gsea_hall_nmf19 <- hypeR(signature = nmf19_genes, genesets = hall)
# gsea_myl_inf <- gsea_hall_nmf19$data[order(gsea_hall_nmf19$data$pval, decreasing = F), ]
ggplot(gsea_myl_inf[which(gsea_myl_inf$pval <= 0.01), ], aes(x=reorder(label, -log10(pval)), y = -log10(pval))) +
  geom_point(stat="identity",color='black') + geom_segment(aes(x=label, xend=label, y=0, yend=-log10(pval))) + geom_hline(yintercept = 2, linetype="dotted") + 
  coord_flip() + xlab('Pathway') + theme_classic() + theme(axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text.y = element_blank())

## Fig. 5E: Flow cytometry analysis of CD14 abundance on neutrophils, AMs, and monocytes
# Geometric mean bar charts
ggplot(myl_facs_data_summary[which(myl_facs_data_summary$celltype == 'neu'), ], aes(x=time, y=GEOMEAN, fill=time)) + 
  geom_boxplot(position = position_dodge(0.8), color='black') + 
  geom_jitter(position = position_dodge(0.8)) + 
  theme_classic() + 
  scale_fill_manual(values=alpha(c('lightcoral','red','maroon'),0.75)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  ylim(c(0,600))

ggplot(myl_facs_data_summary[which(myl_facs_data_summary$celltype == 'alv'), ], aes(x=time, y=GEOMEAN, fill=time)) + 
  geom_boxplot(position = position_dodge(0.8), color='black') + 
  geom_jitter(position = position_dodge(0.8)) + 
  theme_classic() + 
  scale_fill_manual(values=alpha(c('lightcoral','red','maroon'),0.75)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  ylim(c(0,600))

ggplot(myl_facs_data_summary[which(myl_facs_data_summary$celltype == 'mono'), ], aes(x=time, y=GEOMEAN, fill=time)) + 
  geom_boxplot(position = position_dodge(0.8), color='black') + 
  geom_jitter(position = position_dodge(0.8)) + 
  theme_classic() + 
  scale_fill_manual(values=alpha(c('lightcoral','red','maroon'),0.75)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  ylim(c(0,600))

# Density plots
ggplot(myl_facs_data_raw[which(myl_facs_data_raw$celltype == 'neu'), ], aes(x = log10(cd14), color = time)) + 
  geom_density(alpha=0, size=1.25) + 
  theme_classic() + 
  scale_color_manual(values=c('lightcoral','red','maroon')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  xlim(c(2.4,3))

ggplot(myl_facs_data_raw[which(myl_facs_data_raw$celltype == 'alv'), ], aes(x = log10(cd14), color = time)) + 
  geom_density(alpha=0, size=1.25) + 
  theme_classic() + 
  scale_color_manual(values=c('lightcoral','red','maroon')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  xlim(c(2.6,3))

ggplot(myl_facs_data_raw[which(myl_facs_data_raw$celltype == 'mono'), ], aes(x = log10(cd14), color = time)) + 
  geom_density(alpha=0, size=1.25) + 
  theme_classic() + 
  scale_color_manual(values=c('lightcoral','red','maroon')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) 

# Statistics
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'early' & myl_facs_data_raw$celltype == 'alv')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'mid' & myl_facs_data_raw$celltype == 'alv')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'early' & myl_facs_data_raw$celltype == 'alv')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'late' & myl_facs_data_raw$celltype == 'alv')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'mid' & myl_facs_data_raw$celltype == 'alv')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'late' & myl_facs_data_raw$celltype == 'alv')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'early' & myl_facs_data_raw$celltype == 'neu')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'mid' & myl_facs_data_raw$celltype == 'neu')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'early' & myl_facs_data_raw$celltype == 'neu')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'late' & myl_facs_data_raw$celltype == 'neu')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'mid' & myl_facs_data_raw$celltype == 'neu')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'late' & myl_facs_data_raw$celltype == 'neu')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'early' & myl_facs_data_raw$celltype == 'mono')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'mid' & myl_facs_data_raw$celltype == 'mono')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'early' & myl_facs_data_raw$celltype == 'mono')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'late' & myl_facs_data_raw$celltype == 'mono')], alternative = "two.sided")
# wilcox.test(x = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'mid' & myl_facs_data_raw$celltype == 'mono')], y = myl_facs_data_raw$cd14[which(myl_facs_data_raw$time == 'late' & myl_facs_data_raw$celltype == 'mono')], alternative = "two.sided")


## Fig. 5F: PyMT validation cohort monocyte and neutrophil annotation UMAPs and inflammatory subtype proportion bar charts
# seu_fresh_mono <- AddModuleScore(seu_fresh_mono, features = list(nmf19_genes), ctrl=5, name='nmf19_mod')
# seu_fresh_neu <- AddModuleScore(seu_fresh_neu, features = list(nmf19_genes), ctrl=5, name='nmf19_mod')
DimPlot(seu_fresh_mono, group.by = 'subtype', cols=c('dodgerblue','red','steelblue3','navy'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_fresh_mono, 'nmf19_mod1', max.cutoff = 'q95') + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()

DimPlot(seu_fresh_neu, group.by = 'subtype', cols=c('dodgerblue','red','navy','steelblue3'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_fresh_neu, 'nmf19_mod1', max.cutoff = 'q95') + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()

## Mono
# anno_freq <- table(seu_fresh_mono@meta.data$sample, seu_fresh_mono@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$group <- 'unknown'
# anno_freq$group[grep('wt', anno_freq$Var1)] <- 'wt'
# anno_freq$group[grep('pymt_early', anno_freq$Var1)] <- 'early'
# anno_freq$group[grep('pymt_mid', anno_freq$Var1)] <- 'mid'
# anno_freq$group[grep('pymt_late', anno_freq$Var1)] <- 'late'
# anno_freq_summary_mono_fresh <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_mono_fresh) <- c('celltype','group','mean','se')
# anno_freq_summary_mono_fresh$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_mono_fresh$group <- rep(unique(anno_freq$group), 2)
# for (i in 1:nrow(anno_freq_summary_mono_fresh)) {
#   ct <- anno_freq_summary_mono_fresh$celltype[i]
#   cl <- anno_freq_summary_mono_fresh$group[i]
#   anno_freq_summary_mono_fresh$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
#   anno_freq_summary_mono_fresh$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
# }
# anno_freq_summary_mono_fresh$group <- factor(anno_freq_summary_mono_fresh$group, levels=c('wt','early','mid','late'))
# anno_freq_summary_mono_fresh$celltype <- factor(anno_freq_summary_mono_fresh$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_mono_fresh[which(anno_freq_summary_mono_fresh$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics
# sce_mono <- as.SingleCellExperiment(seu_fresh_mono)
# i.mono <- seq_len(ncol(sce_mono))
# boot.mono1 <- sample(i.mono, replace=TRUE)
# boot.mono2 <- sample(i.mono, replace=TRUE)
# sce_mono_rep2 <- sce_mono[,boot.mono1]
# sce_mono_rep3 <- sce_mono[,boot.mono2]
# 
# sample <- c(sce_mono$sample,
#             paste0(sce_mono_rep2$sample,'_boot2'),
#             paste0(sce_mono_rep3$sample,'_boot3'))
# cluster <- c(sce_mono$inf,
#              sce_mono_rep2$inf,
#              sce_mono_rep3$inf)
# group <- c(sce_mono$group,
#            sce_mono_rep2$group,
#            sce_mono_rep3$group)
# allcounts <- cbind(counts(sce_mono),
#                    counts(sce_mono_rep2),
#                    counts(sce_mono_rep3))
# 
# sce_mono_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_mono_boot$sample <- sample
# sce_mono_boot$group <- group
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

## Neutrophil
# anno_freq <- table(seu_fresh_neu@meta.data$sample, seu_fresh_neu@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$group <- 'unknown'
# anno_freq$group[grep('wt', anno_freq$Var1)] <- 'wt'
# anno_freq$group[grep('pymt_early', anno_freq$Var1)] <- 'early'
# anno_freq$group[grep('pymt_mid', anno_freq$Var1)] <- 'mid'
# anno_freq$group[grep('pymt_late', anno_freq$Var1)] <- 'late'
# anno_freq_summary_neu_fresh <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_neu_fresh) <- c('celltype','group','mean','se')
# anno_freq_summary_neu_fresh$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_neu_fresh$group <- rep(unique(anno_freq$group), 2)
# for (i in 1:nrow(anno_freq_summary_neu_fresh)) {
#   ct <- anno_freq_summary_neu_fresh$celltype[i]
#   cl <- anno_freq_summary_neu_fresh$group[i]
#   anno_freq_summary_neu_fresh$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
#   anno_freq_summary_neu_fresh$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
# }
# anno_freq_summary_neu_fresh$group <- factor(anno_freq_summary_neu_fresh$group, levels=c('wt','early','mid','late'))
# anno_freq_summary_neu_fresh$celltype <- factor(anno_freq_summary_neu_fresh$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_neu_fresh[which(anno_freq_summary_neu_fresh$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics
# sce_neu <- as.SingleCellExperiment(seu_fresh_neu)
# i.neu <- seq_len(ncol(sce_neu))
# boot.neu1 <- sample(i.neu, replace=TRUE)
# boot.neu2 <- sample(i.neu, replace=TRUE)
# sce_neu_rep2 <- sce_neu[,boot.neu1]
# sce_neu_rep3 <- sce_neu[,boot.neu2]
# 
# sample <- c(sce_neu$sample,
#             paste0(sce_neu_rep2$sample,'_boot2'),
#             paste0(sce_neu_rep3$sample,'_boot3'))
# cluster <- c(sce_neu$inf,
#              sce_neu_rep2$inf,
#              sce_neu_rep3$inf)
# group <- c(sce_neu$group,
#            sce_neu_rep2$group,
#            sce_neu_rep3$group)
# allcounts <- cbind(counts(sce_neu),
#                    counts(sce_neu_rep2),
#                    counts(sce_neu_rep3))
# 
# sce_neu_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_neu_boot$sample <- sample
# sce_neu_boot$group <- group
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
# prop_test_neu$pval_bin <- 'ns'
# prop_test_neu$pval_bin[which(prop_test_neu$pval <= 0.05)] <- '*'

## Fig. 5G: 4T1 longitudinal atlas monocyte and basophil annotation UMAPs and inflammatory subtype proportion bar charts
# seu_4t1_mono_fin <- AddModuleScore(seu_4t1_mono_fin, features = list(nmf19_genes), ctrl=5, name='nmf19_mod')
# seu_4t1_baso_fin <- AddModuleScore(seu_4t1_baso_fin, features = list(nmf19_genes), ctrl=5, name='nmf19_mod')
DimPlot(seu_4t1_mono_fin, group.by = 'subtype', cols=c('dodgerblue','red','steelblue3','navy'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_4t1_mono_fin, 'nmf19_mod1', max.cutoff = 'q95') + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()

DimPlot(seu_4t1_baso_fin, group.by = 'subtype', cols=c('dodgerblue','steelblue3','navy','red','tan4'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_4t1_baso_fin, 'nmf19_mod1', max.cutoff = 'q95') + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())+NoLegend()

## Mono
# anno_freq <- table(seu_4t1_mono_fin@meta.data$sample_rep, seu_4t1_mono_fin@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$group <- 'unknown'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D8_sham_rep',1:3),paste0('D14_sham_rep',1:3),paste0('D21_sham_rep',1:3),paste0('no_surg_rep',1:3)))] <- 'ctl'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D6_rep',1:3),paste0('D8_rep',1:3)))] <- 'early'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D10_rep',1:3),paste0('D12_rep',1:3)))] <- 'mid'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D14_rep',1:3),paste0('D21_rep',1:3)))] <- 'late'
# 
# anno_freq_summary_mono_inf_4t1 <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_mono_inf_4t1) <- c('celltype','group','mean','se')
# anno_freq_summary_mono_inf_4t1$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_mono_inf_4t1$group <- rep(unique(anno_freq$group), 2)
# for (i in 1:nrow(anno_freq_summary_mono_inf_4t1)) {
#   ct <- anno_freq_summary_mono_inf_4t1$celltype[i]
#   cl <- anno_freq_summary_mono_inf_4t1$group[i]
#   anno_freq_summary_mono_inf_4t1$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
#   anno_freq_summary_mono_inf_4t1$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
# }
# anno_freq_summary_mono_inf_4t1$group <- factor(anno_freq_summary_mono_inf_4t1$group, levels=c('ctl','early','mid','late'))
# anno_freq_summary_mono_inf_4t1$celltype <- factor(anno_freq_summary_mono_inf_4t1$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_mono_inf_4t1[which(anno_freq_summary_mono_inf_4t1$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics
# sce_mono <- as.SingleCellExperiment(seu_4t1_mono_fin)
# i.mono <- seq_len(ncol(sce_mono))
# boot.mono1 <- sample(i.mono, replace=TRUE)
# boot.mono2 <- sample(i.mono, replace=TRUE)
# sce_mono_rep2 <- sce_mono[,boot.mono1]
# sce_mono_rep3 <- sce_mono[,boot.mono2]
# 
# sample <- c(sce_mono$sample_rep,
#             paste0(sce_mono_rep2$sample_rep,'_boot2'),
#             paste0(sce_mono_rep3$sample_rep,'_boot3'))
# cluster <- c(sce_mono$inf,
#              sce_mono_rep2$inf,
#              sce_mono_rep3$inf)
# group <- c(sce_mono$group,
#            sce_mono_rep2$group,
#            sce_mono_rep3$group)
# allcounts <- cbind(counts(sce_mono),
#                    counts(sce_mono_rep2),
#                    counts(sce_mono_rep3))
# 
# sce_mono_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_mono_boot$sample <- sample
# sce_mono_boot$group <- group
# sce_mono_boot$cluster <- cluster
# ind_wt <- grep('ctl',sce_mono_boot$group)
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

## Baso
# anno_freq <- table(seu_4t1_baso_fin@meta.data$sample_rep, seu_4t1_baso_fin@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$group <- 'unknown'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D8_sham_rep',1:3),paste0('D14_sham_rep',1:3),paste0('D21_sham_rep',1:3),paste0('no_surg_rep',1:3)))] <- 'ctl'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D6_rep',1:3),paste0('D8_rep',1:3)))] <- 'early'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D10_rep',1:3),paste0('D12_rep',1:3)))] <- 'mid'
# anno_freq$group[which(anno_freq$Var1 %in% c(paste0('D14_rep',1:3),paste0('D21_rep',1:3)))] <- 'late'
# 
# anno_freq_summary_baso_inf_4t1 <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_baso_inf_4t1) <- c('celltype','group','mean','se')
# anno_freq_summary_baso_inf_4t1$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_baso_inf_4t1$group <- rep(unique(anno_freq$group), 2)
# for (i in 1:nrow(anno_freq_summary_baso_inf_4t1)) {
#   ct <- anno_freq_summary_baso_inf_4t1$celltype[i]
#   cl <- anno_freq_summary_baso_inf_4t1$group[i]
#   anno_freq_summary_baso_inf_4t1$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
#   anno_freq_summary_baso_inf_4t1$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$group == cl)])
# }
# anno_freq_summary_baso_inf_4t1$group <- factor(anno_freq_summary_baso_inf_4t1$group, levels=c('ctl','early','mid','late'))
# anno_freq_summary_baso_inf_4t1$celltype <- factor(anno_freq_summary_baso_inf_4t1$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_baso_inf_4t1[which(anno_freq_summary_baso_inf_4t1$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics
# sce_baso <- as.SingleCellExperiment(seu_4t1_baso_fin)
# i.baso <- seq_len(ncol(sce_baso))
# boot.baso1 <- sample(i.baso, replace=TRUE)
# boot.baso2 <- sample(i.baso, replace=TRUE)
# sce_baso_rep2 <- sce_baso[,boot.baso1]
# sce_baso_rep3 <- sce_baso[,boot.baso2]
# 
# sample <- c(sce_baso$sample_rep,
#             paste0(sce_baso_rep2$sample_rep,'_boot2'),
#             paste0(sce_baso_rep3$sample_rep,'_boot3'))
# cluster <- c(sce_baso$inf,
#              sce_baso_rep2$inf,
#              sce_baso_rep3$inf)
# group <- c(sce_baso$group,
#            sce_baso_rep2$group,
#            sce_baso_rep3$group)
# allcounts <- cbind(counts(sce_baso),
#                    counts(sce_baso_rep2),
#                    counts(sce_baso_rep3))
# 
# sce_baso_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_baso_boot$sample <- sample
# sce_baso_boot$group <- group
# sce_baso_boot$cluster <- cluster
# ind_wt <- grep('ctl',sce_baso_boot$group)
# ind_early <- grep('early',sce_baso_boot$group)
# ind_mid <- grep('mid',sce_baso_boot$group)
# ind_late <- grep('late',sce_baso_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_baso_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_baso_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_baso_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_baso_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_baso_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_baso_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_baso_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_baso_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_baso_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_baso_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_baso_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_baso_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_baso_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_baso_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_baso_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_baso_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_baso_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_baso_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_baso <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_baso) <- c('comparison','celltype','pval')
# prop_test_baso$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_baso$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_baso$pval <- pval_vec
# prop_test_baso$pval_bin <- 'ns'
# prop_test_baso$pval_bin[which(prop_test_baso$pval <= 0.05)] <- '*'

## Fig. 5H: Gonzalez et al MAM source and subtype UMAPs and MAM TLR-NFkB inflammatory signature feature and scatter plots
# GEX space
DimPlot(seu_mam_sub, group.by = 'source', cols=c('goldenrod','steelblue'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_mam_sub, group.by = 'celltype', cols=c('lightcoral','darkred'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())

# nmf19_genes_human <- toupper(nmf19_genes) 
# nmf19_genes_human <- nmf19_genes_human[which(nmf19_genes_human %in% rownames(seu_mam_sub@assays$RNA$counts))]
# seu_mam_sub <- AddModuleScore(seu_mam_sub, features = list(nmf19_genes_human[which(nmf19_genes_human %in% seu_mam_sub@assays$SCT@var.features)]), ctrl=5, name='nmf19')
FeaturePlot(seu_mam_sub, 'nmf19_mod1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

ggplot(seu_mam_sub@meta.data, aes(x=mam_s100a81, y=nmf19_mod1, color=celltype)) + geom_point(size=0.01) +
  theme_classic() + 
  geom_smooth(method=lm, color='black') +
  geom_point(data=seu_mam_sub@meta.data, aes(x=mam_s100a81, y=nmf19_mod1, color=celltype), size=0.01) +
  scale_color_manual(values=c('lightcoral','darkred')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
cor(seu_mam_sub@meta.data$mam_s100a81, seu_mam_sub@meta.data$nmf19_mod1)

ggplot(seu_mam_sub@meta.data, aes(x=mam_apoe1, y=nmf19_mod1, color=celltype)) + geom_point(size=0.01) +
  theme_classic() + 
  geom_smooth(method=lm, color='black') +
  geom_point(data=seu_mam_sub@meta.data, aes(x=mam_apoe1, y=nmf19_mod1, color=celltype), size=0.01) +
  scale_color_manual(values=c('lightcoral','darkred')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
cor(seu_mam_sub@meta.data$mam_apoe1, seu_mam_sub@meta.data$nmf19_mod1)
