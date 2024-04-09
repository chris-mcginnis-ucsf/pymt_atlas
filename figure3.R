#######################################################################
## The temporal progression of immune remodeling during metastasis  ###
## Companion code: Figure 3, Chris McGinnis, PhD, Stanford, 05/2023  ##
#######################################################################
## objects
load('seu_alv.Robj')
load('seu_alv_sub.Robj')
load("nmf_alv.Robj")
load('anno_freq_summary_alv_inf.Robj')
load('prop_test_alv_inf.Robj')
load('gsea_alv_inf.Robj')

#######################################################################################################################
## Figure 3: NMF and GSEA links Cd14+ inflammatory AM signature to TLR-NFκB inflammation and CD14+ ‘activated’ MDSCs ##
#######################################################################################################################
## Fig. 3A: NMF25 feature plots, AM subtype annotations on subsetted data
# A <- seu_alv@assays$RNA@counts[VariableFeatures(seu_alv), ]
# nmf_alv <- RcppML::nmf(A, 30, verbose = F, seed = 1234)
# seu_alv@meta.data[,paste0('nmf',1:30)] <- t(nmf_alv$h)
FeaturePlot(seu_alv, 'nmf8features1', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
DimPlot(seu_alv_sub, group.by = 'subtype_inf', cols=c('maroon','goldenrod','red','lightsalmon'), pt.size = 0.5) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_alv_sub, 'nmf8features1', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. 3B: TLR-NFκB inflammation signature feautre plots
FeaturePlot(seu_alv_sub, 'Cd14', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv_sub, 'Cxcl2', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv_sub, 'Nlrp3', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. 3C: Cd14+ inflammatory AM proportion bar charts
# anno_freq <- table(seu_alv_sub@meta.data$age_rep, seu_alv_sub@meta.data$inf)
# anno_freq <- anno_freq/rowSums(anno_freq)
# anno_freq <- melt(anno_freq)
# anno_freq$clade <- 'unknown'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_wt)] <- 'wt'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_early)] <- 'early'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_mid)] <- 'mid'
# anno_freq$clade[which(anno_freq$Var1 %in% clade_late)] <- 'late'
# 
# anno_freq_summary_alv_sub_inf <- as.data.frame(matrix(0L, nrow=2*4, ncol=4))
# colnames(anno_freq_summary_alv_sub_inf) <- c('celltype','clade','mean','se')
# anno_freq_summary_alv_sub_inf$celltype <- rep(unique(anno_freq$Var2), each=4)
# anno_freq_summary_alv_sub_inf$clade <- rep(unique(anno_freq$clade), 2)
# for (i in 1:nrow(anno_freq_summary_alv_sub_inf)) {
#   ct <- anno_freq_summary_alv_sub_inf$celltype[i]
#   cl <- anno_freq_summary_alv_sub_inf$clade[i]
#   anno_freq_summary_alv_sub_inf$mean[i] <- mean(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
#   anno_freq_summary_alv_sub_inf$se[i] <- stderror(anno_freq$value[which(anno_freq$Var2 == ct & anno_freq$clade == cl)])
# }
# anno_freq_summary_alv_sub_inf$clade <- factor(anno_freq_summary_alv_sub_inf$clade, levels=c('wt','early','mid','late'))
# anno_freq_summary_alv_sub_inf$celltype <- factor(anno_freq_summary_alv_sub_inf$celltype, levels=c('inf','no'))
ggplot(anno_freq_summary_alv_sub_inf[which(anno_freq_summary_alv_sub_inf$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

# Statistics 
# sce_alv_sub <- as.SingleCellExperiment(seu_alv_sub)
# i.alv_sub <- seq_len(ncol(sce_alv_sub))
# boot.alv_sub1 <- sample(i.alv_sub, replace=TRUE)
# boot.alv_sub2 <- sample(i.alv_sub, replace=TRUE)
# sce_alv_sub_rep2 <- sce_alv_sub[,boot.alv_sub1]
# sce_alv_sub_rep3 <- sce_alv_sub[,boot.alv_sub2]
# 
# sample <- c(sce_alv_sub$age_rep,
#             paste0(sce_alv_sub_rep2$age_rep,'_boot2'),
#             paste0(sce_alv_sub_rep3$age_rep,'_boot3'))
# cluster <- c(sce_alv_sub$inf_new,
#              sce_alv_sub_rep2$inf_new,
#              sce_alv_sub_rep3$inf_new)
# group <- c(sce_alv_sub$clade,
#            sce_alv_sub_rep2$clade,
#            sce_alv_sub_rep3$clade)
# allcounts <- cbind(counts(sce_alv_sub),
#                    counts(sce_alv_sub_rep2),
#                    counts(sce_alv_sub_rep3))
# 
# sce_alv_sub_boot <- SingleCellExperiment(assays = list(counts = allcounts))
# sce_alv_sub_boot$sample <- sample
# sce_alv_sub_boot$group <- group
# sce_alv_sub_boot$cluster <- cluster
# ind_wt <- grep('wt',sce_alv_sub_boot$group)
# ind_early <- grep('early',sce_alv_sub_boot$group)
# ind_mid <- grep('mid',sce_alv_sub_boot$group)
# ind_late <- grep('late',sce_alv_sub_boot$group)
# 
# propeller_list <- list()
# propeller_list[[1]] <- propeller(clusters = sce_alv_sub_boot$cluster[c(ind_wt,ind_early)],
#                                  sample = sce_alv_sub_boot$sample[c(ind_wt,ind_early)],
#                                  group = sce_alv_sub_boot$group[c(ind_wt,ind_early)], transform = 'asin')
# propeller_list[[2]] <- propeller(clusters = sce_alv_sub_boot$cluster[c(ind_wt,ind_mid)],
#                                  sample = sce_alv_sub_boot$sample[c(ind_wt,ind_mid)],
#                                  group = sce_alv_sub_boot$group[c(ind_wt,ind_mid)], transform = 'asin')
# propeller_list[[3]] <- propeller(clusters = sce_alv_sub_boot$cluster[c(ind_wt,ind_late)],
#                                  sample = sce_alv_sub_boot$sample[c(ind_wt,ind_late)],
#                                  group = sce_alv_sub_boot$group[c(ind_wt,ind_late)], transform = 'asin')
# propeller_list[[4]] <- propeller(clusters = sce_alv_sub_boot$cluster[c(ind_early,ind_mid)],
#                                  sample = sce_alv_sub_boot$sample[c(ind_early,ind_mid)],
#                                  group = sce_alv_sub_boot$group[c(ind_early,ind_mid)], transform = 'asin')
# propeller_list[[5]] <- propeller(clusters = sce_alv_sub_boot$cluster[c(ind_early,ind_late)],
#                                  sample = sce_alv_sub_boot$sample[c(ind_early,ind_late)],
#                                  group = sce_alv_sub_boot$group[c(ind_early,ind_late)], transform = 'asin')
# propeller_list[[6]] <- propeller(clusters = sce_alv_sub_boot$cluster[c(ind_mid,ind_late)],
#                                  sample = sce_alv_sub_boot$sample[c(ind_mid,ind_late)],
#                                  group = sce_alv_sub_boot$group[c(ind_mid,ind_late)], transform = 'asin')
# 
# prop_test_alv_sub_inf <- as.data.frame(matrix(0L, nrow=nrow(propeller_list[[1]])*6, ncol=3))
# colnames(prop_test_alv_sub_inf) <- c('comparison','celltype','pval')
# prop_test_alv_sub_inf$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), each = nrow(propeller_list[[1]]))
# prop_test_alv_sub_inf$celltype <- rep(rownames(propeller_list[[1]]), 6)
# pval_vec <- NULL
# for (i in 1:6) { pval_vec <- c(pval_vec, propeller_list[[i]][rownames(propeller_list[[1]]), 'P.Value']) }
# prop_test_alv_sub_inf$pval <- pval_vec
# prop_test_alv_sub_inf$pval_bin <- 'ns'
# prop_test_alv_sub_inf$pval_bin[which(prop_test_alv_inf$pval <= 0.05)] <- '*'

## Fig. 3D: NMF25 Hallmark GSEA lollipop plot 
# hall <- msigdb_gsets(species="Mus musculus", category="H")
# nmf8_features <- VariableFeatures(seu_alv)[which(nmf_alv$w[,8] >= 0.005)]
# gsea_hall_nmf8 <- hypeR(signature = nmf8_features, genesets = hall)
# gsea_alv_inf <- gsea_hall_nmf8$data[order(gsea_hall_nmf8$data$pval, decreasing = F), ]
ggplot(gsea_alv_inf[which(gsea_alv_inf$pval <= 0.01), ], aes(x=reorder(label, -log10(pval)), y = -log10(pval))) +
  geom_point(stat="identity",color='black') + 
  geom_segment(aes(x=label, xend=label, y=0, yend=-log10(pval))) + 
  geom_hline(yintercept = 2, linetype="dotted") + 
  coord_flip() + xlab('Pathway') + 
  theme_classic() + 
  theme(axis.text.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

## Fig. 3E: 'Activated' MDSC modeule feature plot
# act_genes <- c('Cd14','Ccl4', 'Ccl3', 'Cxcl2', 'Cxcl3', 'Spp1', 'Il1b', 'Nfkbia', 'Socs3', 'Mif', 'Klf6', 'Atf3', 'Ptgs2', 'Xbp1','Jun','Ccrl2','Saa3','Gad45b','Ninj1','Clec4n','Hcar2','Basp1','Btg1','Il1rn','Ifrd1','Txnip','Ccl9','Ier3','Ier5','Rs1','Thbs1','Cxcl1','Hilpda','Hist1h1c','Srgn','Hspa5','Csf1','Pgts2')
# seu_alv <- AddModuleScore(seu_alv, features = list(act_genes[which(act_genes %in% seu_alv@assays$SCT@var.features)]), ctrl=5, name='mdsc')
FeaturePlot(seu_alv, 'mdsc1', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. 3F: 'Activated' MDSC vs NMF25 module score scatter plot 
# seu_alv <- AddModuleScore(seu_alv, features = list(nmf8_features), ctrl=5, name='nmf8features')
# cor(seu_alv@meta.data$nmf25features1, seu_alv@meta.data$mdsc1)
ggplot(seu_alv@meta.data, aes(x=nmf8features1, y=mdsc1)) + 
  theme_classic() + 
  geom_smooth(method=lm, color='black') + 
  geom_point(data=seu_alv@meta.data, aes(x=nmf8features1, y=mdsc1, color=inf), size=0.01) + 
  scale_color_manual(values=c('red','black')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

