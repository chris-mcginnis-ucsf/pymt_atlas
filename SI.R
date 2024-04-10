###################################################################################
## The temporal progression of immune remodeling during metastasis  ###############
## Companion code: Supplemental Figures, Chris McGinnis, PhD, Stanford, 12/2023  ##
###################################################################################
############################################################################################################################################
## Figure S1. Summaries of MULTI-seq classification and results for cryopreserved PyMT, fresh PyMT, and 4T1 temporal cell atlas datasets. ##
############################################################################################################################################
## objects
load('cellcount_table_pymt_cryo.Robj')
load('cellcount_table_pymt_fresh.Robj')
load('cellcount_table_4t1.Robj')

## Fig. S1A: Sample classification results for cryopreserved PyMT lung atlas
ggplot(cellcount_table_pymt_cryo, aes(x=age_rep, y=count, fill=age)) + 
  geom_col(color='black',width=0.5) + 
  theme_classic() + 
  scale_fill_manual(values=c(viridis(9),'grey50')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1), axis.ticks = element_line(size=1), axis.ticks.length=unit(.15, "cm")) + 
  geom_hline(yintercept = 2915.3, lty=2)

## Fig. S1B: Sample classification results for fresh PyMT lung atlas
ggplot(cellcount_table_pymt_fresh, aes(x=sample, y=count, fill=stage)) + 
  geom_col(color='black',width=0.5) + 
  theme_classic() + 
  scale_fill_manual(values=c(viridis(3)[c(1,3,2)],'grey50')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1), axis.ticks = element_line(size=1), axis.ticks.length=unit(.15, "cm")) + 
  geom_hline(yintercept = 4259.444, lty=2)

## Fig. S1C: Sample classification results for cryopreserved 4T1 lung atlas
ggplot(cellcount_table_4t1, aes(x=sample_rep, y=count, fill=sample)) + 
  geom_col(color='black',width=0.5) + 
  theme_classic() + 
  scale_fill_manual(values=c(viridis(6),rep('grey50',4))) + 
  theme(legend.position = 'none', axis.line = element_line(size=1), axis.ticks = element_line(size=1), axis.ticks.length=unit(.15, "cm")) + 
  geom_hline(yintercept = 2065.5, lty=2)

################################################
## Figure S2. PyMT validation cohort analysis ##
################################################
## objects
load("seu_imm_fresh.Robj")
load('anno_freq_summary_imm_cryo_fresh.Robj')
load('prop_test_imm_fresh.Robj')
load('anno_freq_fresh_vs_cryo.Robj')
load('qpcr_summary_tit.Robj')
load('qpcr_summary_exp.Robj')

## Fig. S2C: PyMT validation cohort histology analysis
hist_summary <- as.data.frame(matrix(0L, nrow=10, ncol=3))
colnames(hist_summary) <- c('group','sample','num_mets')
hist_summary$group <- c(rep('early',3),rep('mid',4),rep('late',3))
hist_summary$sample <- c(1:3,1:4,1:3)
hist_summary$num_mets <- c(rep(0,7),10,15,7)
hist_summary$group <- factor(hist_summary$group, levels=c('early','mid','late'))

ggplot(hist_summary, aes(x=group, y=num_mets, fill=group)) + geom_boxplot() + geom_point() + theme_classic() + scale_fill_manual(values=alpha(c('lightcoral','red','maroon'),0.8)) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

## Fig. 2D: PyMT validation cohort qPCR results with tumor titration experiment
ggplot(qpcr_summary_tit, aes(x=sample, y=exp2_dct_mean)) +
  geom_col(color='black') + theme_classic() +
  geom_errorbar(aes(ymin=exp2_dct_mean-exp2_dct_se, ymax=exp2_dct_mean+exp2_dct_se), width=.2, position=position_dodge(.9)) +
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

ggplot(qpcr_summary_tit[which(qpcr_summary_tit$sample %in% c(paste0('titrate_500_',1:3), paste0('titrate_100_',1:3), paste0('titrate_50_',1:3), paste0('titrate_WT_',1:3))), ], aes(x=sample, y=exp2_dct_mean)) +
  geom_col(color='black') + theme_classic() + geom_hline(yintercept = 0.02, lty = 2, color='red', linewidth = 1) +
  geom_errorbar(aes(ymin=exp2_dct_mean-exp2_dct_se, ymax=exp2_dct_mean+exp2_dct_se), width=.2, position=position_dodge(.9)) +
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), 
        axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggplot(qpcr_summary_exp, aes(x=sample, y=exp2_dct_mean)) +
  geom_col(color='black') + theme_classic() +
  geom_errorbar(aes(ymin=exp2_dct_mean-exp2_dct_se, ymax=exp2_dct_mean+exp2_dct_se), width=.2, position=position_dodge(.9)) +
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

ggplot(qpcr_summary_exp[-grep('late', qpcr_summary_exp$sample), ], aes(x=sample, y=exp2_dct_mean)) +
  geom_col(color='black') + theme_classic() + geom_hline(yintercept = c(0.0938, 0.4017, 1.3696), lty = 2, color='red', linewidth = 1) + 
  geom_errorbar(aes(ymin=exp2_dct_mean-exp2_dct_se, ymax=exp2_dct_mean+exp2_dct_se), width=.2, position=position_dodge(.9)) +
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"),
        axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

## Fig. S2E: Immune cell type annotation UMAP and cell type proportion bar charts
DimPlot(seu_imm_fresh, group.by='celltype', cols=c('skyblue3','lightcoral','darkred','red','turquoise4','turquoise','dodgerblue','lightsalmon','darkorchid4','royalblue3','navy','palegreen3','maroon','turquoise3','grey','pink3'), raster=F, pt.size = 0.01) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

## Fig. S2F: Stage JSD heatmap 
Heatmap(jsd_imm_fresh_all, col = brewer.pal(10,'RdGy'), show_column_names = F, row_names_gp = gpar(fontsize=8))

## Fig. S2G: Cell type proportion scatter line plots
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'b'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'cd4t'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'cd8t'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'treg'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'ilc_gdt'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'nk'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'int'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'alv'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'neu'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'cm'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'intm'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'ncm'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'cdc1'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'cdc2'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'pdc'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'ccr7'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'alv'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_cryo_fresh[which(anno_freq_summary_imm_cryo_fresh$celltype == 'prolif'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('black','goldenrod')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())

## Fig. S2H: Fresh vs cryo cell type proportion comparison
ggplot(anno_freq_fresh_vs_cryo[which(anno_freq_fresh_vs_cryo$group == 'wt'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-4.1,4.1)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
ggplot(anno_freq_fresh_vs_cryo[which(anno_freq_fresh_vs_cryo$group == 'early'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-4.1,4.1)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
ggplot(anno_freq_fresh_vs_cryo[which(anno_freq_fresh_vs_cryo$group == 'mid'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-4.1,4.1)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
ggplot(anno_freq_fresh_vs_cryo[which(anno_freq_fresh_vs_cryo$group == 'late'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-4.1,4.1)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

################################################################################################
## Figure S3. Dot and feature plots for tissue-resident macrophage subtype annotation markers ##
################################################################################################
## objects
load('seu_int.Robj')
load('seu_alv.Robj')
load('int_markers.Robj')
load('alv_markers.Robj')

## Fig. S3A: Interstitial macrophage dot and feature plots
DimPlot(seu_int, group.by = 'subtype', cols = c('darkred','goldenrod','red'), pt.size = 1) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Mrc1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Cd74', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Crip1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Vim', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'Apoe', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_int, 'H2-Ab1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_int, features=unique(int_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(levels(g$data$id)))

## Fig. S3B: Alveolar macrophage dot and feature plots
DimPlot(seu_alv, group.by = 'subtype', cols = c('red','tan4','black','lightsalmon','maroon','pink2','goldenrod','grey'), pt.size = 0.1)+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Tnf', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Cxcl2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Cd63', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Fabp5', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Mt1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Fcer1g', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Sftpc', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Txnip', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Cd74', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Isg15', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_alv, features=unique(alv_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(levels(g$data$id)))

#############################################################################################
## Figure S4. Dot and feature plots for BM-derived myeloid cell subtype annotation markers ##
#############################################################################################
## objects
load('seu_mono.Robj')
load('seu_dc.Robj')
load('seu_neu.Robj')
load('mono_markers.Robj')
load('dc_markers.Robj')
load('neu_markers.Robj')

## Fig. S4A: Monocyte dot and feature plots
DimPlot(seu_mono, group.by = 'subtype', cols = c('dodgerblue','steelblue3','navy')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'Ccr2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'Fn1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'Cd74', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'H2-Ab1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'Ace', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_mono, 'Pou2f2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_mono, features=unique(mono_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(levels(g$data$id)))

## Fig. S4B: Neutrophil dot and feature plots
DimPlot(seu_neu, group.by = 'subtype', cols = c('dodgerblue','dodgerblue4','navy','black','aquamarine4','tan4','deepskyblue3','grey')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_neu, 'Nfkbia', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_neu, 'Ier2', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_neu, 'Retnlg', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_neu, 'Lcn2', max.cutoff = 'q95', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_neu, idents = unique(seu_neu@active.ident)[-8], features=unique(neu_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(c('mature','transition','immature','inf')))

## Fig. S4C: DC dot and feature plots
DimPlot(seu_dc, group.by = 'subtype', cols=c('tan4','deepskyblue3','steelblue3','navy','grey','black')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Clec9a', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Xcr1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Itgam', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Mgl2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Ifitm2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Bst2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Ccr7', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Mki67', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_dc, features=unique(dc_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(levels(g$data$id)))

####################################################################################################################################################
## Figure S5: Metastatic progression is associated with pathological activation and degranulation gene expression signature in mature neutrophils ##
####################################################################################################################################################
## objects
load("seu_fresh3_neu_mature.Robj")
load("seu_4t1_neu_mature_fin.Robj")
load('umap_neu_mature_fresh.Robj')
load('umap_neu_mature_4t1.Robj')
load('mature_neu_fresh_zscores.Robj')
load('mature_neu_4t1_zscores.Robj')

# Metastatic stage UMAPs
ggplot(umap_neu_mature_fresh, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_fresh, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_fresh, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_fresh[which(umap_neu_mature_fresh$clade == 'wt'), ], color = 'black') +
  theme_void() 
ggplot(umap_neu_mature_fresh, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_fresh, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_fresh, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_fresh[which(umap_neu_mature_fresh$clade == 'early'), ], color = 'lightcoral') +
  theme_void() 
ggplot(umap_neu_mature_fresh, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_fresh, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_fresh, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_fresh[which(umap_neu_mature_fresh$clade == 'mid'), ], color = 'red') +
  theme_void() 
ggplot(umap_neu_mature_fresh, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_fresh, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_fresh, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_fresh[which(umap_neu_mature_fresh$clade == 'late'), ], color = 'maroon') +
  theme_void() 

ggplot(umap_neu_mature_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_4t1, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_4t1[which(umap_neu_mature_4t1$clade == 'ctl'), ], color = 'black') +
  theme_void() 
ggplot(umap_neu_mature_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_4t1, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_4t1[which(umap_neu_mature_4t1$clade == 'early'), ], color = 'lightcoral') +
  theme_void() 
ggplot(umap_neu_mature_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_4t1, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_4t1[which(umap_neu_mature_4t1$clade == 'mid'), ], color = 'red') +
  theme_void() 
ggplot(umap_neu_mature_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_neu_mature_4t1, size=3.95, color="black") + 
  geom_point(data=umap_neu_mature_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_neu_mature_4t1[which(umap_neu_mature_4t1$clade == 'late'), ], color = 'maroon') +
  theme_void() 

# Feature plots
FeaturePlot(seu_fresh3_neu_mature, 'Ifitm1', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_fresh3_neu_mature, 'Tspo', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_fresh3_neu_mature, 'Igfbp6', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_fresh3_neu_mature, 'Lrg1', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

FeaturePlot(seu_4t1_neu_mature_fin, 'Ifitm1', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_4t1_neu_mature_fin, 'Tspo', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_4t1_neu_mature_fin, 'Igfbp6', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_4t1_neu_mature_fin, 'Lrg1', max.cutoff = 'q95', pt.size = 0.5) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

# Z-score heatmaps
Heatmap(mature_neu_fresh_zscores, cluster_rows = F, show_heatmap_legend = F)
Heatmap(mature_neu_4t1_zscores, cluster_rows = F, show_heatmap_legend = F)

#############################################################################################################
## Figure S6: Cd14+ inflammatory myeloid cells, excluding DCs, express Ccrl2 and are associated with NMF19 ##
#############################################################################################################
## objects
load("seu_neu.Robj")
load('seu_alv.Robj')
load('seu_dc.Robj')
load('seu_int.Robj')
load('seu_mono.Robj')

## Fig. S6A: Neutrophil and AM inf subset UMAPs and NMF19 feature plots
DimPlot(seu_neu, group.by = 'inf', cols=c('red','grey80'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_alv, group.by = 'inf_total', cols=c('red','grey80'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
FeaturePlot(seu_neu, 'nmf19_mod1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'nmf19_mod1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. S6B: DC NMF19 and TLR-NFkB marker feature plots 
FeaturePlot(seu_dc, 'nmf19_mod1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Cd14', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Cxcl2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_dc, 'Nlrp3', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Fig. S6C: Ccrl2 violin plots for (non-)inflammatory neutrophils, AMs, and monocytes
VlnPlot(seu_alv, 'Ccrl2', pt.size = 0, group.by = 'inf', cols=c('darkred','grey')) + theme(plot.title = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank()) + NoLegend() + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_mono, 'Ccrl2', pt.size = 0, group.by = 'inf', cols=c('darkred','grey')) + theme(plot.title = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank())+NoLegend() + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_neu, 'Ccrl2', pt.size = 0, group.by = 'inf', cols=c('darkred','grey')) + theme(plot.title = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank())+NoLegend() + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_int, 'Ccrl2', pt.size = 0, group.by = 'inf', cols=c('darkred','grey')) + theme(plot.title = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank())+NoLegend() + stat_summary(fun.y=mean, geom="point", size=1, color="black") 

################################################################################################################
## Fig. S8: Analysis of TLR-NFκB inflammation signature in PyMT validation cohort tissue-resident macrophages ##
################################################################################################################
## objects
load('seu_fresh_alv.Robj')
load('seu_fresh_int.Robj')

## Annotation UMAPs
DimPlot(seu_fresh_alv, group.by = 'subtype', cols=c('goldenrod','black','red','grey'), pt.size = 2)+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_fresh_int, group.by = 'subtype', cols=c('maroon','goldenrod','red','black'), pt.size = 2)+NoLegend()+NoAxes()+theme(plot.title = element_blank())

## Feature plots
FeaturePlot(seu_fresh_alv, 'nmf19_mod1', max.cutoff = 'q95', pt.size = 2) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_fresh_int, 'nmf19_mod1', max.cutoff = 'q95', pt.size = 2) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

## Inflammatory proportion bar charts 
ggplot(anno_freq_summary_alv_inf_fresh[which(anno_freq_summary_alv_inf_fresh$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

ggplot(anno_freq_summary_int_inf_fresh[which(anno_freq_summary_int_inf_fresh$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

##############################################
## Fig. S9: 4T1 longitudinal atlas analysis ##
##############################################
## objects
load('pt_size_4t1.Robj')
load('jsd_4t1_sample_all.Robj')
load('jsd_4t1_sample.Robj')
load('seu_4t1_imm_fin.Robj')
load("anno_freq_summary_imm_pymt_4t1.Robj")
load('prop_test_4t1.Robj')
load('anno_freq_summary_4t1_dc.Robj')
load('anno_freq_4t1_vs_pymt.Robj')
load('umap_baso_4t1.Robj')
load('baso_markers.Robj')
load('seu_4t1_baso_fin.Robj')
load('seu_4t1_alv_fin.Robj')
load('seu_4t1_int_fin.Robj')
load('seu_4t1_neu_mature_fin.Robj')

## Fig. S9A: 4T1 primary tumor growth dynamics
ggplot(pt_size_4t1, aes(x=time, y=vol)) + geom_point(aes(color=time), size=1.5) + geom_point(shape=1, size=1.5, color="black") + theme_classic() + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + scale_color_manual(values=viridis(8))

## Fig. S9B: 4T1 timepoint JSD heatmaps 
Heatmap(jsd_4t1_sample_all, col = brewer.pal(10,'RdGy'), show_column_names = F, row_names_gp = gpar(fontsize=8))
Heatmap(jsd_4t1_sample, col = brewer.pal(10,'RdGy'), show_column_names = F, row_names_gp = gpar(fontsize=8))

## Fig. S9D: 4T1 immune cell type annotation UMAP 
DimPlot(seu_4t1_imm_fin, group.by='celltype', cols=c('skyblue3','lightcoral','goldenrod','black','darkred','red','turquoise4','turquoise','dodgerblue','tan4','lightsalmon','darkorchid4','royalblue3','navy','palegreen3','maroon','turquoise3','grey','pink3'), raster=F, pt.size = 0.01) + NoAxes() + NoLegend() + theme(plot.title = element_blank())

## Fig. S9E: 4T1 cell type proportion scatter line plots
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'b'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'cd4t'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
gplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'cd8t'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'treg'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'ilc_gdt'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'nk'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'int'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'alv'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'neu'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'cm'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'intm'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'ncm'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'cdc1'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'cdc2'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'pdc'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'ccr7'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())
ggplot(anno_freq_summary_imm_pymt_4t1[which(anno_freq_summary_imm_pymt_4t1$celltype == 'prolif'), ], aes(x=group_vec, y=mean, color=dataset)) + geom_point(size=2) + geom_line(linewidth=1) + theme_classic() + scale_color_manual(values=c('tan4','black')) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), axis.text = element_blank(), axis.title = element_blank())

## Fig. S9F: 4T1 DC subtype proportion bar charts
ggplot(anno_freq_summary_dc_4t1, aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

## Fig. S9G: Fresh vs cryo cell type proportion comparison
ggplot(anno_freq_4t1_vs_pymt[which(anno_freq_4t1_vs_pymt$group == 'wt'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-5.3,5.3)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
ggplot(anno_freq_4t1_vs_pymt[which(anno_freq_4t1_vs_pymt$group == 'early'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-5.3,5.3)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
ggplot(anno_freq_4t1_vs_pymt[which(anno_freq_4t1_vs_pymt$group == 'mid'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-5.3,5.3)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
ggplot(anno_freq_4t1_vs_pymt[which(anno_freq_4t1_vs_pymt$group == 'late'), ], aes(x=log2fc, y=celltype, color=change)) + geom_point(size=2) + theme_classic() + geom_vline(xintercept = c(-1,0,1), lty=c(1,2,1)) + xlim(c(-5.3,5.3)) + scale_color_manual(values=c('steelblue3','black','red')) +   theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

## Fig. S9H: Basophil metastastic stage UMAPs and marker dot plots
ggplot(umap_baso_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_baso_4t1, size=3.95, color="black") + 
  geom_point(data=umap_baso_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_baso_4t1[which(umap_baso_4t1$clade == 'ctl'), ], color = 'black') +
  theme_void() 
ggplot(umap_baso_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_baso_4t1, size=3.95, color="black") + 
  geom_point(data=umap_baso_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_baso_4t1[which(umap_baso_4t1$clade == 'early'), ], color = 'lightcoral') +
  theme_void() 
ggplot(umap_baso_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_baso_4t1, size=3.95, color="black") + 
  geom_point(data=umap_baso_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_baso_4t1[which(umap_baso_4t1$clade == 'mid'), ], color = 'red') +
  theme_void() 
ggplot(umap_baso_4t1, aes(x=umap_1,y=umap_2)) + 
  geom_point(data=umap_baso_4t1, size=3.95, color="black") + 
  geom_point(data=umap_baso_4t1, size=3.55, color="gray90") +
  geom_point(data = umap_baso_4t1[which(umap_baso_4t1$clade == 'late'), ], color = 'maroon') +
  theme_void() 

g <- DotPlot(seu_4t1_baso_fin, features=unique(baso_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(c('ctl','early','mid','late')))

## Fig. S9I: 4T1 longitudinal atlas mature neutrophil, IM, and AM annotation UMAPs, inflammatory marker/module violin plots, and inflammatory subtype proportion bar charts
# UMAPs
DimPlot(seu_4t1_alv_fin, group.by = 'subtype', cols=c('maroon','black','red','lightsalmon','grey'), pt.size = 1)+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_4t1_int_fin, group.by = 'subtype', cols=c('darkred','goldenrod','red','black'), pt.size = 2)+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_4t1_neu_mature_fin, group.by = 'inf', cols=c('red','navy'), pt.size = 0.25)+NoLegend()+NoAxes()+theme(plot.title = element_blank())

# Violin plots
VlnPlot(seu_4t1_neu_mature_fin, group.by = 'inf', 'Cd14', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 4) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_alv_fin, group.by = 'inf', 'Cd14', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 4) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_int_fin, group.by = 'inf', 'Cd14', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 4) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 

VlnPlot(seu_4t1_neu_mature_fin, group.by = 'inf', 'Cxcl2', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 6) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_alv_fin, group.by = 'inf', 'Cxcl2', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 6) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_int_fin, group.by = 'inf', 'Cxcl2', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 6) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 

VlnPlot(seu_4t1_neu_mature_fin, group.by = 'inf', 'Nlrp3', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 3) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_alv_fin, group.by = 'inf', 'Nlrp3', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 3) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_int_fin, group.by = 'inf', 'Nlrp3', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 3) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 

VlnPlot(seu_4t1_neu_mature_fin, group.by = 'inf', 'nmf19_mod1', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 2) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_alv_fin, group.by = 'inf', 'nmf19_mod1', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 2) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 
VlnPlot(seu_4t1_int_fin, group.by = 'inf', 'nmf19_mod1', pt.size = 0, cols = c('darkred','grey'), sort=T, y.max = 2) + NoLegend() + theme(plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + stat_summary(fun.y=mean, geom="point", size=1, color="black") 

# Bar charts
ggplot(anno_freq_summary_neu_inf_4t1[which(anno_freq_summary_neu_inf_4t1$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

ggplot(anno_freq_summary_alv_inf_4t1[which(anno_freq_summary_alv_inf_4t1$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

ggplot(anno_freq_summary_int_inf_4t1[which(anno_freq_summary_int_inf_4t1$celltype == 'inf'), ], aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_fill_manual(values=c('grey50','lightcoral','red','maroon'))

################################################################################
## Figure S10. Dot and feature plots for lymphocyte subtype annotation markers ##
################################################################################
## objects
load('seu_nk.Robj')
load('seu_t.Robj')
load('seu_b.Robj')
load('seu_t_prolif.Robj')
load('nk_markers.Robj')
load('t_markers.Robj')
load('b_markers.Robj')
load('anno_freq_t_all_vs_prolif.Robj')

## Fig. S10A: NK cell dot and feature plots
DimPlot(seu_nk, group.by = 'subtype', cols=c('forestgreen','goldenrod','palegreen3','black'), pt.size = 1) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_nk, 'Ctla2a', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_nk, 'Prf1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_nk, 'Ccl4', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_nk, 'Mcm2', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_nk, features=unique(nk_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(c('mod','cyto','cyto_ccl','prolif')))

## Fig. S10B: T cell dot and feature plots
DimPlot(seu_t, group.by = 'subtype', cols = c('olivedrab3','darkseagreen4','navy','steelblue3','deepskyblue3','grey','black','forestgreen')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_t, 'Lef1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_t, 'Il2rb', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_t, 'Nkg7', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_t, 'Tmem176a', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_t, 'Itgb1', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_t, 'Ikzf2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_t, features=unique(t_markers), cols='RdGy', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(c('cd4_naive','cd8_naive','cd8_memory','cd8_effector','treg','cd4_cd29','gdT_ILC','prolif')))

## Fig. S10C: Proliferative T cell feature plots and frequency bar chart
DimPlot(seu_t_prolif, group.by = 'subtype2', cols = c('olivedrab3','navy','grey','darkseagreen4','forestgreen'), pt.size = 2) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

ggplot(anno_freq_t_all_vs_prolif, aes(x=subset, y=value, fill=subtype)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('olivedrab3','navy','grey','darkseagreen4','forestgreen')) + 
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))

## Fig. S10D: B cell dot and feature plots
DimPlot(seu_b, group.by = 'subtype', cols=c('forestgreen','olivedrab3','palegreen3','darkseagreen4','grey'), pt.size = 0.1) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_b, 'Plac8', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_b, 'Ccnd2', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_b, 'Ighd', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_b, 'Ighm', max.cutoff = 'q95') + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

g <- DotPlot(seu_b, features=unique(b_markers), cols='RdGy', dot.scale = 3, idents = c('ighm','naive','breg','crip1')) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"))
g$data$id <- factor(g$data$id, levels=rev(c('naive','ighm','breg')))

####################################################################################################
## Fig. S11: Flow cytometry analysis of natural killer cell subtypes during metastatic progression ##
####################################################################################################
## objects
load('nk_facs_gating.Robj')

## Fig. S11B: Quadrant CD27xCD11b gating strategy and validation using markers of immunomodulation (THY1) and cytotoxicity (KLRG1)
ggplot(nk_facs_gating, aes(x = log10(thy1), color = gate)) + 
  geom_density(alpha=0, size=1.25) + 
  theme_classic() + 
  scale_color_manual(values=rev(c('cadetblue2','dodgerblue','navy'))) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) 

ggplot(nk_facs_gating, aes(x = log10(klrg1), color = gate)) + 
  geom_density(alpha=0, size=1.25) + 
  theme_classic() + 
  scale_color_manual(values=c('cadetblue2','dodgerblue','navy')) + 
  theme(legend.position = 'none', axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  xlim(c(2.6,3.1))

#####################################################################################################################
## Fig. S12: Interrogating cell-cell communication predictions associated with the TLR-NFκB inflammatory signature ##
#####################################################################################################################
## objects
load('seu_imm_fin_2.Robj')
load("cellchat_list.Robj")
load('seu_neu.Robj')
load('seu_alv.Robj')
load('seu_fresh_neu.Robj')
load('seu_fresh_alv.Robj')
load('seu_4t1_alv_fin.Robj')
load('seu_4t1_neu_fin.Robj')
load('seu_4t1_baso_fin.Robj')
load('alv_fresh_ccl6_deg_table.Robj')
load('neu_fresh_ccl6_deg_table.Robj')
load('alv_neu_ccl6_deg_table.Robj')
load('neu_neu_ccl6_deg_table.Robj')

## Fig. S12A: Tnf-Tnfrsf1a/Tnfrsf1b networks
pathways.show <- 'TNF'
weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
netVisual_individual(cellchat_list[[1]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1A', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[2]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1A', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[3]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1A', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[4]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1A', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)

netVisual_individual(cellchat_list[[1]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1B', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[2]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1B', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[3]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1B', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[4]], edge.weight.max = weight.max[1], signaling = 'TNF', pairLR.use = 'TNF_TNFRSF1B', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)

## Fig. S12B: Cell type x clade IL1b/Il1r1 violin plots, AM/neutrophil Il1rn feature and violin plots
# Cell type x clade violin plots
# ident_vec <- paste(rep(c('neu','alv','mono','dc','int','t','b','nk'), each=4), rep(c('wt','early','mid','late'), 8), sep="_")
g <- VlnPlot(seu_imm_fin_2, idents = ident_vec, features = rev(c('Il1b','Il1r1')), stack = T, fill.by = 'ident', cols=rev(rep(c('palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4','darkred','indianred1','maroon'), each=4)), adjust = 1.5) + 
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), strip.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_x_continuous(position = "top") + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=rev(ident_vec))
print(g)

# AM/neutrophil feature and violin plots
FeaturePlot(seu_neu, 'Il1rn', max.cutoff = 'q90', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())
FeaturePlot(seu_alv, 'Il1rn', max.cutoff = 'q90', pt.size = 0.1) + NoLegend() + NoAxes() + scale_color_viridis() + theme(plot.title = element_blank())

VlnPlot(seu_neu, sort=T, features = 'Il1rn', group.by = 'inf', cols=c('darkred','grey'), adjust = 2, pt.size = 0) + 
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), plot.title = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 

VlnPlot(seu_alv, sort=T, features = 'Il1rn', group.by = 'inf', cols=c('darkred','grey'), adjust = 2, pt.size = 0) + 
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), plot.title = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 

## Fig. S12C: PyMT fresh AM/Neutrophil subtype x clade Ccl6 violin plots
ident_vec <- paste(rep(c('homeo','inf','antigen'), each=4), rep(c('wt','early','mid','late'), 3), sep="_")
g <- VlnPlot(seu_fresh_alv, features = 'Ccl6', idents = ident_vec, cols=rep(alpha('skyblue3',0.75), 16), adjust = 1.5, pt.size = 0, assay = 'RNA') + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

ident_vec <- paste(rep(c('mature','transition','immature','inf'), each=4), rep(c('wt','early','mid','late'), 4), sep="_")
g <- VlnPlot(seu_fresh_neu, features = 'Ccl6', idents = ident_vec, cols=rep(alpha('palegreen3',0.75), 16), adjust = 1.5, pt.size = 0, assay = 'RNA') + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

Heatmap(neu_fresh_ccl6_deg_table, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))
Heatmap(alv_fresh_ccl6_deg_table, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))

## Fig. S12D: 4T1 AM/Neutrophil/Basophil subtype x clade Ccl6 violin plots
ident_vec <- paste(rep(c('homeo','inf','allergy','lipid'), each=4), rep(c('ctl','early','mid','late'), 3), sep="_")
g <- VlnPlot(seu_4t1_alv_fin, features = 'Ccl6', idents = ident_vec, cols=rep(alpha('skyblue3',0.75), 16), adjust = 1.5, pt.size = 0, assay = 'RNA') + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

ident_vec <- paste(rep(c('mature','transitional','immature','inf'), each=4), rep(c('ctl','early','mid','late'), 4), sep="_")
g <- VlnPlot(seu_4t1_neu_fin, features = 'Ccl6', idents = ident_vec, cols=rep(alpha('palegreen3',0.75), 16), adjust = 1.5, pt.size = 0, assay = 'RNA') + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

VlnPlot(seu_4t1_baso_fin, features = 'Ccl6', group.by = 'group', cols=c('goldenrod','goldenrod','goldenrod','goldenrod'), adjust = 1.5, pt.size = 0, assay = 'RNA') + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 

Heatmap(neu_4t1_ccl6_deg_table, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))
Heatmap(alv_4t1_ccl6_deg_table, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))

