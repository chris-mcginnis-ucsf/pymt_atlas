#######################################################################
## The temporal progression of immune remodeling during metastasis  ###
## Companion code: Figure 7, Chris McGinnis, PhD, Stanford, 05/2023  ##
#######################################################################
## objects
load('seu_imm_fin_2.Robj')
load('seu_int.Robj')
load('seu_alv.Robj')
load('seu_neu.Robj')
load('cellchat_list.Robj')
load('neu_ccl6_deg_table.Robj')
load('alv_ccl6_deg_table')
load('int_igf1_deg_table.Robj')

###########################################################################################################################
## Figure 7: Intercellular communication modeling reveals metastasis-associated changes in immune cell signaling network ##
###########################################################################################################################
## Fig. 7A: Cxcl2-Cxcr2 networks
# pathways.show <- 'CXCL'
# weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show) 
netVisual_individual(cellchat_list[[1]], edge.weight.max = weight.max[1], signaling = 'CXCL', pairLR.use = 'CXCL2_CXCR2', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[2]], edge.weight.max = weight.max[1], signaling = 'CXCL', pairLR.use = 'CXCL2_CXCR2', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[3]], edge.weight.max = weight.max[1], signaling = 'CXCL', pairLR.use = 'CXCL2_CXCR2', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)
netVisual_individual(cellchat_list[[4]], edge.weight.max = weight.max[1], signaling = 'CXCL', pairLR.use = 'CXCL2_CXCR2', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, graphics.init = F)

## Fig. 7B: Cell type x clade Igf1/Igf1r violin plots, IM Igf1 feature plots, and IM subset Igf1 violin plots
# Cell type x clade violin plots
# ident_vec <- paste(rep(c('neu','alv','mono','dc','int','t','b','nk'), each=4), rep(c('wt','early','mid','late'), 8), sep="_")
g <- VlnPlot(seu_imm_fin_2, idents = ident_vec, features = rev(c('Igf1','Igf1r')), stack = T, fill.by = 'ident', cols=rev(rep(c('palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4','darkred','indianred1','maroon'), each=4)), adjust = 1.5) + 
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), strip.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  scale_x_continuous(position = "top") + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=rev(ident_vec))
print(g)

# IM violin plots
VlnPlot(seu_int, sort=T, group.by = 'subtype', features = 'Igf1', cols=c('red','goldenrod','darkred','grey'), pt.size = 0)+NoLegend()+theme(plot.title = element_blank()) + 
  stat_summary(fun=mean, geom="point", size=1, color="black") + 
  theme(axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), plot.title = element_blank())

VlnPlot(seu_int, 'Igf1', idents = paste("mrc1",c('wt','early','mid','late'),sep='_'), sort=T, adjust = 1, pt.size = 0, cols=c('grey','lightcoral','red','maroon')) + NoLegend() + 
  stat_summary(fun=mean, geom="point", size=1, color="black") + 
  theme(axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm"), plot.title = element_blank())

# Statistics
# int_igf1_deg_table <- as.data.frame(matrix(1L, nrow = 6*1, ncol = 3))
# colnames(int_igf1_deg_table) <- c('subtype','comparison','pval')
# int_igf1_deg_table$subtype <- rep(c('mrc1'), each = 6)
# int_igf1_deg_table$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), 1)
# 
# for (i in 1:nrow(int_igf1_deg_table)) {
#   subtype <- int_igf1_deg_table$subtype[i]
#   comp <- unlist(strsplit(int_igf1_deg_table$comparison[i], split='_'))
#   ident1 <- paste(subtype,comp[1],sep='_')
#   ident2 <- paste(subtype,comp[2],sep='_')
#   temp <- FindMarkers(seu_int, ident.1 = ident1, ident.2 = ident2, features = 'Igf1')
#   if (nrow(temp) > 0) { int_igf1_deg_table$pval[i] <- as.numeric(temp$p_val) }
# }
# 
# int_igf1_deg_table$pval_bin <- 'ns'
# int_igf1_deg_table$pval_bin[which(int_igf1_deg_table$pval <= 0.05)] <- '*'

## Fig. 8C: Ccl6-Ccr1/Ccr2 networks, cell type x clade Ccl6 violin plot, and AM/Neutrophil subtype x clade Ccl6 violin plots
# Networks
# pathways.show <- 'CCL'
# weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show) 
netVisual_individual(cellchat_list[[3]], edge.weight.max = weight.max[1], signaling = 'CCL', pairLR.use = 'CCL6_CCR1', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)
netVisual_individual(cellchat_list[[3]], edge.weight.max = weight.max[1], signaling = 'CCL', pairLR.use = 'CCL6_CCR2', color.use=c('indianred1','darkred','maroon','palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)

# Cell type x clade violin plot
# ident_vec <- paste(rep(c('neu','alv','mono','dc','int','t','b','nk'), each=4), rep(c('wt','early','mid','late'), 8), sep="_")
g <- VlnPlot(seu_imm_fin_2, features = 'Ccl6', idents = ident_vec, cols=rep(c('palegreen3','skyblue3','dodgerblue','turquoise4','darkorchid4','darkred','indianred1','maroon'), each=4), adjust = 1.5, pt.size = 0) + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

# Neutrophil and AM subtype x clade violin plots
ident_vec <- paste(rep(c('mature','transition','immature','inf'), each=4), rep(c('wt','early','mid','late'), 4), sep="_")
g <- VlnPlot(seu_neu, features = 'Ccl6', idents = ident_vec, cols=rep(alpha('palegreen3',0.75), 16), adjust = 1.5, pt.size = 0) + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

ident_vec <- paste(rep(c('homeo','inf','lipid','allergy'), each=4), rep(c('wt','early','mid','late'), 4), sep="_")
g <- VlnPlot(seu_alv, features = 'Ccl6', idents = ident_vec, cols=rep(alpha('skyblue3',0.75), 16), adjust = 1.5, pt.size = 0) + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_line(size=1.25), axis.ticks = element_line(size=1.25), axis.ticks.length=unit(.2, "cm")) + 
  stat_summary(fun.y=mean, geom="point", size=1, color="black") 
g$data$ident <- factor(g$data$ident, levels=ident_vec)
print(g)

# Statistics
# neu_ccl6_deg_table <- as.data.frame(matrix(1L, nrow = 6*4, ncol = 3))
# colnames(neu_ccl6_deg_table) <- c('subtype','comparison','pval')
# neu_ccl6_deg_table$subtype <- rep(c('mature','transition','immature','inf'), each = 6)
# neu_ccl6_deg_table$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), 4)
# 
# for (i in 10:nrow(neu_ccl6_deg_table)) {
#   subtype <- neu_ccl6_deg_table$subtype[i]
#   comp <- unlist(strsplit(neu_ccl6_deg_table$comparison[i], split='_'))
#   ident1 <- paste(subtype,comp[1],sep='_')
#   ident2 <- paste(subtype,comp[2],sep='_')
#   temp <- FindMarkers(seu_neu, ident.1 = ident1, ident.2 = ident2, features = 'Ccl6')
#   if (nrow(temp) > 0) { neu_ccl6_deg_table$pval[i] <- as.numeric(temp$p_val) }
# }
# 
# neu_ccl6_deg_table$pval_bin <- 'ns'
# neu_ccl6_deg_table$pval_bin[which(neu_ccl6_deg_table$pval <= 0.05)] <- '*'
# neu_ccl6_deg_table$comparison <- factor(neu_ccl6_deg_table$comparison, levels = c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'))
# neu_ccl6_deg_table$subtype <- factor(neu_ccl6_deg_table$subtype, levels=c('mature','transition','immature','inf'))
# neu_ccl6_deg_table$pval_bin <- factor(neu_ccl6_deg_table$pval_bin, levels=c('ns','*'))
# neu_ccl6_deg_table <- acast(neu_ccl6_deg_table, comparison~subtype, value.var="pval_bin")
Heatmap(neu_ccl6_deg_table, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))

# alv_ccl6_deg_table <- as.data.frame(matrix(1L, nrow = 6*4, ncol = 3))
# colnames(alv_ccl6_deg_table) <- c('subtype','comparison','pval')
# alv_ccl6_deg_table$subtype <- rep(c('homeo','inf','lipid','allergy'), each = 6)
# alv_ccl6_deg_table$comparison <- rep(c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'), 4)
# 
# for (i in 1:nrow(alv_ccl6_deg_table)) {
#   subtype <- alv_ccl6_deg_table$subtype[i]
#   comp <- unlist(strsplit(alv_ccl6_deg_table$comparison[i], split='_'))
#   ident1 <- paste(subtype,comp[1],sep='_')
#   ident2 <- paste(subtype,comp[2],sep='_')
#   temp <- FindMarkers(seu_alv, ident.1 = ident1, ident.2 = ident2, features = 'Ccl6')
#   if (nrow(temp) > 0) { alv_ccl6_deg_table$pval[i] <- as.numeric(temp$p_val) }
# }
# 
# alv_ccl6_deg_table$pval_bin <- 'ns'
# alv_ccl6_deg_table$pval_bin[which(alv_ccl6_deg_table$pval <= 0.05)] <- '*'
# alv_ccl6_deg_table$comparison <- factor(alv_ccl6_deg_table$comparison, levels = c('wt_early','wt_mid','wt_late','early_mid','early_late','mid_late'))
# alv_ccl6_deg_table$celltype <- factor(alv_ccl6_deg_table$subtype, levels=c('homeo','inf','lipid','allergy'))
# alv_ccl6_deg_table$pval_bin <- factor(alv_ccl6_deg_table$pval_bin, levels=c('ns','*'))
# alv_ccl6_deg_table <- acast(alv_ccl6_deg_table, comparison~celltype, value.var="pval_bin")
Heatmap(alv_ccl6_deg_table, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))

