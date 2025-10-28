scRNA4SeJeAP@meta.data$cellType2<- 'other'
Idents(scRNA4SeJeAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNA4SeJeAP,idents='6')
Idents(scRNA4SeJeAP)='cellType2'
Idents(scRNA4SeJeAP, cells = filtered_cellsAP6) <- 'AP6'
table(Idents(scRNA4SeJeAP))
scRNA4SeJeAP@meta.data$cellType2=Idents(scRNA4SeJeAP)

library(limma)
library(edgeR)
data=scRNA4SeJeAP@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))

group2<-  factor(scRNA4SeJeAP@meta.data$cellType2,levels = c("other", "AP6"))
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~group2)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='AP6'
sample2='other'

logFC <- lrt$table$logFC
lrt$table$SE <- abs(lrt$table$logFC) / sqrt(lrt$table$LR)
lrt$table$CI.L <- logFC - 1.96 * SE
lrt$table$CI.R <- logFC + 1.96 * SE
results <- data.frame(logFC, CI.L, CI.R, p.value = lrt$table$PValue)
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)
d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)
cut_off_FDR = 0.05  
cut_off_logFC = 0.5           
d1$change<- 'Stable'
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')

p <- ggplot(
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP6 v.s. other')
# p
ggsave(filename = "diffgene2nolabel.pdf", height = 6, width = 8, plot = p)


scRNA4SeJeAPFBOEC=subset(x = scRNA4SeJe, subset = (cellType == "AP"|cellType == "FBO"|cellType == "EC"))
library(dplyr)
 
scRNA4SeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAP@meta.data), var = "rowname")
scRNA4SeJeAP_tibble$rowname=rownames(scRNA4SeJeAP@meta.data)
scRNA4SeJeAPFBOEC_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAPFBOEC@meta.data), var = "rowname")
scRNA4SeJeAPFBOEC_tibble$rowname=rownames(scRNA4SeJeAPFBOEC@meta.data)
merged_df <- left_join(scRNA4SeJeAPFBOEC_tibble, scRNA4SeJeAP_tibble, by = "rowname")
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)

# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP6','other','FBO','EC'))
scRNA4SeJeAPFBOEC@meta.data$celltype2=merged_df$celltype2

dot=DotPlot(scRNA4SeJeAPFBOEC, features = c(genes),group.by='celltype2')
dotdata=dot$data
write.csv('scRNA4SeJeAPFBOEC_expression.csv',dotdata)

# filt expression result 
# 1) FBO%/AP*% > 2; 2) EC%/AP*% < 0.8 (to reduce the likelihood of endothelial transdifferentiation.); Secondary filters: 1) AP* expression ≥ 15% of cells; 2) FBO%/AP*% > 3.
library(readxl)
library(ggplot2)
library(ggrepel)
df <- read_excel("gene93.xlsx")
df <- df[,2:4]

p <- ggplot(df, aes(x = `pct.exp-AP6`, y = perFCFBOvsAP6)) +
  geom_point(
    aes(color = (`pct.exp-AP6` > 15 & perFCFBOvsAP6 > 3)), 
    size = 3, 
    alpha = 0.7) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red") 
  ) +guides(color = "none") +
  geom_vline(xintercept = 15,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 3,lty=4,col="black",lwd=0.8) +
  labs(
    x = "pct.exp-AP6 (%)", 
    y = "perFC(FBO vs AP6)",
    title = "") +
  scale_x_continuous(limits = c(0, 40))+
  scale_y_continuous(limits = c(0,6))+
  theme_bw() +theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  ) +geom_text_repel(
    data = df,
    aes(label = label),
    box.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20
  )
ggsave(filename = 'point93gene2.pdf',plot=p,width = 6,height = 6)


# beige
degs_bat <- FindMarkers(
  scRNA4SeJeAP,
  ident.1 = "6",  
  ident.2 = NULL,         
  min.pct = 0.25,
  logfc.threshold = 0.5
)

sig_genes <- rownames(subset(degs_bat, p_val_adj < 0.05 & avg_log2FC > 0.5))
scRNA4SiJi=subset(x = scRNA4, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNA4SiJiAP=subset(x = scRNA4SiJi, subset = (cellType == "AP"))
expr_matrix <- GetAssayData(scRNA4SiJiAP, slot = "data") 
gene_sets <- list(BAT_signature = sig_genes)  

cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  nCores = 50,      
  plotStats = TRUE
)

cells_AUC <- AUCell_calcAUC(
  gene_sets,
  cells_rankings,
  aucMaxRank = ceiling(nrow(cells_rankings)) 
)

auc_scores <- getAUC(cells_AUC)
rownames(auc_scores) <- names(gene_sets) 
scRNA4SiJiAP[["BAT_AUC"]] <- CreateAssayObject(data = auc_scores) 
scRNA4SiJiAP@assays$BAT_AUC@key <- "rna_"
DefaultAssay(scRNA4SiJiAP) <- "BAT_AUC"
scRNA4SiJiAP@meta.data$orig.ident=factor(scRNA4SiJiAP@meta.data$orig.ident,levels =c('Ji','Si'))
scRNA4SiJiAP$BAT_AUC_score <- auc_scores["BAT_signature", ]

high_auc_cells <- colnames(scRNA4SiJiAP)[scRNA4SiJiAP$BAT_AUC_score > quantile(scRNA4SiJiAP$BAT_AUC_score, 0.8)]
scRNA4SiJiAP$BAT_like <- ifelse(colnames(scRNA4SiJiAP) %in% high_auc_cells, "BAT-like", "Other")

DefaultAssay(scRNA4SiJiAP) <- "RNA"
# DefaultAssay(scRNA4SiJiAP) <- "integrated"
scRNA4SiJiAP <- NormalizeData(scRNA4SiJiAP, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA4SiJiAP <- FindVariableFeatures(scRNA4SiJiAP,nfeatures = 3000)
scRNA4SiJiAP <- ScaleData(scRNA4SiJiAP,vars.to.regress = 'percent.mt')
scRNA4SiJiAP <- RunPCA(scRNA4SiJiAP)
scRNA4SiJiAP <- FindNeighbors(scRNA4SiJiAP, reduction = "pca", dims = 1:30)
scRNA4SiJiAP <- FindClusters(scRNA4SiJiAP, resolution = 0.3)
scRNA4SiJiAP <- RunUMAP(scRNA4SiJiAP, reduction = "pca", dims = 1:30) 
library(harmony)
scRNA4SiJiAP <- RunHarmony(scRNA4SiJiAP, group.by.vars = "orig.ident")
scRNA4SiJiAP  <- FindNeighbors(scRNA4SiJiAP , dims = 1:30 , reduction = "harmony")
scRNA4SiJiAP  <- FindClusters(scRNA4SiJiAP, save.snn=T , resolution = 0.3)
scRNA4SiJiAP  <- RunUMAP(scRNA4SiJiAP , dims=1:30,reduction='harmony')
saveRDS(scRNA4SiJiAP,'scRNA4SiJiAP_harmony.rds')

DefaultAssay(scRNA4SiJiAP) <- "BAT_AUC"
pdf(file=sprintf("scRNA4SiJiAP_genesets.AUCdotplot.pdf"),height = 3, width = 5)
DotPlot(scRNA4SiJiAP,group.by = 'seurat_clusters', features =c("BAT-signature")) + coord_flip()
dev.off()

scRNA4SiJiAPFBOEC=subset(x = scRNA4SiJi, subset = (cellType == "AP"|cellType == "FBO"|cellType == "EC"))
scRNA4SiJiAP@meta.data$celltype2<- 'other'
Idents(scRNA4SiJiAP)='seurat_clusters'
filtered_cellsAP4 <- WhichCells(scRNA4SiJiAP,idents='4')
Idents(scRNA4SiJiAP)='celltype2'
Idents(scRNA4SiJiAP, cells = filtered_cellsAP4) <- 'AP4'
table(Idents(scRNA4SiJiAP))
scRNA4SiJiAP@meta.data$celltype2=Idents(scRNA4SiJiAP)

library(dplyr)
scRNA4SiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SiJiAP@meta.data), var = "rowname")
scRNA4SiJiAP_tibble$rowname=rownames(scRNA4SiJiAP@meta.data)
scRNA4SiJiAPFBOEC_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SiJiAPFBOEC@meta.data), var = "rowname")
scRNA4SiJiAPFBOEC_tibble$rowname=rownames(scRNA4SiJiAPFBOEC@meta.data)
merged_df <- left_join(scRNA4SiJiAPFBOEC_tibble, scRNA4SiJiAP_tibble, by = "rowname")
merged_df$celltype2 <- as.character(merged_df$celltype2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
merged_df$celltype_f <- ifelse(is.na(merged_df$celltype2), merged_df$cellType.x, merged_df$celltype2)
 
# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype_f<- factor(merged_df$celltype_f, levels = c('AP4','other','FBO','EC'))
scRNA4SiJiAPFBOEC@meta.data$celltype2=merged_df$celltype_f
table(scRNA4SiJiAPFBOEC@meta.data$celltype2)
saveRDS(scRNA4SiJiAPFBOEC,'scRNA4SiJiAPFBOEC.rds')

library(limma)
library(edgeR)
data=scRNA4SiJiAP@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)
group2<-  factor(scRNA4SiJiAP@meta.data$celltype2,levels = c("other", "AP4"))
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~group2)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='SiJiAP4'
sample2='SiJiother'
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)
cut_off_FDR = 0.05 
cut_off_logFC = 0.5 
d1$change<- 'Stable'
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
write.csv(file='sijiAP4diff.csv',d1)
p <- ggplot(
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP4 v.s. other')
# p
ggsave(filename = "sijidiffgene2nolabelAP4.pdf", height = 6, width = 6, plot = p)

genes=rownames(d1)[d1$change=='Up']
DefaultAssay(scRNA4SiJiAPFBOEC)='RNA'
dot=DotPlot(scRNA4SiJiAPFBOEC, features = c(genes),group.by='celltype2')
dotdata=dot$data
write.csv('scRNA4SiJiAPFBOEC_expression.csv',dotdata)

library(readxl)
library(ggplot2)
library(ggrepel)

data <- read_excel("gene49.xlsx")
data <- data[,2:4]

df <- data
p <- ggplot(df, aes(x = `pct.exp-AP4`, y = perFCFBOvsAP4)) +
  geom_point(
    aes(color = (`pct.exp-AP4` > 15 & perFCFBOvsAP4 > 3)),
    size = 3, 
    alpha = 0.7
  ) +scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  guides(color = "none") +
  geom_vline(xintercept = 15,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 3,lty=4,col="black",lwd=0.8) +
  labs(
    x = "pct.exp-AP4 (%)", 
    y = "perFC(FBO vs AP4)",
    title = ""
  ) + scale_x_continuous(limits = c(0, 40))+
  scale_y_continuous(limits = c(0,6))+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  ) +geom_text_repel(
    data = df,
    aes(label = label),
    box.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20
  )

ggsave(filename = 'point49gene.pdf',plot=p,width = 6,height = 6)

