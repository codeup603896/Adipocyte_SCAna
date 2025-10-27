# 棕色脂肪组织通路活性分析
scRNA4SeJe=subset(x = scRNA4, subset = (orig.ident == "Se"|orig.ident == "Je"))

library(AUCell)
library(clusterProfiler)

fatty <- read.gmt("fat.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)
features=rownames(unique(fatty$term))
features=gsub( "_", "-",features)#后续要加入到seurat，行名不能含有下划线

sample1='Je'
sample2='Se'

cells_rankings <- AUCell_buildRankings(scRNA4SeJe@assays$RNA@data, plotStats=TRUE) 
# cells_rankings
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- getAUC(cells_AUC)
# 将aucs加入seurat
# aucs=t(aucs)
scRNA4SeJe_AUC <- CreateAssayObject(counts = aucs)
scRNA4SeJe@assays$AUC <- scRNA4SeJe_AUC
scRNA4SeJe@assays$AUC@key <- "rna_"
DefaultAssay(scRNA4SeJe) <- "AUC"
scRNA4SeJe@meta.data$orig.ident=factor(scRNA4SeJe@meta.data$orig.ident,levels =c(sample1,sample2))

# 作差异分析
library(limma)
library(edgeR)
scRNA4SeJeAUCdata=scRNA4SeJe@assays$AUC@data
scRNA4SeJedd_CM<- as.data.frame(as.matrix(scRNA4SeJeAUCdata))
# colnames(scRNA4SeJedd_CM)
# table(Idents(scRNA4SeJe))
scRNA4SeJegroup2<-  factor(c(rep('Je', length(rownames(scRNA4@meta.data[scRNA4@meta.data$orig.ident=='Je',]))),rep('Se', length(rownames(scRNA4@meta.data[scRNA4@meta.data$orig.ident=='Se',])))))#顺序靠后的Se vs 顺序靠前的Je
dgelist <- DGEList(counts =scRNA4SeJedd_CM, group = scRNA4SeJegroup2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~scRNA4SeJegroup2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
scRNA4SeJelrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(scRNA4SeJelrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

# 根据差异分析结果选择差异明显的通路进行作图
features=c("BIOCARTA-VOBESITY-PATHWAY","WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG","VERNOCHET-ADIPOGENESIS","GOBP-FAT-PAD-DEVELOPMENT","GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION")
Idents(scRNA4SeJe)='orig.ident'
pdf(file=sprintf("%s-ActivepathwayFC0.4.pdf",'SeJe'),height = 6, width = 12)
DotPlot(scRNA4SeJe,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
dev.off()
saveRDS(scRNA4SeJe,'scRNA4SeJe_AUC.rds')

# 米色脂肪组织通路活性分析
scRNA4SiJi=subset(x = scRNA4, subset = (orig.ident == "Si"|orig.ident == "Ji"))

fatty <- read.gmt("fat.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)
features=rownames(unique(fatty$term))
features=gsub( "_", "-",features)#后续要加入到seurat，行名不能含有下划线

sample1='Ji'
sample2='Si'

cells_rankings <- AUCell_buildRankings(scRNA4SiJi@assays$RNA@data, plotStats=TRUE) 
# cells_rankings
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- getAUC(cells_AUC)
# 将aucs加入seurat
# aucs=t(aucs)
scRNA4SiJi_AUC <- CreateAssayObject(counts = aucs)
scRNA4SiJi@assays$AUC <- scRNA4SiJi_AUC
scRNA4SiJi@assays$AUC@key <- "rna_"
DefaultAssay(scRNA4SiJi) <- "AUC"
scRNA4SiJi@meta.data$orig.ident=factor(scRNA4SiJi@meta.data$orig.ident,levels =c(sample1,sample2))

# 作差异分析

scRNA4SiJiAUCdata=scRNA4SiJi@assays$AUC@data
scRNA4SiJidd_CM<- as.data.frame(as.matrix(scRNA4SiJiAUCdata))
# colnames(scRNA4SiJidd_CM)
# table(Idents(scRNA4SiJi))
scRNA4SiJigroup2<-  factor(c(rep('Ji', length(rownames(scRNA4@meta.data[scRNA4@meta.data$orig.ident=='Ji',]))),rep('Si', length(rownames(scRNA4@meta.data[scRNA4@meta.data$orig.ident=='Si',])))))
dgelist <- DGEList(counts =scRNA4SiJidd_CM, group = scRNA4SiJigroup2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~scRNA4SiJigroup2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
scRNA4SiJilrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(scRNA4SiJilrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

# 根据差异分析结果选择差异明显的通路进行作图
features=c("BIOCARTA-VOBESITY-PATHWAY","WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG","VERNOCHET-ADIPOGENESIS","GOBP-FAT-PAD-DEVELOPMENT","GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION")
Idents(scRNA4SiJi)='orig.ident'
pdf(file=sprintf("%s-ActivepathwayFC0.4.pdf",'SiJi'),height = 6, width = 12)
DotPlot(scRNA4SiJi,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
dev.off()
saveRDS(scRNA4SiJi,'scRNA4SiJi_AUC.rds')
