library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)


# matrix path
dir4 = c(
  "../Je_Ana/step3/filtered_feature_bc_matrix",
  "../Ji_Ana/step3/filtered_feature_bc_matrix",
  "../Se_Ana/step3/filtered_feature_bc_matrix",
  "../Si_Ana/step3/filtered_feature_bc_matrix"
)

scRNAlist4 <- list()  # creat a empty list
names(dir4) = c("Je", "Ji",'Se','Si')
sample_name = c("Je", "Ji",'Se','Si')

for (i in 1:length(dir4)){
  counts <- Read10X(data.dir = dir4[i])

  scRNAlist4[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist4[[i]] <- RenameCells(scRNAlist4[[i]],add.cell.id = sample_name[i])
}

# merge sample object
seurat_obj <- merge(scRNAlist4[[1]],y = c(scRNAlist4[[2]],scRNAlist4[[3]],scRNAlist4[[4]]))

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# QC cutoff
min_genes <- 200    
max_genes <- 4000    
max_mt <- 5        

qc_plots <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, group.by = "orig.ident") & 
  theme(plot.title = element_text(size=10))
# qc_plots
ggsave(filename = "qc_plots.pdf", height = 6, width = 6, plot = qc_plots)


pre_filter_counts <- table(seurat_obj$orig.ident)
cat("cell number before qc:\n")
print(pre_filter_counts)

seurat_obj_filtered <- subset(seurat_obj, 
                              subset = nFeature_RNA > min_genes & 
                                nFeature_RNA < max_genes & 
                                percent.mt < max_mt)

post_filter_counts <- table(seurat_obj_filtered$orig.ident)
cat("\ncell number after qc:\n")
print(post_filter_counts)

filter_stats <- data.frame(
  Sample = names(pre_filter_counts),
  Before_Filter = as.numeric(pre_filter_counts),
  After_Filter = as.numeric(post_filter_counts[match(names(pre_filter_counts), names(post_filter_counts))])
)

filter_stats$Percentage_Remaining <- round(filter_stats$After_Filter / filter_stats$Before_Filter * 100, 2)
filter_stats$Cells_Removed <- filter_stats$Before_Filter - filter_stats$After_Filter
filter_stats$Percentage_Removed <- 100 - filter_stats$Percentage_Remaining

cat("\nstats:\n")
print(filter_stats)

filter_stats_long <- gather(filter_stats, key = "Filter_Status", value = "Cell_Count", Before_Filter, After_Filter)

p_compare=ggplot(filter_stats_long, aes(x = Sample, y = Cell_Count, fill = Filter_Status)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Cell_Count), vjust = -0.3, position = position_dodge(0.9), size = 3) +
  labs(title = "Comparison", y = "cell num", x = "samples") +
  theme_bw() +theme(panel.border = element_blank(), axis.line = element_line(),panel.grid=element_blank()) +
scale_y_continuous(expand = c(0,0),limits = c(0,12000)) 
ggsave(filename = "p_compare.pdf", height = 4, width =6, plot =p_compare)

p_qc_filter=VlnPlot(seurat_obj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident") & 
  theme(plot.title = element_text(size=10))
ggsave(filename = "p_qc_filter.pdf", height = 6, width = 6, plot =p_qc_filter)

