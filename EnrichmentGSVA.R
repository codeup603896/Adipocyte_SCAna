# gsva分析
library(progeny)
library(readr)
library(pheatmap)
library(tibble)
library(GSVA)

# 棕色 myoFB
genelist<- list.files('./gmtps/')
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}

scRNA4SeJemyoFB=subset(x=scRNA4SeJe,subset=(cellType == "myoFB"))
data=scRNA4SeJemyoFB@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='scRNA4SeJemyoFBkegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNA4SeJemyoFB)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNA4SeJemyoFB)), 
                            CellType = as.character(Idents(scRNA4SeJemyoFB)),
                            stringsAsFactors = FALSE)
kegg_scores_df <- as.data.frame(t(kegg)) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell)
kegg_scores_df <- inner_join(kegg_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- kegg_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_kegg_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor, 
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SeJe-myoFB-gsva_hmap.pdf',progeny_hmap,height=6,width=10)

# 棕色 EC
genelist<- list.files('./gmtps/')

fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
}

scRNA4SeJeEC=subset(x=scRNA4SeJe,subset=(cellType == "EC"))
data=scRNA4SeJeEC@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='scRNA4SeJeECkegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNA4SeJeEC)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNA4SeJeEC)), 
                            CellType = as.character(Idents(scRNA4SeJeEC)),
                            stringsAsFactors = FALSE)
kegg_scores_df <- as.data.frame(t(kegg)) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell)
kegg_scores_df <- inner_join(kegg_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- kegg_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_kegg_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor,  
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SeJe-EC-gsva_hmap.pdf',progeny_hmap,height=6,width=10)

# 米色 myoFB
genelist<- list.files('./gmtps/')
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}

scRNA4SiJimyoFB=subset(x=scRNA4SiJi,subset=(cellType == "myoFB"))
data=scRNA4SiJimyoFB@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='scRNA4SiJimyoFBkegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNA4SiJimyoFB)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNA4SiJimyoFB)), 
                            CellType = as.character(Idents(scRNA4SiJimyoFB)),
                            stringsAsFactors = FALSE)
kegg_scores_df <- as.data.frame(t(kegg)) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell)
kegg_scores_df <- inner_join(kegg_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- kegg_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_kegg_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor, 
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SiJi-myoFB-gsva_hmap.pdf',progeny_hmap,height=6,width=10)

# # 米色 EC
genelist<- list.files('./gmtps/')

fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
}

scRNA4SiJiEC=subset(x=scRNA4SiJi,subset=(cellType == "EC"))
data=scRNA4SiJiEC@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='scRNA4SiJiECkegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNA4SiJiEC)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNA4SiJiEC)), 
                            CellType = as.character(Idents(scRNA4SiJiEC)),
                            stringsAsFactors = FALSE)
kegg_scores_df <- as.data.frame(t(kegg)) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell)
kegg_scores_df <- inner_join(kegg_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- kegg_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_kegg_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor,  
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SiJi-EC-gsva_hmap.pdf',progeny_hmap,height=6,width=10)

