# 细胞通讯
library(CellChat)
options(stringsAsFactors = FALSE)

# 读取数据
scRNA4@meta.data$samples=scRNA4@meta.data$orig.ident
# 拆分样本
Idents(scRNA4)='cellType'
scRNA4Se=subset(x = scRNA4, subset = (orig.ident == "Se"))
scRNA4Je=subset(x = scRNA4, subset = (orig.ident == "Je"))
scRNA4Si=subset(x = scRNA4, subset = (orig.ident == "Si"))
scRNA4Ji=subset(x = scRNA4, subset = (orig.ident == "Ji"))

i=1
samplenames=c('Se','Si','Je','Ji')
searchs=c('Secreted Signaling')
for (seuratobj_raw in c(scRNA4Se,scRNA4Si,scRNA4Je,scRNA4Ji)){
  seuratobj <- seuratobj_raw
  samplename=samplenames[i]

  data.input = seuratobj@assays$RNA@data # normalized data matrix
  meta = seuratobj@meta.data # a dataframe with rownames containing cell mata data
  ###创建CellChat 对象
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  # CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data 
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
  showDatabaseCategory(CellChatDB)

  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)

  for (search in searchs){
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # set the used database in the object
    cellchat@DB <- CellChatDB.use

    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multisession", workers = 50) # do parallel

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.mouse)
    # smoothData
    # cellchat <- smoothData(cellchat, adj = PPI.mouse)
    ###细胞通信网络的推断
    cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    groupSize <- as.numeric(table(cellchat@idents))
    save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))
    }

  i=i+1
}

searchsabbr=c('SS')
# 比较分组的样本名称，用于保存文件
comparegroup=c('SeJe','SiJi')
# 比较分组的样本在object.list所在的位置进行分组，用于后续分析
groups=list(c(1,3),c(2,4))

for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    # 整合Secreted Signaling cellchat
    load(sprintf('Se_cellchat_%s.RData',search))
    Se_cellchat=cellchat
    load(sprintf('Si_cellchat_%s.RData',search))
    Si_cellchat=cellchat
    load(sprintf('Je_cellchat_%s.RData',search))
    Je_cellchat=cellchat
    load(sprintf('Ji_cellchat_%s.RData',search))
    Ji_cellchat=cellchat

    object.list <- list(Je=Je_cellchat,Ji=Ji_cellchat,Se=Se_cellchat,Si =Si_cellchat)
    #run netAnalysis_computeCentrality
    object.list<- lapply(object.list,function(x){
        x=netAnalysis_computeCentrality(x)})
    cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

    ####从宏观角度预测细胞通讯####
    #比较交互总数和交互强度

    for (num in 1:length(groups)){
        compare=comparegroup[num]
        group=groups[[num]]
        # compare=comparegroup[3]
        # group=groups[[3]]
        pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        compareInteractions(cellchat_m, show.legend = F, group = (1:4),size.text = 20)#group颜色向量 默认measure='count'
        print(compareInteractions(cellchat_m, show.legend = F, group = (1:4),size.text = 20))
        compareInteractions(cellchat_m, show.legend = F, group = (1:4), measure = "weight",size.text = 20)
        print(compareInteractions(cellchat_m, show.legend = F, group = (1:4), measure = "weight",size.text = 20))
        dev.off()
        #不同细胞群之间的相互作用数量或强度的差异 circle
        pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group)#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
        dev.off() 

        #AP为中心不同细胞群之间的相互作用数量或强度的差异 circle
        pdf(file=sprintf("%sChatNumCellDIFF-AP_%s.pdf",compare,abbr),height = 6, width = 6)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group,sources.use='AP')#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置,compare[2]-compare[1]
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group,sources.use='AP')
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group,targets.use='AP')#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group,targets.use='AP')
        dev.off()

        #bubble        
        pdf(file=sprintf("%sChatbubble_%s-AME.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_bubble(cellchat_m, sources.use =c('AP','myoFB','EC'),targets.use = c('AP','myoFB','EC'),comparison = group, angle.x = 45, font.size = 16) +
        theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
        print(netVisual_bubble(cellchat_m, sources.use =c('AP','myoFB','EC'),targets.use = c('AP','myoFB','EC'),comparison =  group, angle.x = 45, font.size = 16)) +
        theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
        dev.off()
      }
}

cellchat_m.net <- subsetCommunication(cellchat_m)
write.xlsx(cellchat_m.net, "cellchat_m_net_lr.xlsx")
