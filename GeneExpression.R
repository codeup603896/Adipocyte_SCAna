DefaultAssay(scRNA4SeJe)='RNA'
gene='Nucb2'
# myoFB Nucb2 Ifg1r基因表达比例占比
scRNA4SeJemyoFB=subset(x=scRNA4SeJe,subset=(cellType == "myoFB"))
filtered_cells <- WhichCells(scRNA4SeJemyoFB, expression = Nucb2 > 0)
ra=length(filtered_cells)
filtered_cells <- WhichCells(scRNA4SeJemyoFB, expression = Igf1r > 0)
rb=length(filtered_cells)
filtered_cells <- WhichCells(scRNA4SeJemyoFB, expression = Igf1r > 0&Nucb2 > 0)
rc=length(filtered_cells)
filtered_cells <- WhichCells(scRNA4SeJemyoFB, expression = Igf1r == 0&Nucb2 == 0)
rd=length(filtered_cells)
a=ra-rc 
b=rb-rc 

df=data.frame(name=sprintf(c('myoFB','myoFB','myoFB','myoFB'),gene,gene,gene,gene),Type=sprintf(c('%s+Igf1r-','%s-Igf1r+','%s+Igf1r+','%s-Igf1r-'),gene,gene,gene,gene),Num=c(a,b,rc,rd))
df$Type=factor(df$Type,levels=c(sprintf(c('%s-Igf1r-','%s-Igf1r+','%s+Igf1r-','%s+Igf1r+'),gene,gene,gene,gene)))
myoFBNIdata = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(myoFBNIdata,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
ggsave(filename = sprintf("scRNA4SeJemyoFB_%s_Igf1r-cor.ratio.pdf",gene), plot = p, height = 6, width = 6)

# EC Nucb2 Vegfa基因表达比例占比
scRNA4SeJeEC=subset(x=scRNA4SeJe,subset=(cellType == "EC"))
filtered_cells <- WhichCells(scRNA4SeJeEC, expression =Nucb2> 0)
ra=length(filtered_cells)
filtered_cells <- WhichCells(scRNA4SeJeEC, expression = Vegfa > 0)
rb=length(filtered_cells)
filtered_cells <- WhichCells(scRNA4SeJeEC, expression = Vegfa > 0&Nucb2 > 0)
rc=length(filtered_cells)
filtered_cells <- WhichCells(scRNA4SeJeEC, expression = Vegfa == 0&Nucb2 == 0)
rd=length(filtered_cells)
a=ra-rc 
b=rb-rc 

df=data.frame(name=sprintf(c('EC','EC','EC','EC'),gene,gene,gene,gene),Type=sprintf(c('%s+Vegfa-','%s-Vegfa+','%s+Vegfa+','%s-Vegfa-'),gene,gene,gene,gene),Num=c(a,b,rc,rd))
df$Type=factor(df$Type,levels=c(sprintf(c('%s-Vegfa-','%s-Vegfa+','%s+Vegfa-','%s+Vegfa+'),gene,gene,gene,gene)))
ECNVdata = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(ECNVdata,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
ggsave(filename = sprintf("scRNA4SeJeEC_%s_Vegfa-cor.ratio.pdf",gene), plot = p, height = 6, width = 6)
