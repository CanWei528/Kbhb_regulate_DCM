############################################################
setwd("03.results")
############################################################
# 配色
# 分析图中db组命名为DCM，dbm组命名为Control
color_samp = c('DCM'='#D96244','Control'='#5EB5C8')
if(!dir.exists('05.kbhb_status')){
  dir.create('05.kbhb_status')
}
############################################################
# 读取数据
data_total = openxlsx::read.xlsx('01.data/kbhb报告/4-Differentially_expressed_protein/T-test_analysis/Differentially_annot.xlsx',sheet = 'dbvsdbm')
norm_exp = data_total[,c(1:3,8,17:22)]
Intensity_exp = data_total[,c(1:3,8,11:16)]
de_res = data_total[,c(1:5,8,23)]
rm(data_total)
#----------------------------------------------------
# 预处理数据
colnames(de_res)[4:5]=c('FC','pvalue')
de_res$kbhb_site = paste0(de_res$Protein.accession,'_',de_res$Amino.acid,de_res$Position,'_',de_res$Gene.name)
de_res$log2FC=log2(de_res$FC)
de_res$regulate = ifelse(
  de_res$log2FC>log2(1.2)&de_res$pvalue<0.05,'up',
  ifelse(de_res$log2FC< log2(5/6)&de_res$pvalue<0.05,'down','nosig')
)
de_res$FC=NULL
de_res$padjust = p.adjust(de_res$pvalue, method='fdr')

# 处理原始强度值
Intensity_exp$kbhb_site = paste0(Intensity_exp$Protein.accession,'_',Intensity_exp$Amino.acid,Intensity_exp$Position,'_',Intensity_exp$Gene.name)
Intensity_exp[,1:4]=NULL
rownames(Intensity_exp)=Intensity_exp$kbhb_site
Intensity_exp$kbhb_site=NULL
colnames(Intensity_exp) = gsub('Intensity.dbm','Control_',colnames(Intensity_exp))
colnames(Intensity_exp) = gsub('Intensity.db','DCM_',colnames(Intensity_exp))

# 处理相对表达值
norm_exp$kbhb_site = paste0(norm_exp$Protein.accession,'_',norm_exp$Amino.acid,norm_exp$Position,'_',norm_exp$Gene.name)
norm_exp[,1:4]=NULL
rownames(norm_exp)=norm_exp$kbhb_site
norm_exp$kbhb_site=NULL
colnames(norm_exp) = gsub('dbm','Control_',colnames(norm_exp))
colnames(norm_exp) = gsub('db','DCM_',colnames(norm_exp))

group=c(rep('Control',3),rep('DCM',3))
########################################################3333
# 合并表达数据和差异数据
# 相对表达结果  norm_exp
# 差异数据  de_res
# 原始强度数据  Intensity_exp
identical(rownames(Intensity_exp), rownames(norm_exp))
identical(de_res$kbhb_site, rownames(norm_exp))

norm_exp1 = norm_exp
Intensity_exp1 = Intensity_exp
colnames(norm_exp1) = paste0(colnames(norm_exp1),'.normalized')
colnames(Intensity_exp1) = paste0(colnames(Intensity_exp1),'.raw.intensity')
data_merged = cbind(de_res,norm_exp1,Intensity_exp1)
rm(norm_exp1,Intensity_exp1)
data_merged = data_merged[,c(1:3,5:7,4,10,8,9,11:22)]

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
sheet1=openxlsx::addWorksheet(wb, sheetName = 'Kbhb table and expressison')
openxlsx::writeDataTable(wb=wb,sheet=sheet1,x=data_merged,bandedRows = T)
nwid=nchar(colnames(data_merged))
nwid=nwid+3
nwid[7]=22
openxlsx::setColWidths(wb=wb,sheet=sheet1,cols=1:ncol(data_merged),widths=nwid)
openxlsx::freezePane(wb,sheet=sheet1, firstRow = TRUE)
openxlsx::saveWorkbook(wb, '05.kbhb_status/all_Kbhb_sites_table_and_expressison.xlsx', overwrite = TRUE)
#############################################################
#Atp5f1a修饰位点的FC和P值表格
need_gene = de_res[de_res$Gene.name=='Atp5f1a',]

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
sheet1=openxlsx::addWorksheet(wb, sheetName = 'Atp5f1a Differential Results')
openxlsx::writeDataTable(wb=wb,sheet=sheet1,x=need_gene,bandedRows = T)
nwid=nchar(colnames(need_gene))
nwid=nwid+5
nwid[7]=22
openxlsx::setColWidths(wb=wb,sheet=sheet1,cols=1:ncol(need_gene),widths=nwid)
openxlsx::freezePane(wb,sheet=sheet1, firstRow = TRUE)
openxlsx::saveWorkbook(wb, '05.kbhb_status/Atp5f1a_gene_differential_results.xlsx', overwrite = TRUE)

##################################################################
# 所有样本丁酰化蛋白水平和Kbhb位点丁酰化水平相关性
# 处理蛋白数据
data_path = "01.data/DIA报告/4-Differentially_expressed_protein/T-test_analysis/Differentially_annot.xlsx"
data_protein_DIA = openxlsx::read.xlsx(data_path,sheet = 'dbvsdbm') 
protein_de_res = data_protein_DIA[,c(1,3:5)]
colnames(protein_de_res)=c('protein','gene','FC','pvalue')
protein_de_res$log2FC = log2(protein_de_res$FC)
protein_de_res$regulate = ifelse(
  protein_de_res$log2FC>= log2(1.5) & protein_de_res$pvalue<0.05,'up',
  ifelse(protein_de_res$log2FC <= log2(2/3) & protein_de_res$pvalue<0.05,'down','nosig')
)
protein_de_res$ids = paste0(protein_de_res$protein,'__',protein_de_res$gene)
protein_de_res=protein_de_res[!is.na(protein_de_res$pvalue),]
rm(data_protein_DIA)
# 处理修饰数据
plot_data = de_res[,c(1:5,8,9)]
colnames(plot_data)[6]='logFC.kbhb'
plot_data$logFC.protein=protein_de_res[match(plot_data$Protein.accession,protein_de_res$protein),'log2FC']


fig__corr = cor.test(plot_data$logFC.kbhb,plot_data$logFC.protein,method='spearman')
fig__corr_pvalue=signif(fig__corr$p.value,2)
fig__corr_rho=signif(fig__corr$estimate[['rho']],3)
fig__label=paste0('R=',fig__corr_rho,', P=',fig__corr_pvalue)
plot_data$regulate=factor(plot_data$regulate,levels=c('down','nosig','up'))

pdf('05.kbhb_status/kbhb_protein_exp_corr.pdf',width = 6,height = 3.5)
ggplot()+
  geom_point(data=plot_data[plot_data$regulate=='nosig'|is.na(plot_data$regulate),],aes(x=logFC.kbhb,y=logFC.protein),color='grey')+
  geom_point(data=plot_data[plot_data$regulate!='nosig'&!is.na(plot_data$regulate),],aes(x=logFC.kbhb,y=logFC.protein,color=regulate),size=2)+
  geom_smooth(data=plot_data,aes(x=logFC.kbhb,y=logFC.protein),method='lm',color='red',se=F)+
  theme_test()+
  ylim(c(-4,5))+xlim(c(-4,6))+
  annotate('text',x=-2.3,y=4.5,label=fig__label,size=3,color='black')+
  theme(legend.title = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.7,0.9),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.text=element_text(size=8)
  )+
  scale_color_manual(values=c('down'=color_samp[['Control']],'nosig'='grey','up'=color_samp[['DCM']]))+
  guides(color=guide_legend(direction = 'horizontal'))+
  labs(x=bquote(log[2]~(DCM/Control~Ratio)~of~Kbhb),y=bquote(log[2]~(DCM/Control~Ratio)~of~Proteome))
dev.off()

############################################################
sig_kbhb_sites = de_res[!is.na(de_res$regulate)&de_res$regulate!='nosig',]
write.table(sig_kbhb_sites,file='05.kbhb_status/5_1_DE_kbhb_sites_results.txt',sep='\t',quote=F,row.names = F)
DCM_UP_sites = de_res[!is.na(de_res$regulate) & de_res$regulate=='up','kbhb_site']
############################################################
# 对DCM鉴定到的总体β羟基丁酰化修饰位点（kbhb）和蛋白数量进行描述（图A） 
library(ggplot2)

# 统计每个样本上检测到的位点个数， 仅是DCM组样本
kbhb_number_DCM = colSums(!is.na(norm_exp[,4:6]))
protein_num_DCM = c()# 从乳酸化数据中获得蛋白数量

for(i in 4:6){
  sites = rownames(norm_exp)[!is.na(norm_exp[,i])]
  sam = colnames(norm_exp)[i]
  protein_nm =length(unique(de_res[de_res$kbhb_site%in% sites, 'Protein.accession']))
  protein_num_DCM = c(protein_num_DCM,protein_nm)
  names(protein_num_DCM)[length(protein_num_DCM)]=sam
}
fig_data_5_1 = data.frame(
  sample=names(kbhb_number_DCM),
  Kbhb=as.numeric(kbhb_number_DCM),
  Protein=as.numeric(protein_num_DCM)
)
fig_data_5_1 = reshape2::melt(fig_data_5_1,id.vars='sample',value.name = 'count',variable.name='type')

pdf('05.kbhb_status/fig_5_1.pdf',width = 5,height = 3)
ggplot(fig_data_5_1,aes(x=sample,y=count,fill=type))+
  geom_bar(stat='identity',position = 'dodge2',width=0.8)+
  theme_test()+labs(y='No. of Identified Kbhb(Protein)')+
  geom_text(aes(label=count),vjust=-0.5,position = position_dodge2(0.9))+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size=10),
        legend.text = element_text(size=12))+
  scale_y_continuous(expand=c(0,0),limits=c(0,5800))
dev.off()
#-----------------------------------------------------------
# 统计每个样本上检测到的位点个数， 包括对照组
kbhb_number_DCM2 = colSums(!is.na(norm_exp[,1:6]))
protein_num_DCM2 = c()# 从乳酸化数据中获得蛋白数量

for(i in 1:6){
  sites = rownames(norm_exp)[!is.na(norm_exp[,i])]
  sam = colnames(norm_exp)[i]
  protein_nm =length(unique(de_res[de_res$kbhb_site%in% sites, 'Protein.accession']))
  protein_num_DCM2 = c(protein_num_DCM2,protein_nm)
  names(protein_num_DCM2)[length(protein_num_DCM2)]=sam
}
fig_data_5_1.1 = data.frame(
  sample=names(kbhb_number_DCM2),
  Kbhb=as.numeric(kbhb_number_DCM2),
  Protein=as.numeric(protein_num_DCM2)
)
fig_data_5_1.1 = reshape2::melt(fig_data_5_1.1,id.vars='sample',value.name = 'count',variable.name='type')

fig_data_5_1.1$sample = factor(fig_data_5_1.1$sample,levels=c(paste0('DCM_',1:3),paste0('Control_',1:3)))

pdf('05.kbhb_status/fig_5_1.full.pdf',width = 6,height = 2.5)
ggplot(fig_data_5_1.1,aes(x=sample,y=count,fill=type))+
  geom_bar(stat='identity',position = 'dodge2',width=0.8)+
  theme_test()+labs(y='No. of Identified Kbhb(Protein)')+
  geom_text(aes(label=count),vjust=-0.5,position = position_dodge2(0.9))+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size=10),
        legend.text = element_text(size=12))+
  scale_y_continuous(expand=c(0,0),limits=c(0,5850))
dev.off()


############################################################
# 样本每个蛋白值平均kbhb平均数量展示（图C）
fig_data_5_3 = de_res[,c(1,7)]
fig_data_5_3 = split(fig_data_5_3$kbhb_site,fig_data_5_3$Protein.accession)
kbhb_sites_num_per_pro=data.frame(
  protein = names(fig_data_5_3),
  kbhb_num=0,
  row.names = names(fig_data_5_3)
)

for(name in names(fig_data_5_3)){
  sites = fig_data_5_3[[name]]
  kbhb_sites_num_per_pro[kbhb_sites_num_per_pro$protein==name,'kbhb_num']=length(sites)
}
kbhb_sites_num_per_pro$type=kbhb_sites_num_per_pro$kbhb_num
kbhb_sites_num_per_pro$type=as.character(kbhb_sites_num_per_pro$type)
kbhb_sites_num_per_pro[kbhb_sites_num_per_pro$kbhb_num>15,'type']='16+'
#------------------------------------------
plot_data_5_3 = table(kbhb_sites_num_per_pro$type)
plot_data_5_3 = as.data.frame(plot_data_5_3)
colnames(plot_data_5_3)=c('type','num')
plot_data_5_3$type = factor(plot_data_5_3$type,levels=c(as.character(1:15),'16+'))


pdf('05.kbhb_status/fig_5_3_bar_kbhb_number_in_protein.pdf',width=5,height = 3)
ggplot(plot_data_5_3,aes(x=type,y=num))+
  geom_bar(stat='identity',width=0.8,fill='#D96244')+
  theme_classic()+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12)
  )+
  labs(y='Number of Proteins',x="No. of kbhb sites per protein")+
  scale_y_continuous(expand=c(0,0),limits = c(0,950))+
  geom_text(aes(label=num),vjust=-0.5,position = position_dodge2(0.9))
dev.off()

############################################################
# 疾病vs对照组中差异β羟基丁酰化修饰位点火山图展示（图D） 
fig_data_5_4 = de_res
fig_data_5_4 = fig_data_5_4[!is.na(fig_data_5_4$regulate),]
fig_data_5_4 = fig_data_5_4[!is.na(fig_data_5_4$pvalue),]

aa = table(fig_data_5_4$regulate)
names(aa)=c('DCM_Down','No_Sig','DCM_Up')
aa = as.data.frame(aa)
fig_data_5_4$regulate=factor(fig_data_5_4$regulate,levels=c('down','nosig','up'))

label_data = fig_data_5_4
label_data=label_data[label_data$regulate!='nosig',]
label_data=label_data[abs(label_data$log2FC)>=1.5 | -log10(label_data$pvalue)>=3,]
label_data= label_data[!(label_data$regulate=='down'&label_data$log2FC > -2.5),]


pdf("05.kbhb_status/fig_5_4_de_kbhb_sites_volcano.pdf",width = 7,height = 5)
ggplot(fig_data_5_4,aes(x=log2FC,y=-log10(pvalue)))+
  geom_point(size=2,aes(color=regulate))+
  labs(x=expression(log[2]~of~DCM/Control~ratio),y=expression(-log[10]~of~Pvalue))+
  geom_vline(xintercept = c(log2(5/6),log2(1.2)),lty=2,col ="#696969",lwd=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=2,col = "#696969",lwd=0.8)+ 
  theme_test()+theme(
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    legend.position = 'top',
    axis.text = element_text(size=12),
    axis.title=element_text(size=14),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(2,'mm'),
    legend.text = element_text(size=12)
  )+scale_color_manual(values = c('down'=color_samp[['Control']], 'nosig'='grey','up' =color_samp[['DCM']]),
                       breaks=c('down','nosig','up'),
                       labels=paste0(aa$Var1,':',aa$Freq))+
  geom_text_repel(data=label_data,aes(x=log2FC,y=-log10(pvalue),label=kbhb_site),show.legend = F)+
  geom_point(data=label_data,aes(x=log2FC,y=-log10(pvalue)),size=2.5,color='black',shape=21)
dev.off()

# fig_data_5_4[fig_data_5_4$Gene.name%in%c('Cpt1a','Hmgcs2','Acsm5'),]
# nosig

###########################################################
# 上调和下调差异β羟基丁酰化修饰位点、蛋白数量统计（图E）
fig_data_5_5 = fig_data_5_4[,c(1,9)]
fig_data_5_5=fig_data_5_5[fig_data_5_5$regulate!='nosig',]
fig_data_5_5$regulate=factor(fig_data_5_5$regulate,levels=c('up','down'))

protein_stat=as.data.frame.matrix(table(fig_data_5_5$Protein.accession,fig_data_5_5$regulate))
aa.pp=colSums(protein_stat!=0)
aa.kbhb=as.data.frame(table(fig_data_5_5$regulate))

plot_data_5_5=data.frame(
  number=c(99,149,90,129),
  regulated=c('Control_UP','DCM_UP','Control_UP','DCM_UP'),
  type=c('Kbhb','Kbhb','Protein','Protein')
)
plot_data_5_5$regulated=factor(plot_data_5_5$regulated,levels=c('DCM_UP','Control_UP'))

pdf('05.kbhb_status/fig_5_5_de_kbhb_protein.pdf',width=3,height = 3)
ggplot(plot_data_5_5,aes(x=type,y=number,fill=regulated))+
  geom_bar(stat='identity',position = 'dodge2',width=0.8)+
  theme_test()+labs(y='Number of Kbhb(Proteins)',fill=NULL)+
  geom_text(aes(label=number),vjust=-0.5,position = position_dodge2(0.9))+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = 'top'
  )+
  scale_y_continuous(expand = c(0,0),limits=c(0,180))+
  scale_fill_manual(values=c('Control_UP'=color_samp[['Control']], 'DCM_UP' =color_samp[['DCM']]))
dev.off()

##############################################################
# 上调和下调差异β羟基丁酰化修饰位点对应蛋白亚细胞定位展示
library(ggforce)
library(ggplot2)
library(stringr)

fig_data_5_6 = de_res[,c(1,4,6,7,9)]
fig_data_5_6$Regulated=ifelse(
  fig_data_5_6$regulate=='up','Up in DCM',
  ifelse(fig_data_5_6$regulate=='down','Down in DCM','nosig')
)
fig_data_5_6 = fig_data_5_6[!is.na(fig_data_5_6$regulate),]
fig_data_5_6 = fig_data_5_6[fig_data_5_6$regulate !='nosig',]
fig_data_5_6[!is.na(fig_data_5_6$Subcellular.localization)&fig_data_5_6$Subcellular.localization=='extracellular,plasma membrane','Subcellular.localization']='extracellular'
fig_data_5_6[!is.na(fig_data_5_6$Subcellular.localization)&fig_data_5_6$Subcellular.localization=='cytoplasm,nucleus','Subcellular.localization']='cytoplasm'

kbhb_site_stat = table(fig_data_5_6$Subcellular.localization,fig_data_5_6$Regulated)
kbhb_site_stat = as.data.frame.matrix(kbhb_site_stat)

protein_stat_from_kbhb=NULL
fig_data_5_6_1=fig_data_5_6[fig_data_5_6$Subcellular.localization!='',]
fig_data_5_6_1=fig_data_5_6_1[!is.na(fig_data_5_6_1$Subcellular.localization),]

for(loc_t in unique(fig_data_5_6_1$Subcellular.localization)){
  #loc_t = 'cytoplasm'
  tmp_Data = fig_data_5_6_1[fig_data_5_6_1$Subcellular.localization==loc_t,]
  up=length(unique(tmp_Data[tmp_Data$Regulated=='Up in DCM','Protein.accession']))
  down=length(unique(tmp_Data[tmp_Data$Regulated=='Down in DCM','Protein.accession']))
  protein_stat_from_kbhb = rbind(protein_stat_from_kbhb,
                                c(loc_type=loc_t,Up_in_DCM=up,Down_in_DCM=down)
  )
}
#---------------------------------------------
d=c(89,10,17,28,32,11)
d=d/sum(d)+0.6

u=c(87,4,18,84,50,22)
u=u/sum(u)+0.6

rr = c()
for(ii in 1:6){
  rr=c(rr,d[ii],u[ii])
}

fill_color=c('#5EB5C8', '#20A18D','#D96244',  '#E8A48D', '#59276E', '#47628E')
names(fill_color)=letters[1:6]

test=data.frame(
  x0=rep(seq(1,11,2),each=2),
  y0=rep(c(1,3),6),
  r=rr,
  fill=rep(letters[1:6],each=2)
)


text_size=3
pdf('05.kbhb_status/fig_5_6_circle_from_kbhb_data.pdf',width=5.5,height = 3)
ggplot(test)+geom_circle(aes(x0=x0,y0=y0,r=r,fill=fill),
                         alpha=0.8,color=NA)+
  scale_fill_manual(values=fill_color)+
  theme_classic()+
  scale_x_continuous(breaks=seq(1,11,2),labels = str_to_sentence(rownames(kbhb_site_stat)))+
  scale_y_continuous(breaks=c(1,3),labels=c('Down in DCM','Up in DCM'))+
  theme(legend.position = 'none',
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=45,hjust=1),
        axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth = 1),
        axis.title = element_blank())+
  annotate('text',x=1,y=1,label=expression(paste(frac(46,43))),size=text_size)+
  annotate('text',x=1,y=3,label=expression(paste(frac(44,43))),size=text_size)+
  annotate('text',x=3,y=1,label=expression(paste(frac(5,5))),size=text_size)+
  annotate('text',x=3,y=3,label=expression(paste(frac(2,2))),size=text_size)+
  annotate('text',x=5,y=1,label=expression(paste(frac(9,8))),size=text_size)+
  annotate('text',x=5,y=3,label=expression(paste(frac(9,9))),size=text_size)+
  annotate('text',x=7,y=1,label=expression(paste(frac(15,13))),size=text_size)+
  annotate('text',x=7,y=3,label=expression(paste(frac(45,39))),size=text_size)+
  annotate('text',x=9,y=1,label=expression(paste(frac(17,15))),size=text_size)+
  annotate('text',x=9,y=3,label=expression(paste(frac(26,24))),size=text_size)+
  annotate('text',x=11,y=1,label=expression(paste(frac(6,5))),size=text_size)+
  annotate('text',x=11,y=3,label=expression(paste(frac(11,11))),size=text_size)
dev.off()


############################################################
# 基于差异β羟基丁酰化修饰位点对应蛋白进行GO功能富集
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
if(!dir.exists('06.DE_kbhb_enrich')){
  dir.create('06.DE_kbhb_enrich')
}

DCM_UP_genes = unique(de_res[!is.na(de_res$regulate)&de_res$regulate=='up','Gene.name'])
DCM_UP_genes_ids=bitr(DCM_UP_genes,fromType = 'SYMBOL',toType ='ENTREZID',OrgDb = 'org.Mm.eg.db')

# GO 富集 up
ego_ALL_up <- enrichGO(gene = as.numeric(DCM_UP_genes_ids$ENTREZID), 
                       OrgDb = 'org.Mm.eg.db',
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2,
                       readable = TRUE) 
ego_ALL_up=as.data.frame(ego_ALL_up)
write.table(ego_ALL_up,file='06.DE_kbhb_enrich/6_1_enrich_GO_DCM_UP_sites.txt',sep='\t',quote=F,row.names = F)
# KEGG 富集 up
enrich_kegg_up <- enrichKEGG(gene = DCM_UP_genes_ids$ENTREZID,
                             organism = 'mmu',
                             pvalueCutoff = 0.05,
                             pAdjustMethod = 'BH',
)
enrich_kegg_up = DOSE::setReadable(enrich_kegg_up,OrgDb='org.Mm.eg.db',keyType = 'ENTREZID')
enrich_kegg_up = as.data.frame(enrich_kegg_up)
write.table(enrich_kegg_up,file='06.DE_kbhb_enrich/6_2_enrich_KEGG_DCM_UP_sites.txt',sep='\t',quote=F,row.names = F)

# GO可视化 up
library(dplyr)
library(stringr)
library(Hmisc)
plot_top_go_up = ego_ALL_up %>% 
  group_by(ONTOLOGY) %>% 
  slice_head(n=5)

go_fatty_relation = ego_ALL_up[grepl('fatty',ego_ALL_up$Description)&nchar(ego_ALL_up$Description)<35,]
plot_top_go_up = as.data.frame(plot_top_go_up)
plot_top_go_up = plot_top_go_up[plot_top_go_up$ONTOLOGY!='BP',]
plot_top_go_up = plot_top_go_up[plot_top_go_up$ONTOLOGY!='MF',]
plot_top_go_up = rbind(plot_top_go_up,go_fatty_relation)
plot_top_go_up = plot_top_go_up[!plot_top_go_up$ID%in% c('GO:1902001','GO:0035337','GO:0015909','GO:0015908'),]


plot_top_go_up$Description2 = capitalize(plot_top_go_up$Description)
plot_top_go_up=plot_top_go_up %>%
  arrange(ONTOLOGY,desc(Count))
plot_top_go_up$Description2=factor(plot_top_go_up$Description2,levels=rev(plot_top_go_up$Description2))

pdf('06.DE_kbhb_enrich/6_1_enrich_GO_DCM_UP_top5.pdf',width = 6,height = 3.5)
ggplot(plot_top_go_up,aes(x=Count,y=Description2,color=-log10(p.adjust),shape=ONTOLOGY))+
  geom_point(size=3)+xlim(c(min(plot_top_go_up$Count)-2,max(plot_top_go_up$Count)+2))+theme_bw()+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.5,'cm'),
        axis.text.y=element_text(size=12))+
  labs(shape='Ontology',color=bquote(-log[10]~(FDR)),x='Count')+
  scale_y_discrete(position = 'right',label=function(x){stringr::str_wrap(x,width = 60)} )+
  scale_color_gradient(low='blue',high='red')+
  guides(shape=guide_legend(order=1))
dev.off() 

# KEGG可视化 up
library(tidyverse)
library(cowplot)

plot_up_kegg = head(enrich_kegg_up,15)
plot_up_kegg = plot_up_kegg[plot_up_kegg$ID!='mmu05208',]
plot_up_kegg = rbind(plot_up_kegg,enrich_kegg_up[enrich_kegg_up$ID=='mmu00061',])

plot_up_kegg$Description2=gsub('\\s-\\sMus\\smusculus\\s\\(house\\smouse\\)','',plot_up_kegg$Description)

plot_up_kegg = plot_up_kegg[order(-log10(plot_up_kegg$pvalue),decreasing = F),]
plot_up_kegg$Description2 = factor(plot_up_kegg$Description2,levels=plot_up_kegg$Description2)
label_col = rep('black',15)
label_col[grepl('fatty',plot_up_kegg$Description2,ignore.case = T)]='red'

pdf('06.DE_kbhb_enrich/6_2_enrich_KEGG_DCM_UP_top15.pdf',width = 6,height = 3.5)
ggplot(plot_up_kegg,aes(x=-log10(pvalue),y=Description2,color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+
  xlim(sort(-log10(plot_up_kegg$pvalue))[c(1,15)]+c(-0.2,0.5))+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.4,'cm'),
        axis.text.y=element_text(size=12,color=label_col))+
  labs(x=bquote(-log[10]~(pvalue)))+
  scale_y_discrete(position = 'right')+
  scale_color_gradient(low='blue',high='red')+
  guides(color='none')
dev.off()
#------------------------------------------------------
DCM_DOWN_genes = unique(de_res[!is.na(de_res$regulate)&de_res$regulate=='down','Gene.name'])
DCM_DOWN_genes_ids=bitr(DCM_DOWN_genes,fromType = 'SYMBOL',toType ='ENTREZID',OrgDb = 'org.Mm.eg.db')


# GO 富集 down
ego_ALL_down <- enrichGO(gene = as.numeric(DCM_DOWN_genes_ids$ENTREZID), 
                         OrgDb = 'org.Mm.eg.db',
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.2,
                         readable = TRUE) 
ego_ALL_down=as.data.frame(ego_ALL_down)
write.table(ego_ALL_down,file='06.DE_kbhb_enrich/6_1_enrich_GO_DCM_DOWN_sites.txt',sep='\t',quote=F,row.names = F)
# KEGG 富集 down
enrich_kegg_down <- enrichKEGG(gene = as.numeric(DCM_DOWN_genes_ids$ENTREZID),
                               organism = 'mmu',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',use_internal_data = F)
enrich_kegg_down = DOSE::setReadable(enrich_kegg_down,OrgDb='org.Mm.eg.db',keyType = 'ENTREZID')
enrich_kegg_down = as.data.frame(enrich_kegg_down)
write.table(enrich_kegg_down,file='06.DE_kbhb_enrich/6_2_enrich_KEGG_DCM_DOWN_sites.txt',sep='\t',quote=F,row.names = F)
# GO可视化 down

plot_top_go_down = ego_ALL_down %>% 
  group_by(ONTOLOGY) %>% 
  slice_head(n=5)
plot_top_go_down = as.data.frame(plot_top_go_down)
plot_top_go_down$Description2 = capitalize(plot_top_go_down$Description)
plot_top_go_down=plot_top_go_down %>%
  arrange(ONTOLOGY,desc(Count))
plot_top_go_down$Description2=factor(plot_top_go_down$Description2,levels=rev(plot_top_go_down$Description2))

pdf('06.DE_kbhb_enrich/6_1_enrich_GO_DCM_DOWN_top5.pdf',width = 7,height = 3.5)
ggplot(plot_top_go_down,aes(x=Count,y=Description2,color=-log10(p.adjust),shape=ONTOLOGY))+
  geom_point(size=3)+xlim(c(min(plot_top_go_down$Count)-2,max(plot_top_go_down$Count)+2))+theme_bw()+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.5,'cm'),
        axis.text.y=element_text(size=12))+
  labs(shape='Ontology',color=bquote(-log[10]~(FDR)),x='Count')+
  scale_y_discrete(position = 'right',label=function(x){stringr::str_wrap(x,width = 60)} )+
  scale_color_gradient(low='blue',high='red')+guides(shape=guide_legend(order=1))
dev.off()

# KEGG可视化 down
plot_down_kegg = head(enrich_kegg_down,15)
plot_down_kegg$Description2=gsub('\\s-\\sMus\\smusculus\\s\\(house\\smouse\\)','',plot_down_kegg$Description)

plot_down_kegg = plot_down_kegg[order(-log10(plot_down_kegg$pvalue),decreasing = F),]
plot_down_kegg$Description2 = factor(plot_down_kegg$Description2,levels=plot_down_kegg$Description2)


pdf('06.DE_kbhb_enrich/6_2_enrich_KEGG_DCM_DOWN.pdf',width = 5,height = 2.5)
ggplot(plot_down_kegg,aes(x=-log10(pvalue),y=Description2,color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+
  xlim(sort(-log10(plot_down_kegg$pvalue))[c(1,10)]+c(-0.2,0.5))+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.4,'cm'),
        axis.text.y=element_text(size=12))+
  labs(x=bquote(-log[10]~(pvalue)))+
  scale_y_discrete(position = 'right')+
  scale_color_gradient(low='blue',high='red')+
  guides(color='none')
dev.off()

###############################################################################
setwd("03.results/")
save.image('05.kbhb_status_GO_KEGG.0302.RData')

###################################################################################