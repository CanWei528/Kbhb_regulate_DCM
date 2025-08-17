###############################################################################
setwd('03.results')

# 七、X蛋白的β羟基丁酰化修饰在DCM中起了关键作用
if(!dir.exists('07.find_key_kbhb')){
  dir.create('07.find_key_kbhb')
}
color_samp = c('DCM'='#D96244','Control'='#5EB5C8')
################################################################################
# 1、基于差异β羟基丁酰化位点对应蛋白进行蛋白互作网络构建，筛选关键β羟基丁酰化蛋白X（图A）
de_info =read.table('05.kbhb_status/5_1_DE_kbhb_sites_results.txt',header=T,sep='\t')
DCM_up_genes= unique(de_info[de_info$regulate=='up','Gene.name'])
# DCM_down_genes= unique(de_info[de_info$regulate=='down','Gene.name'])

write(DCM_up_genes,file='07.find_key_kbhb/7_1_DCM_UP_genes.txt',ncolumns = 1)
# write(DCM_down_genes,file='07.find_key_kbhb/7_1_DCM_DOWN_genes.txt',ncolumns = 1)

################################################################################
key_genes = 'Atp5f1a'

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

################################################################################
# 2、条形图展示X蛋白在疾病vs对照组中β羟基丁酰化修饰水平（图B）
library(ggplot2)
data_7_2 = de_res[de_res$Gene.name==key_genes,]
data_7_2$sites2 = paste0(data_7_2$Amino.acid,data_7_2$Position)
data_7_2 = data_7_2[!is.na(data_7_2$regulate),]
data_7_2 = data_7_2[data_7_2$log2FC>0.2,]

pdf('07.find_key_kbhb/7_2_bar_key_kbhb_sites_ratio_of_Atp5f1a.pdf',width = 4,height = 3)
ggplot(data_7_2,aes(x=reorder(sites2,FC,decreasing=T),y=FC,fill=FC))+
  geom_bar(stat='identity',position = 'dodge2',width=0.5)+
  geom_hline(yintercept = 1,linetype='dotted',linewidth=1,color='grey')+
  theme_classic()+labs(x='Kbhb sites on Atp5f1a',y='DCM/Control Ratio')+
  theme(axis.title.x=element_text(size=14),
        axis.text = element_text(size=12),
        legend.position = 'none',
        plot.margin = unit(c(10,5,10,5),'mm'))+
  scale_y_continuous(expand=c(0,0),limits=c(0,2.5))+
  scale_fill_gradient(high=color_samp[['DCM']],low='#FFE4E1')
dev.off()
###############################################################################
# 3、热图展示X蛋白在疾病vs对照组中不同的β羟基丁酰化修饰位点β羟基丁酰化水平（图C）
# 热图
key_sites_de = de_res[de_res$Gene.name==key_genes,]
key_sites_de$sites2 = paste0(key_sites_de$Gene.name,'_',key_sites_de$Amino.acid,key_sites_de$Position)
key_sites_de = key_sites_de[!is.na(key_sites_de$regulate),]
key_sites_de = key_sites_de[key_sites_de$log2FC>0.2,]

data_7_3 = norm_exp[key_sites_de$kbhb_site,]
rownames(data_7_3) = key_sites_de$sites2

data_7_3_scale = apply(data_7_3,1,scale)
data_7_3_scale = t(data_7_3_scale)
colnames(data_7_3_scale) = colnames(data_7_3)

group2 = data.frame(
  row.names = colnames(data_7_3_scale),
  group=group
)

pdf('07.find_key_kbhb/7_3_heatmap.de_key_sites.pdf',height = 3,width = 5)
pheatmap::pheatmap(data_7_3_scale,scale='none',
                   annotation_col = group2,
                   cluster_cols = F,
                   annotation_colors=list(group=color_samp),
                   treeheight_row = 15,
                   cellwidth = 15,
                   cellheight = 15,
                   show_colnames = F,
                   color = colorRampPalette(c("navy","white", "firebrick3"))(50)
)
dev.off()
#---------------------------------------------------------------------------------
# complexheatmap 绘图
library(ComplexHeatmap)
library(circlize)
col_fun = colorRampPalette(c("#191970","white", "#8B0000"))(50)


show_sites= rownames(data_7_3_scale)
col_labels = rep('black',length(show_sites))
col_labels[3]='#D96244'
# at 是在原矩阵中的行索引。
# anno_meta = rowAnnotation(
#   foo = anno_mark(
#     at=match(show_sites,rownames(data_7_3_scale)),
#     labels=show_sites,
#     labels_gp = gpar(cex=1,col=col_labels)
#   )
# )# important_meta

df = data.frame(
  group=group,
  row.names = colnames(data_7_3_scale)
)

top_anno = HeatmapAnnotation(
  df=df,
  show_annotation_name = F,
  col=list(group=color_samp),
  annotation_legend_param = list(title=NULL)
)

ht = Heatmap(as.matrix(data_7_3_scale),
             cluster_rows = T,
             column_title = NULL,
             row_title = NULL,
             #show_parent_dend_line = F,
             row_names_max_width = unit(3,'cm'),
             column_dend_height = unit(0.5,'cm'), # 设置列聚类树的高度
             row_dend_width = unit(0.5,'cm'), # 设置行聚类树的宽度
             show_row_names = T,
             row_names_gp = gpar(cex=1,col=col_labels),
             #right_annotation = anno_meta,
             top_annotation = top_anno,
             show_column_names = F,
             heatmap_legend_param = list(title=NULL),
             heatmap_height = unit(40,'mm'),
             heatmap_width =unit(80,'mm') ,
             col=col_fun)


pdf('07.find_key_kbhb/7_3_heatmap_desites_ComplexHeatmap.pdf',width = 5,height = 4)
draw(ht)
dev.off()

################################################################################
# 4、箱线图展示X蛋白不同的β羟基丁酰化修饰位点β羟基丁酰化水平在疾病vs对照组中差异（图D）
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

fig_data_7_4 = cbind(kbhb=rownames(data_7_3_scale),data_7_3_scale)
fig_data_7_4 = as.data.frame(fig_data_7_4)
fig_data_7_4=melt(fig_data_7_4,value.name ='exp' ,id.vars = 'kbhb',variable.name = 'sample')
fig_data_7_4$group = gsub('_.*$','',fig_data_7_4$sample)
fig_data_7_4$group = factor(fig_data_7_4$group,levels=c('DCM','Control'))
fig_data_7_4$sites = gsub('^.*_','',fig_data_7_4$kbhb)
fig_data_7_4$exp=as.numeric(fig_data_7_4$exp)

pdf('07.find_key_kbhb/7_4_key_kbhb_compare.pdf',width = 5,height = 3.5)
ggboxplot(fig_data_7_4,x='sites',y='exp',color = 'group',add='jitter',width = 0.5)+
  scale_color_manual(values=color_samp)+
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.direction = 'horizontal')+
  stat_compare_means(aes(group=group),label.y=2.6,method='anova',label='p.signif')+
  ylim(c(-2,3))+
  labs(x='Kbhb sites on Atp5f1a',y='Relative Expression')
dev.off()
################################################################################
save.image('07.find_key_kbhb.0305.RData')
