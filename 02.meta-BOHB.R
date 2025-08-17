################################################################################
setwd('03.results')
# db/ dbm
# db = DCM
# dbm = Control
if(!dir.exists('03.meta.DCM')){
  dir.create('03.meta.DCM')
}
color_samp = c('DCM'='#D96244','Control'='#5EB5C8')

library(ggplot2)
############################################################################
# 三、基于代谢组数据分析糖尿病心肌病中的β羟基丁酰化相关代谢异常
# 使用数据
# 01.data\代谢-小鼠(MJ20240401212) - 10 - 1\interaction_results\4.ExpDiff\01.TwoGroupExpDiff\ExpDiff_s1s2s3s4s5_d1d2d3d4d5\mix\DiffTest\db_vs_dbm.diff.exp.xls
data_2 = read.table('01.data/代谢-小鼠(MJ20240401212) - 10 - 1/interaction_results/4.ExpDiff/01.TwoGroupExpDiff/ExpDiff_s1s2s3s4s5_d1d2d3d4d5/mix/DiffTest/db_vs_dbm.diff.exp.xls',header=T,sep='\t',quote='"',comment.char = '!')

de_res = data_2[,1:6]
de_res = de_res[!grepl('^metab_',de_res$Metabolite),]
de_res$logFC = log2(de_res$FC)
de_res$regulate = ifelse(
  de_res$Vip_oplsda>=1&de_res$P_value<0.05,'sig','nosig'
)
de_res$regulate2 = ifelse(
  de_res$logFC>=0.2&de_res$P_value<0.05,'DCM_Up',
  ifelse(de_res$logFC<= -0.2&de_res$P_value<0.05,'DCM_Down','NoSig')
)

exp_data = data_2[,c(1,7:16)]
exp_data = exp_data[!grepl('^metab_',exp_data$Metabolite),]
rownames(exp_data)=exp_data$Metabolite
exp_data$Metabolite=NULL
colnames(exp_data) = gsub('dbm_s','Control_',colnames(exp_data))
colnames(exp_data) = gsub('db_d','DCM_',colnames(exp_data))
group=c(rep('Control',5),rep('DCM',5))
################################################################################
# 数据处理
# 图 AB 数据
#\代谢-小鼠(MJ20240401212) - 10 - 1\interaction_results\2.SampleComp\04.ExpPLSDA\ExpPLSDA_20240621_140650
fig_data_3_1 = read.table('01.data/代谢-小鼠(MJ20240401212) - 10 - 1/interaction_results/2.SampleComp/04.ExpPLSDA/ExpPLSDA_20240621_140650/PLS-DA.sites.xls',header = T,sep='\t',row.names = 1)
fig_data_3_1 = fig_data_3_1[,-c(3,4)]
fig_data_3_1 = cbind(sample=rownames(fig_data_3_1),fig_data_3_1)
fig_data_3_1$sample = gsub('dbm_s','Control_',fig_data_3_1$sample)
fig_data_3_1$sample = gsub('db_d','DCM_',fig_data_3_1$sample)
fig_data_3_1$group=gsub('_.*$','',fig_data_3_1$sample)
fig_data_3_1$group[1:3]='QC'

fig_data_3_1$group = factor(fig_data_3_1$group,levels=c('DCM','Control','QC'))

pdf('03.meta.DCM/3_1_A_PLSDA_PCA.pdf',width = 6,height = 5)
ggplot(fig_data_3_1,aes(x=p1,y=p2,color=group))+
  geom_hline(yintercept = 0,color='grey',linetype=81,linewidth=1)+
  geom_vline(xintercept = 0,color='grey',linetype=81,linewidth=1)+
  geom_point(size=3)+
  geom_point(size=3.5,shape=21,color='black')+
  theme_test()+
  scale_y_continuous(breaks=seq(-60,60,10),labels = seq(-60,60,10),limits=c(-60,60))+
  scale_x_continuous(breaks=seq(-100,140,20),labels = seq(-100,140,20),limits=c(-100,140))+
  stat_ellipse(type='norm',linetype=1,linewidth=1,alpha=0.5,level = 0.75,show.legend = F)+
  scale_color_manual(values=c(color_samp,'QC'='yellow'))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'inside',
        legend.position.inside = c(0.85,0.15),
        plot.title = element_text(size=14,face='bold',hjust=0.5)
        )+
  labs(x='Component 1(35.6%)',y='Component 2(10.5%)',color=NULL)+
  ggtitle('Score (PLS-DA) plot')
dev.off()
################################################################################
# 疾病vs对照组间差异代谢产物火山图展示
#de_res[de_res$Metabolite%in% c('3-Hydroxybutanoic Acid','Isocitric Acid'),]
# 绘图 logFC 和 pvalue 不加 VIP值
library(ggplot2)
library(ggrepel)


de_res$regulate2 = factor(de_res$regulate2,levels=c('DCM_Down','NoSig','DCM_Up'))
aa = table(de_res$regulate2)
aa = as.data.frame(aa)
show_metas =c('3-Hydroxybutanoic Acid','Isocitric Acid',
              'Paxilline','Lefamulin','NNAL-N-glucuronide')

pdf("03.meta.DCM/3_2_Volcano_full.pdf",width = 4,height = 4)
ggplot(de_res,aes(x=logFC,y=-log10(P_value)))+
  geom_point(size=1.5,aes(color=regulate2))+
  labs(x='log2(Fold Change)',y="-log10(Pvalue)")+
  geom_vline(xintercept = c(-0.2,0.2),lty=2,col ="grey",lwd=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=2,col = "grey",lwd=0.8)+ 
  theme_test()+theme(
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    legend.position = 'top',
    axis.text = element_text(size=12),
    axis.title=element_text(size=12),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(1.5,'mm'),
    legend.text = element_text(size=10)
  )+scale_color_manual(values = c('DCM_Down'= color_samp[['Control']], 'NoSig'="grey52",'DCM_Up' =color_samp[['DCM']]),
                       breaks=c('DCM_Down','NoSig','DCM_Up'),
                       labels=paste0(aa$Var1,':',aa$Freq))+
  geom_text_repel(data=de_res[de_res$Metabolite%in% show_metas,],
                  aes(x=logFC,y=-log10(P_value),label=Metabolite),
                  show.legend = F,size=3)+
  geom_point(data=de_res[de_res$Metabolite%in% show_metas,],
             aes(x=logFC,y=-log10(P_value)),color='black',size=2.5,shape=21,show.legend = F)+
  scale_x_continuous(breaks = seq(-4,4,1),labels=seq(-4,4,1),limits=c(-3,3))+
  scale_y_continuous(breaks=seq(1,9,1),labels=seq(1,9,1),limits=c(0,9))
dev.off()
#----------------------------------------------------------------------
# 绘图 logFC 和 pvalue 加 VIP值 用作点的大小
pdf("03.meta.DCM/3_2_Volcano_VIP.pdf",width = 5.5,height = 4)
ggplot(de_res,aes(x=logFC,y=-log10(P_value)))+
  geom_point(aes(color=regulate2,size=Vip_oplsda*0.8))+
  labs(x='log2(Fold Change)',y="-log10(Pvalue)",size='VIP')+
  geom_vline(xintercept = c(-0.2,0.2),lty=2,col ="grey",lwd=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=2,col = "grey",lwd=0.8)+ 
  theme_test()+theme(
    axis.text = element_text(size=12),
    axis.title=element_text(size=12),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(1.5,'mm'),
    legend.text = element_text(size=10)
  )+scale_color_manual(values = c('DCM_Down'= color_samp[['Control']], 'NoSig'="grey52",'DCM_Up' =color_samp[['DCM']]),
                       breaks=c('DCM_Down','NoSig','DCM_Up'),
                       labels=paste0(aa$Var1,':',aa$Freq))+
  geom_text_repel(data=de_res[de_res$Metabolite%in% show_metas,],
                  aes(x=logFC,y=-log10(P_value),label=Metabolite),
                  show.legend = F,size=3)+
  geom_point(data=de_res[de_res$Metabolite%in% show_metas,],
             aes(x=logFC,y=-log10(P_value),size=Vip_oplsda+0.8),color='black',shape=21,show.legend = F)+
  scale_x_continuous(breaks = seq(-4,4,1),labels=seq(-4,4,1),limits=c(-3,3))+
  scale_y_continuous(breaks=seq(1,9,1),labels=seq(1,9,1),limits=c(0,9))+
  guides(color=guide_legend(title=NULL))
dev.off()
################################################################################
# 疾病vs对照组间差异代谢产物热图展示
# 热图
exp_use_scale = apply(exp_data,1,scale)
exp_use_scale = t(exp_use_scale)
colnames(exp_use_scale) = colnames(exp_data)

de_res=de_res[order(de_res$logFC,decreasing = T),]
demetas = c(
  de_res[1:20,'Metabolite'],
  tail(de_res[,'Metabolite'],20)
)
heat_data = exp_use_scale[demetas,]

group2 = data.frame(
  row.names = colnames(heat_data),
  group=group
)


pdf('03.meta.DCM/3_3_heatmap.demeta.pdf',height = 6,width = 9)
pheatmap::pheatmap(heat_data,scale='none',
                   annotation_col = group2,
                   cluster_cols = F,
                   annotation_colors=list(group=color_samp),
                   treeheight_row = 15,
                   show_colnames = F,
                   color = colorRampPalette(c("navy","white", "firebrick3"))(50)
)
dev.off()
#colours = colorRampPalette(c("navy", "white", "firebrick3"))(50)
#-------------------------------------------
# complexheatmap 绘图
library(ComplexHeatmap)
library(circlize)
heat_data2 = exp_use_scale[c('Isocitric Acid','3-Hydroxybutanoic Acid',demetas),]
col_fun = colorRampPalette(c("#191970","white", "#8B0000"))(50)

show_meta= unique(c(demetas[nchar(demetas)<15],'Isocitric Acid','3-Hydroxybutanoic Acid'))

col_labels = rep('black',length(show_meta))
col_labels[15:16]='#D96244'
# at 是在原矩阵中的行索引。
anno_meta = rowAnnotation(
  foo = anno_mark(
    at=match(show_meta,rownames(heat_data2)),
    labels=show_meta,
    labels_gp = gpar(cex=0.6,col=col_labels)
  )
)# important_meta

df = data.frame(
  group=group,
  row.names = colnames(heat_data2)
)

top_anno = HeatmapAnnotation(
  df=df,
  show_annotation_name = F,
  col=list(group=color_samp),
  annotation_legend_param = list(title=NULL)
)

ht = Heatmap(as.matrix(heat_data2),
             cluster_rows = T,
             column_km = 2,
             column_title = NULL,
             row_title = NULL,
             show_parent_dend_line = F,
             row_km = 2,
             row_names_max_width = unit(3,'cm'),
             column_dend_height = unit(0.5,'cm'), # 设置列聚类树的高度
             row_dend_width = unit(0.5,'cm'), # 设置行聚类树的宽度
             show_row_names = F,
             right_annotation = anno_meta,
             top_annotation = top_anno,
             show_column_names = F,
             heatmap_legend_param = list(title=NULL),
             col=col_fun)


pdf('03.meta.DCM/3_3_heatmap_demeta_ComplexHeatmap.pdf',width = 5,height = 4)
draw(ht)
dev.off()

################################################################################
# 疾病vs对照组间β羟基丁酸、Isocitric Acid差异水平箱线图
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

fig_data_3_4 = exp_data[rownames(exp_data) %in% c('3-Hydroxybutanoic Acid','Isocitric Acid'),]
fig_data_3_4 = cbind(meta=rownames(fig_data_3_4),fig_data_3_4)
fig_data_3_4=melt(fig_data_3_4,value.name ='exp' ,id.vars = 'meta',variable.name = 'sample')
fig_data_3_4$group = gsub('_.*$','',fig_data_3_4$sample)
fig_data_3_4$group = factor(fig_data_3_4$group,levels=c('DCM','Control'))


pdf('03.meta.DCM/3_4_meta_exp_compare.pdf',width = 4.5,height = 3.5)
ggboxplot(fig_data_3_4,x='meta',y='exp',color = 'group',add='jitter',width = 0.5)+
  scale_color_manual(values=color_samp)+
  theme(axis.title.x = element_blank())+
  stat_compare_means(aes(group=group),label.y=7.5,method='anova',label='p.signif')+
  labs(y='Expression')
dev.off()
################################################################################
# 疾病vs对照组间差异代谢产物通路富集分析
# 注释数据
# 01.data\代谢-小鼠(MJ20240401212) - 10 - 1\interaction_results\1.Preprocess\Metab_table_20240714_103500306\mix\metab_desc.txt
anno_meta = read.table('01.data/代谢-小鼠(MJ20240401212) - 10 - 1/interaction_results/1.Preprocess/Metab_table_20240714_103500306/mix/metab_desc.txt',sep='\t',header=T,quote='"')
DCM_up_Metabolites = de_res[de_res$regulate2=='DCM_Up','Metabolite']
# logFC 筛选条件
ko_ids = anno_meta[anno_meta$Metabolite%in% DCM_up_Metabolites,'KEGG.Compound.ID']
ko_ids = ko_ids[ko_ids!='-']
ko_ids = unique(unlist(strsplit(ko_ids,split=';')))
write(ko_ids,file='03.meta.DCM/3_5_DCM_up_KO_ids.txt',ncolumns = 1)
# ----画图
data_3_5 = read.table('03.meta.DCM/3_5_DCM_UP_Download/msea_ora_result.csv',sep=',',header = T,row.names = 1)
data_3_5 = cbind(path=rownames(data_3_5),data_3_5)


data_3_5$path = factor(data_3_5$path,levels=rev(data_3_5$path))


pdf('03.meta.DCM/3_5_KEGG_enrich_DCM_UP_META.pdf',width = 5,height = 2)
ggplot(data_3_5,aes(x=path,y=-log2(Raw.p)))+
  geom_bar(stat='identity',width=0.3,fill='#8B0000')+coord_flip()+
  theme_test()+
  scale_y_continuous(expand = c(0,0),limits=c(0,7))+
  theme(axis.title.y = element_blank(),
        axis.title.x=element_text(size=14),
        plot.title = element_text(size=14,hjust=0.5,face='bold'),
        axis.text=element_text(size=12))+
  geom_hline(yintercept = -log2(0.05),color='grey',linetype='dotted',linewidth=1)+
  ggtitle('KEGG Pathway')+
  labs(y='-log2(pvalue)')
dev.off()

################################################################################

save.image(paste0('03.meta.DCM.',Sys.Date(),'.RData'))
################################################################################
# de_res[de_res$Metabolite %in% c('3-Hydroxybutanoic Acid','Isocitric Acid'),]

################################################################################


