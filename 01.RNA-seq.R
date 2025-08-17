################################################################################
setwd('03.results')
# db/ dbm
# db = DCM
# dbm = Control
################################################################################
# RNAseq表达矩阵 及差异分析结果
# 预处理数据
data_path = "01.data\\小鼠(MJ20240226006)-6-RNA-seq\\interaction_results\\01 Diff_Express\\DiffExpress_G_RSEM_DESeq2_20240312_111343\\dbdb_vs_dbm.deseq2.xls"
data_rnaseq = read.table(data_path,header=T,sep='\t',comment.char = '@',quote='"') 
de_res = data_rnaseq[,c(2,6:10)]
count_data = data_rnaseq[,c(2,11:16)]
tpm_data = data_rnaseq[,c(2,17:22)]

# 处理差异结果
colnames(de_res)[2]='log2FC'
de_res$significant=NULL
de_res$regulate = ifelse(
  de_res$log2FC>= log2(1.5) & de_res$padjust<0.05,'up',
    ifelse(de_res$log2FC <= log2(2/3) & de_res$padjust<0.05,'down','nosig')
)
de_res$pvalue=NULL
de_res = de_res[de_res$gene_name!='',]
# 处理表达数据
colnames(count_data) = gsub('dbm','Control',colnames(count_data))
colnames(count_data) = gsub('db','DCM',colnames(count_data))
colnames(count_data) = gsub('_count','',colnames(count_data))
colnames(tpm_data) = gsub('dbm','Control',colnames(tpm_data))
colnames(tpm_data) = gsub('db','DCM',colnames(tpm_data))
colnames(tpm_data) = gsub('_tpm','',colnames(tpm_data))
count_data = count_data[count_data$gene_name!='',]
tpm_data = tpm_data[tpm_data$gene_name!='',]

rownames(count_data) = make.unique(count_data$gene_name)
count_data$gene_name=NULL
rownames(tpm_data) = make.unique(tpm_data$gene_name)
tpm_data$gene_name = NULL

# 保存处理结果
if(!dir.exists('02.RNAseq')){
  dir.create('02.RNAseq')
}
write.table(de_res,file='02.RNAseq/00.de_results.txt',sep='\t',row.names = F,quote=F)
write.table(count_data,file='02.RNAseq/00.count_data.txt',sep='\t',quote=F)
write.table(tpm_data,file='02.RNAseq/00.tpm_data.txt',sep='\t',quote=F)

group = c(rep('DCM',3),rep('Control',3))
color_samp = c('DCM'='#D96244','Control'='#5EB5C8')
################################################################################
#  热图
library(pheatmap)
fig_data_2_1 = log2(count_data+1)
dat=apply(fig_data_2_1,1,scale)
rownames(dat)=colnames(fig_data_2_1)
dat=t(dat)
degs = de_res[de_res$regulate!='nosig','gene_name']
dat = dat[degs,]

anno_col=data.frame(
  row.names = colnames(dat),
  group=group
)
anno_colors=list(
  group=color_samp
)

mat_dist = dist(dat)
hclust_1 = hclust(mat_dist)
test=rev(hclust_1$labels[hclust_1$order])

pdf('02.RNAseq/fig_2_1_heatmap.pdf',width=4,height=2.5)
pheatmap(dat[test,],scale='none',
         show_rownames = F,cluster_rows = F,
         show_colnames = F,
         cluster_cols = F,
         cellwidth = 15,
         annotation_col = anno_col,annotation_colors = anno_colors,
         annotation_names_col = F,
         color=colorRampPalette(c('#191970','white','#FF4500'))(100),
         border_color = NA
)
dev.off()
################################################################################
# 差异基因 GO  KEGG富集
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
#------------------------------------------------------
DCM_UP_genes = unique(de_res[de_res$regulate=='up','gene_name'])
DCM_UP_genes_ids=bitr(DCM_UP_genes,fromType = 'SYMBOL',toType ='ENTREZID',OrgDb = 'org.Mm.eg.db')

DCM_DOWN_genes = unique(de_res[de_res$regulate=='down','gene_name'])
DCM_DOWN_genes_ids=bitr(DCM_DOWN_genes,fromType = 'SYMBOL',toType ='ENTREZID',OrgDb = 'org.Mm.eg.db')


# GO 富集 up
ego_ALL_up <- enrichGO(gene = as.numeric(DCM_UP_genes_ids$ENTREZID), 
                       OrgDb = 'org.Mm.eg.db',
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2,
                       readable = TRUE) 
ego_ALL_up=as.data.frame(ego_ALL_up)
write.table(ego_ALL_up,file='02.RNAseq/2_2_enrich_GO_DCM_UP.txt',sep='\t',quote=F,row.names = F)
# GO 富集 down
ego_ALL_down <- enrichGO(gene = as.numeric(DCM_DOWN_genes_ids$ENTREZID),
                       OrgDb = 'org.Mm.eg.db',
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
ego_ALL_down=as.data.frame(ego_ALL_down)
write.table(ego_ALL_down,file='02.RNAseq/2_2_enrich_GO_DCM_DOWN.txt',sep='\t',quote=F,row.names = F)


# KEGG 富集 up
enrich_kegg_up <- enrichKEGG(gene = DCM_UP_genes_ids$ENTREZID,
                             organism = 'mmu',
                             pvalueCutoff = 0.05,
                             pAdjustMethod = 'none',
)
enrich_kegg_up = DOSE::setReadable(enrich_kegg_up,OrgDb='org.Mm.eg.db',keyType = 'ENTREZID')
enrich_kegg_up = as.data.frame(enrich_kegg_up)
write.table(enrich_kegg_up,file='02.RNAseq/2_2_enrich_KEGG_DCM_UP.txt',sep='\t',quote=F,row.names = F)

# KEGG 富集 down
enrich_kegg_down <- enrichKEGG(gene = DCM_DOWN_genes_ids$ENTREZID,
                             organism = 'mmu',
                             pvalueCutoff = 0.05,
                             pAdjustMethod = 'none',
)
enrich_kegg_down = DOSE::setReadable(enrich_kegg_down,OrgDb='org.Mm.eg.db',keyType = 'ENTREZID')
enrich_kegg_down = as.data.frame(enrich_kegg_down)
write.table(enrich_kegg_down,file='02.RNAseq/2_2_enrich_KEGG_DCM_DOWN.txt',sep='\t',quote=F,row.names = F)


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
plot_top_go_up = rbind(plot_top_go_up,go_fatty_relation)
plot_top_go_up = plot_top_go_up[!plot_top_go_up$ID%in% c('GO:0046873','GO:0015081','GO:0005244'),]
plot_top_go_up = rbind(plot_top_go_up,
                       ego_ALL_up[ego_ALL_up$ID%in% c('GO:0071617','GO:0008374'),]
                       )
plot_top_go_up$Description2 = capitalize(plot_top_go_up$Description)
plot_top_go_up=plot_top_go_up %>%
  arrange(ONTOLOGY,desc(Count))
plot_top_go_up$Description2=factor(plot_top_go_up$Description2,levels=rev(plot_top_go_up$Description2))

pdf('02.RNAseq/2_2_enrich_GO_DCM_UP_top5.pdf',width = 6,height = 3.5)
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
# ----------------------------------------
# GO 可视化 down
library(dplyr)
library(stringr)
library(Hmisc)
library(ggplot2)
plot_top_go_down = ego_ALL_down %>% 
  group_by(ONTOLOGY) %>% 
  slice_head(n=5)

plot_top_go_down$Description2 = capitalize(plot_top_go_down$Description)
plot_top_go_down=plot_top_go_down %>%
  arrange(ONTOLOGY,desc(Count))
plot_top_go_down$Description2=factor(plot_top_go_down$Description2,levels=rev(plot_top_go_down$Description2))

pdf('02.RNAseq/2_2_enrich_GO_DCM_DOWN_top5.pdf',width = 5.5,height = 3.5)
ggplot(plot_top_go_down,aes(x=Count,y=Description2,color=-log10(p.adjust),shape=ONTOLOGY))+
  geom_point(size=3)+xlim(c(min(plot_top_go_down$Count)-0.5,max(plot_top_go_down$Count)+0.5))+theme_bw()+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.5,'cm'),
        axis.text.y=element_text(size=12))+
  labs(shape='Ontology',color=bquote(-log[10]~(FDR)),x='Count')+
  scale_y_discrete(position = 'right',label=function(x){stringr::str_wrap(x,width = 60)} )+
  scale_color_gradient(low='blue',high='red')+
  guides(shape=guide_legend(order=1))
dev.off() 

#------------------------------------------
# KEGG可视化 up
library(tidyverse)
library(cowplot)
plot_up_kegg = head(enrich_kegg_up,15)
plot_up_kegg$Description2=gsub('\\s-\\sMus\\smusculus\\s\\(house\\smouse\\)','',plot_up_kegg$Description)

plot_up_kegg = plot_up_kegg[order(-log10(plot_up_kegg$pvalue),decreasing = F),]
plot_up_kegg$Description2 = factor(plot_up_kegg$Description2,levels=plot_up_kegg$Description2)
label_col = rep('black',15)
label_col[grepl('fatty',plot_up_kegg$Description2,ignore.case = T)]='red'

pdf('02.RNAseq/2_2_enrich_KEGG_DCM_UP_top15.pdf',width = 5.5,height = 3.5)
ggplot(plot_up_kegg,aes(x=-log10(pvalue),y=Description2,color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+
  xlim(sort(-log10(plot_up_kegg$pvalue))[c(1,15)]+c(-1,0.5))+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.4,'cm'),
        axis.text.y=element_text(size=12,color=label_col))+
  labs(x=bquote(-log[10]~(pvalue)))+
  scale_y_discrete(position = 'right')+
  scale_color_gradient(low='blue',high='red')+
  guides(color='none')
dev.off()
#--------------------------------------------------------
# KEGG可视化 down
library(tidyverse)
library(cowplot)
plot_down_kegg = enrich_kegg_down
plot_down_kegg = plot_down_kegg[plot_down_kegg$Count>2,]
plot_down_kegg = head(plot_down_kegg,15)
plot_down_kegg$Description2=gsub('\\s-\\sMus\\smusculus\\s\\(house\\smouse\\)','',plot_up_kegg$Description)

plot_down_kegg = plot_down_kegg[order(-log10(plot_down_kegg$pvalue),decreasing = F),]
plot_down_kegg$Description2 = factor(plot_down_kegg$Description2,levels=plot_down_kegg$Description2)

pdf('02.RNAseq/2_2_enrich_KEGG_DCM_DOWN_top15.pdf',width = 5.5,height = 3.5)
ggplot(plot_down_kegg,aes(x=-log10(pvalue),y=Description2,color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+
  xlim(sort(-log10(plot_down_kegg$pvalue))[c(1,15)]+c(-0.5,0.5))+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.4,'cm'),
        axis.text.y=element_text(size=12))+
  labs(x=bquote(-log[10]~(pvalue)))+
  scale_y_discrete(position = 'right')+
  scale_color_gradient(low='blue',high='red')+
  guides(color='none')
dev.off()
################################################################################
# 基因表达展示
de_res[de_res$gene_name%in%c('Cpt1a','Hmgcs2','Acsm5'),]

fig_data_2_3 = tpm_data[rownames(tpm_data) %in% c('Cpt1a','Hmgcs2','Acsm5'),]
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

fig_data_2_3 = cbind(gene=rownames(fig_data_2_3),fig_data_2_3)
fig_data_2_3=melt(fig_data_2_3,value.name ='tpm' ,id.vars = 'gene',variable.name = 'sample')
fig_data_2_3$gene=factor(fig_data_2_3$gene,levels=c('Acsm5','Cpt1a','Hmgcs2'))
fig_data_2_3$group = gsub('_.*$','',fig_data_2_3$sample)
fig_data_2_3$group = factor(fig_data_2_3$group,levels=c('DCM','Control'))
fig_data_2_3$logtpm = log2(fig_data_2_3$tpm+1)

pdf('02.RNAseq/2_3_gene_exp_compare.pdf',width = 5,height = 3.5)
ggboxplot(fig_data_2_3,x='gene',y='logtpm',color = 'group',add='jitter')+
  scale_color_manual(values=color_samp)+
  theme(axis.title.x = element_blank())+
  stat_compare_means(aes(group=group),label.y=6.5,method='anova',label='p.signif')+
  labs(y='Expression log2(TPM)')
dev.off()
# p.signif p.format
################################################################################
# GSVA 
# 获得β羟基丁酰化(Kbhb)相关过程通路
library(msigdbr)
library(dplyr)

# 获取小鼠的 基因集
gobp_set=msigdbr(species = 'Mus musculus',category = 'C5',subcategory = 'GO:BP')
pathway_set = msigdbr(species = 'Mus musculus',category = 'C2')
all_pathway_bp_process = rbind(gobp_set,pathway_set)

Kbhb_key_words="_FATTY_ACID_OXIDATION|FATTY_ACID_BETA_OXIDATION|GLUCONEOGENESIS|GLYCOLYSIS|CITRATE_CYCLE|TCA_CYCLE|_KETONE_"
Kbhb_relation_set=all_pathway_bp_process[grepl(Kbhb_key_words,all_pathway_bp_process$gs_name),]
path_set=pathway_set[grepl('_FATTY_ACID',pathway_set$gs_name),]
path_set = path_set[path_set$gs_subcat=='CP:KEGG',]

Kbhb_relation_set = rbind(Kbhb_relation_set, path_set)
Kbhb_relation_set= Kbhb_relation_set %>%
  dplyr::select(gs_name, gene_symbol)
Kbhb_relation_set = as.data.frame(Kbhb_relation_set)
kbhb_set = split(Kbhb_relation_set$gene_symbol, Kbhb_relation_set$gs_name)

# 获取人的基因集
gobp_set=msigdbr(species = 'Homo sapiens',category = 'C5',subcategory = 'GO:BP')
pathway_set = msigdbr(species = 'Homo sapiens',category = 'C2')
all_pathway_bp_process = rbind(gobp_set,pathway_set)

Kbhb_relation_set=all_pathway_bp_process[grepl(Kbhb_key_words,all_pathway_bp_process$gs_name),]
path_set=pathway_set[grepl('_FATTY_ACID',pathway_set$gs_name),]
path_set = path_set[path_set$gs_subcat=='CP:KEGG',]

Kbhb_relation_set = rbind(Kbhb_relation_set, path_set)
Kbhb_relation_set= Kbhb_relation_set %>%
  dplyr::select(gs_name, gene_symbol)
Kbhb_relation_set = as.data.frame(Kbhb_relation_set)

kbhb_set_human = split(Kbhb_relation_set$gene_symbol, Kbhb_relation_set$gs_name)

#----GSVA分析
library(GSVA)

data_2_4=log2(tpm_data+1)
gp=gsvaParam(as.matrix(data_2_4), kbhb_set,
             maxDiff =FALSE,
             kcdf='Gaussian',
             minSize = 3)
res_2_4_gsva <- gsva(gp,BPPARAM=BiocParallel::SnowParam(workers=6))
res_2_4_gsva = as.data.frame(res_2_4_gsva)

# 绘图展示GSVA score
library(pheatmap)
pdf('02.RNAseq/2_4_kbhb_gsva_score_heatmap.pdf',width = 10,height = 8)
pheatmap(res_2_4_gsva, show_colnames = T, 
         scale = "none",angle_col = "45",
         cluster_row = T,cluster_col = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

#找差异GSVA得分条目
library(limma)
group_list=factor(group,levels=c('DCM','Control'))

design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(res_2_4_gsva)
contrast.matrix<-makeContrasts("DCM-Control",levels = design)

fit <- lmFit(res_2_4_gsva,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯
DE_term <- topTable(fit2, coef = 1,n=Inf,adjust.method = 'none')

# 输出GSVA表
res_2_4_gsva1=cbind(pathway=rownames(res_2_4_gsva),res_2_4_gsva)
res_2_4_gsva1=res_2_4_gsva1[order(res_2_4_gsva1$pathway),]
res_2_4_gsva2=cbind(res_2_4_gsva1,DE_term[match(rownames(res_2_4_gsva1),rownames(DE_term)),])
res_2_4_gsva2$B=NULL
res_2_4_gsva2$t=NULL

wb = openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,'GSVA_Kbhb')
openxlsx::writeData(wb,'GSVA_Kbhb',x=res_2_4_gsva2)
openxlsx::addFilter(wb,1,1,cols=1:ncol(res_2_4_gsva2))
nwid=rep(15,ncol(res_2_4_gsva2))
nwid[1]=50
openxlsx::setColWidths(wb,1,cols=1:ncol(res_2_4_gsva2),widths = nwid)
openxlsx::freezePane(wb,1,firstRow = T)
openxlsx::saveWorkbook(wb,'02.RNAseq/2_4_GSVA_pathway_table_RNAseq.xlsx',overwrite = T)

# 显著的过程通路
DE_set <- DE_term[DE_term$logFC > 0.3,]
dat <- res_2_4_gsva[match(rownames(DE_set),rownames(res_2_4_gsva)),]
library(ComplexHeatmap)
  
top_anno = HeatmapAnnotation(
  group=group,
  annotation_legend_param = list(
    title='group'
  ),
  col=list(group=color_samp),
  annotation_name_gp = gpar(col='black')
)
ht = Heatmap(as.matrix(dat),col=colorRampPalette(c("navy", "white", "firebrick3"))(50),row_dend_width = unit(5,'mm'),top_annotation = top_anno,show_column_names = F,cluster_columns = F,heatmap_legend_param = list(title=""),row_names_max_width =unit(20,'cm'))

pdf('02.RNAseq/2_4_kbhb_gsva_score_heatmap_sig.pdf',width = 11,height = 2.5)
draw(ht,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()
#------------------------------------ GSEA
library(clusterProfiler)
library(enrichplot)
enrich_data_2_4 = de_res$log2FC
names(enrich_data_2_4) = de_res$gene_name
enrich_data_2_4 = sort(enrich_data_2_4,decreasing = T)
enrich_data_2_4=enrich_data_2_4[enrich_data_2_4!=0]

#GSEA 富集分析
enrich_2_4_kbhb = GSEA(geneList = enrich_data_2_4,TERM2GENE = Kbhb_relation_set,seed=123,pvalueCutoff = 0.05,minGSSize = 3,maxGSSize = 200,pAdjustMethod = 'none')
#gseaplot2(enrich_2_4_kbhb,geneSetID = enrich_2_4_kbhb@result$Description)
write.table(enrich_2_4_kbhb@result,file='02.RNAseq/2_4_GSEA_kbhb_sig.txt',sep='\t',quote=F,row.names = F)

#  火山图-----------------------------------------------------------
library(ggrepel)
label_genes = unique(c('Cpt1a','Hmgcs2','Acsm5',Kbhb_relation_set[Kbhb_relation_set$gs_name %in% rownames(DE_set),'gene_symbol']))

aa = table(de_res$regulate)
names(aa)=c('DCM_Down','No_Sig','DCM_Up')
aa = as.data.frame.table(aa)

labeldata = de_res[de_res$regulate!='nosig',]
label_genes = label_genes[label_genes%in% labeldata$gene_name]
labeldata = labeldata[labeldata$gene_name %in% label_genes,]


pdf("02.RNAseq/2_1_volcano.pdf",width = 8,height = 5)
ggplot(de_res,aes(x=log2FC,y=-log10(padjust)))+
  geom_point(size=2,aes(color=regulate))+
  labs(x='log2(Fold Change)',y="-log10(padjust)")+
  geom_vline(xintercept = c(log2(2/3),log2(1.5)),lty=2,col ="grey30",lwd=0.5)+
  geom_hline(yintercept=-log10(0.05),lty=2,col = "grey30",lwd=0.5)+ 
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    legend.position = 'top',
    axis.text = element_text(size=12),
    axis.title=element_text(size=14),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(2,'mm'),
    legend.text = element_text(size=12)
  )+scale_color_manual(values = c('down'=color_samp[['Control']], 'nosig'="grey",'up' =color_samp[['DCM']]),
                       breaks=c('down','nosig','up'),
                       labels=paste0(aa$Var1,':',aa$Freq))+
  geom_text_repel(data=labeldata,aes(x=log2FC,y=-log10(padjust),label=gene_name),size=3.5,max.overlaps = 30)+
  geom_point(data=labeldata,aes(x=log2FC,y=-log10(padjust)),shape=21,size=2.5,col='red')

dev.off()


################################################################################
# 公开数据表达比较
de_res[de_res$gene_name%in%c('Hdac1','Hdac2','Hdac6','Ep300'),]
library(reshape2)
library(ggpubr)
library('homologene')
human_gene=homologene(c('Hdac1','Hdac2','Hdac6','Ep300'),inTax=10090,outTax=9606)
# 表达数据处理
#-----------------------------------------
# GSE197850（DCM vs CM) 处理表达居住 人的
data_2_5_GSE197850.tpm = read.table('02.RNAseq/GEO_DATA/GSE197850_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz',header = T,sep='\t')
data_2_5_GSE197850.tpm = data_2_5_GSE197850.tpm[,c('GeneID',paste0('GSM59313',65:72))]
colnames(data_2_5_GSE197850.tpm)[2:9]=c(paste0('CM',1:4),paste0('DCM',1:4))

ids_2_5 = data_2_5_GSE197850.tpm$GeneID
library(clusterProfiler)
library(org.Hs.eg.db)
ids_2_name = bitr(ids_2_5,fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
data_2_5_GSE197850.tpm$genename=ids_2_name[match(data_2_5_GSE197850.tpm$GeneID,ids_2_name$ENTREZID),'SYMBOL']
data_2_5_GSE197850.tpm = data_2_5_GSE197850.tpm[!is.na(data_2_5_GSE197850.tpm$genename),]
data_2_5_GSE197850.tpm = data_2_5_GSE197850.tpm[!duplicated(data_2_5_GSE197850.tpm$genename),]
rownames(data_2_5_GSE197850.tpm)=data_2_5_GSE197850.tpm$genename
data_2_5_GSE197850.tpm$genename=NULL
data_2_5_GSE197850.tpm$GeneID=NULL
group_1 = c(rep('CM',4),rep('DCM',4))

#-----------------------------------------
# GSE211106（DCM vs cotrol） 小鼠
data_2_5_GSE211106.fpkm = read.table('02.RNAseq/GEO_DATA/GSE211106_RNA-seq_fpkm.csv.gz',sep=',',header = T)
colnames(data_2_5_GSE211106.fpkm)[2:7]=c(paste0('Control_',1:3),paste0('DCM',1:3))
rownames(data_2_5_GSE211106.fpkm)=data_2_5_GSE211106.fpkm$GeneID
data_2_5_GSE211106.fpkm$GeneID=NULL
group_2 = c(rep('Control',3),rep('DCM',3))
#-----------------------------------------
# GSE161931（Db/db vs Bks）小鼠
data_2_5_GSE161931.fpkm = read.table('02.RNAseq/GEO_DATA/GSE161931_geneFpkm.csv.gz',sep=',',header = T)
colnames(data_2_5_GSE161931.fpkm)[2:11]=c(paste0('Bks_',1:5),paste0('DCM_',1:5))

library(biomaRt)
mart = useMart('ENSEMBL_MART_ENSEMBL',dataset = 'mmusculus_gene_ensembl',host = 'https://sep2019.archive.ensembl.org')
all_id = data_2_5_GSE161931.fpkm$gene_id
id_trans = getBM(attributes = c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',values = all_id,mart = mart)
common_id = intersect(id_trans$ensembl_gene_id,all_id)
id_trans = id_trans[match(common_id,id_trans$ensembl_gene_id),]
data_2_5_GSE161931.fpkm = data_2_5_GSE161931.fpkm[match(common_id,data_2_5_GSE161931.fpkm$gene_id),]

identical(data_2_5_GSE161931.fpkm$gene_id,id_trans$ensembl_gene_id)
data_2_5_GSE161931.fpkm$gene_name=id_trans$mgi_symbol
data_2_5_GSE161931.fpkm$gene_id=NULL
data_2_5_GSE161931.fpkm = data_2_5_GSE161931.fpkm[data_2_5_GSE161931.fpkm$gene_name!='',]

rownames(data_2_5_GSE161931.fpkm) = make.unique(data_2_5_GSE161931.fpkm$gene_name)
data_2_5_GSE161931.fpkm$gene_name=NULL
group_3 = c(rep('Bks',5),rep('DCM',5))
#-----------------------------------------
# GSE244904（dbdb AAV9 null vs control） 小鼠
data_2_5_GSE244904 = read.table('02.RNAseq/GEO_DATA/GSE244904_Expressed_gene_reads.txt.gz',sep='\t',header=T)
data_2_5_GSE244904[,2:6]=NULL
data_2_5_GSE244904[,8:13]=NULL
colnames(data_2_5_GSE244904)[2:4]=paste0('DCM_',1:3)
data_2_5_GSE244904 = data_2_5_GSE244904[!duplicated(data_2_5_GSE244904$Geneid),]
rownames(data_2_5_GSE244904)=data_2_5_GSE244904$Geneid
data_2_5_GSE244904$Geneid=NULL
group_4 = c(rep('DCM',3),rep('Control',3))

#-----------------------------------------
color_samp2=c(color_samp,'CM'='#5EB5C8','Bks'='#5EB5C8')

fig_2_5_plot = function(exp_data,curgene,curgroup,filename,compare_methods='anova',ytitle='Expression'){
  tmp_data = exp_data[rownames(exp_data) %in% curgene,]
  tmp_data = cbind(gene=rownames(tmp_data),tmp_data)
  tmp_data = reshape2::melt(tmp_data,id.vars='gene',value.name = 'exp',variable.name='sample')
  anno_sample=data.frame(
    sample=colnames(exp_data),
    group=curgroup
  )
  tmp_data$group=anno_sample[match(tmp_data$sample,anno_sample$sample),'group']
  tmp_data$sample=NULL
  tmp_data$gene = factor(tmp_data$gene,levels=sort(curgene))
  require(ggpubr)
  require(ggplot2)
  tmp_data$logexp = log2(tmp_data$exp+1)
  top_y = max(tmp_data$logexp)+1
  
  fig_pl = ggboxplot(tmp_data,x='gene',y='logexp',color = 'group',add='jitter')+
    scale_color_manual(values=color_samp2)+
    theme(axis.title.x = element_blank())+
    stat_compare_means(aes(group=group),label.y=top_y,method=compare_methods,label='p.signif')+
    labs(y=ytitle)
  
  if(!dir.exists(dirname(filename))){
    dir.create(dirname(filename),recursive = T)
  }
  pdf(filename,width=6,height = 4)
  print(fig_pl)
  dev.off()
}
#------------
fig_2_5_plot(exp_data = data_2_5_GSE197850.tpm,curgene = c('HDAC1','HDAC2','HDAC6','EP300'),curgroup = group_1,filename = '02.RNAseq/2_5_GSE197850/kbhb_key_enzyme_exp_in_DCM.pdf',ytitle = 'log2(TPM Expression)')
fig_2_5_plot(exp_data = data_2_5_GSE211106.fpkm,curgene = c('Hdac1','Hdac2','Hdac6','Ep300'),curgroup = group_2,filename = '02.RNAseq/2_5_GSE211106/kbhb_key_enzyme_exp_in_DCM.pdf',ytitle = 'log2(FPKM Expression)')
fig_2_5_plot(exp_data = data_2_5_GSE161931.fpkm,curgene = c('Hdac1','Hdac2','Hdac6','Ep300'),curgroup = group_3,filename = '02.RNAseq/2_5_GSE161931/kbhb_key_enzyme_exp_in_DCM.pdf',ytitle = 'log2(FPKM Expression)')
fig_2_5_plot(exp_data = data_2_5_GSE244904,curgene = c('Hdac1','Hdac2','Hdac6','Ep300'),curgroup = group_4,filename = '02.RNAseq/2_5_GSE244904/kbhb_key_enzyme_exp_in_DCM.pdf',ytitle = 'log2(FPKM Expression)')

################################################################################
# 公开数据 疾病 vs对照组间β羟基丁酰化相关通路GSVA活性差异
library(limma)
library(ComplexHeatmap)
fig_2_6_gsva = function(exp_data,gsset,outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  require(GSVA)
  data_2_6=log2(exp_data+1)
  gp1=gsvaParam(as.matrix(data_2_6), gsset,
             maxDiff =FALSE,
             kcdf='Gaussian',
             minSize = 3)
  res_2_6_gsva <- gsva(gp1,BPPARAM=BiocParallel::SnowParam(workers=6))
  res_2_6_gsva = as.data.frame(res_2_6_gsva)
  # 绘图展示GSVA score
  pdf(paste0(outdir,'/gsva_res_full_set_heatmap.pdf'),width = 10,height = 8)
  pheatmap::pheatmap(res_2_6_gsva, show_colnames = T, 
         scale = "none",angle_col = "45",
         cluster_row = T,cluster_col = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  dev.off()
  return(res_2_6_gsva)
}
#---------------------------------------------------------
# GSE197850
gsva_res_GSE197850 = fig_2_6_gsva(exp_data = data_2_5_GSE197850.tpm,gsset=kbhb_set_human,outdir = '02.RNAseq/2_5_GSE197850')

#找差异GSVA得分条目
group_list=factor(group_1,levels=c('CM','DCM'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(gsva_res_GSE197850)
contrast.matrix<-makeContrasts("DCM-CM",levels = design)

fit <- lmFit(gsva_res_GSE197850,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯
DE_term <- topTable(fit2, coef = 1,n=Inf,adjust.method = 'none')

# 显著的过程通路
DE_set_1 <- DE_term[DE_term$logFC > 0.26,]
dat <- gsva_res_GSE197850[match(rownames(DE_set_1),rownames(gsva_res_GSE197850)),]

dat = dat[,c(paste0('DCM',1:4),paste0('CM',1:4))]

top_anno = HeatmapAnnotation(
  group=rev(group_1),
  annotation_legend_param = list(
    title='group'
  ),
  col=list(group=color_samp2),
  annotation_name_gp = gpar(col='black')
)
ht = Heatmap(as.matrix(dat),col=colorRampPalette(c("navy", "white", "firebrick3"))(50),heatmap_height = unit(30,'mm'),row_dend_width = unit(5,'mm'),top_annotation = top_anno,show_column_names = F,cluster_columns = F,heatmap_legend_param = list(title=""),row_names_max_width =unit(20,'cm'))
pdf('02.RNAseq/2_5_GSE197850/2_6_kbhb_gsva_score_heatmap_sig.pdf',width = 11,height = 2)
draw(ht,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()

#---------------------------------------------------------
# GSE211106
gsva_res_GSE211106 = fig_2_6_gsva(exp_data = data_2_5_GSE211106.fpkm,gsset=kbhb_set,outdir = '02.RNAseq/2_5_GSE211106/')
#找差异GSVA得分条目
group_list=factor(group_2,levels=c('DCM','Control'))

design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(gsva_res_GSE211106)
contrast.matrix<-makeContrasts("DCM-Control",levels = design)

fit <- lmFit(gsva_res_GSE211106,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯
DE_term <- topTable(fit2, coef = 1,n=Inf,adjust.method = 'none')

# 显著的过程通路
DE_set_2 <- DE_term[DE_term$logFC > 0.58,]
dat <- gsva_res_GSE211106[match(rownames(DE_set_2),rownames(gsva_res_GSE211106)),]

dat = dat[,c(paste0('DCM',1:3),paste0('Control_',1:3))]

top_anno = HeatmapAnnotation(
  group=rev(group_2),
  annotation_legend_param = list(
    title='group'
  ),
  col=list(group=color_samp2),
  annotation_name_gp = gpar(col='black')
)
ht = Heatmap(as.matrix(dat),col=colorRampPalette(c("navy", "white", "firebrick3"))(50),row_dend_width = unit(5,'mm'),top_annotation = top_anno,show_column_names = F,cluster_columns = F,heatmap_legend_param = list(title=""),row_names_max_width =unit(20,'cm'))

pdf('02.RNAseq/2_5_GSE211106/2_6_kbhb_gsva_score_heatmap_sig.pdf',width = 11,height = 5)
draw(ht,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()
#---------------------------------------------------------
# GSE161931
gsva_res_GSE161931 = fig_2_6_gsva(exp_data = data_2_5_GSE161931.fpkm,gsset=kbhb_set,outdir = '02.RNAseq/2_5_GSE161931/')
#找差异GSVA得分条目
group_list=factor(group_3,levels=c('DCM','Bks'))

design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(gsva_res_GSE161931)
contrast.matrix<-makeContrasts("DCM-Bks",levels = design)

fit <- lmFit(gsva_res_GSE161931,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯
DE_term <- topTable(fit2, coef = 1,n=Inf,adjust.method = 'none')

# 显著的过程通路
DE_set_3 <- DE_term[DE_term$logFC > 0.3,]
dat <- gsva_res_GSE161931[match(rownames(DE_set_3),rownames(gsva_res_GSE161931)),]

dat = dat[,c(paste0('DCM_',1:5),paste0('Bks_',1:5))]


top_anno = HeatmapAnnotation(
  group=rev(group_3),
  annotation_legend_param = list(
    title='group'
  ),
  col=list(group=color_samp2),
  annotation_name_gp = gpar(col='black')
)
ht = Heatmap(as.matrix(dat),col=colorRampPalette(c("navy", "white", "firebrick3"))(50),row_dend_width = unit(5,'mm'),top_annotation = top_anno,show_column_names = F,cluster_columns = F,heatmap_legend_param = list(title=""),row_names_max_width =unit(20,'cm'))

pdf('02.RNAseq/2_5_GSE161931/2_6_kbhb_gsva_score_heatmap_sig.pdf',width = 11,height = 3)
draw(ht,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()
#---------------------------------------------------------
# GSE244904
gsva_res_GSE244904 = fig_2_6_gsva(exp_data = data_2_5_GSE244904,gsset=kbhb_set,outdir = '02.RNAseq/2_5_GSE244904/')
#找差异GSVA得分条目
group_list=factor(group_4,levels=c('DCM','Control'))

design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(gsva_res_GSE244904)
contrast.matrix<-makeContrasts("DCM-Control",levels = design)

fit <- lmFit(gsva_res_GSE244904,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯
DE_term <- topTable(fit2, coef = 1,n=Inf,adjust.method = 'none')

# 显著的过程通路
DE_set_4 <- DE_term[DE_term$logFC > 0.3,]
dat <- gsva_res_GSE244904[match(rownames(DE_set_4),rownames(gsva_res_GSE244904)),]

top_anno = HeatmapAnnotation(
  group=group_4,
  annotation_legend_param = list(
    title='group'
  ),
  col=list(group=color_samp2),
  annotation_name_gp = gpar(col='black')
)
ht = Heatmap(as.matrix(dat),col=colorRampPalette(c("navy", "white", "firebrick3"))(50),row_dend_width = unit(5,'mm'),top_annotation = top_anno,show_column_names = F,cluster_columns = F,heatmap_legend_param = list(title=""),row_names_max_width =unit(20,'cm'))

pdf('02.RNAseq/2_5_GSE244904/2_6_kbhb_gsva_score_heatmap_sig.pdf',width = 11,height = 4)
draw(ht,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()

################################################################################
setwd('E:\\项目\\2.滕雁波\\10.多组学联合分析-修饰\\03.results')
save.image(paste0('02.RNAseq.',Sys.Date(),'.RData'))
################################################################################
