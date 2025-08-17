################################################################################
setwd('03.results')
# db/ dbm
# db = DCM
# dbm = Control
group = c(rep('Control',3),rep('DCM',3))
color_samp = c('DCM'='#D96244','Control'='#5EB5C8')
################################################################################
# 蛋白质表达矩阵 及差异分析结果
# 预处理数据
data_path = "01.data/DIA报告/4-Differentially_expressed_protein/T-test_analysis/Differentially_annot.xlsx"
data_protein_DIA = openxlsx::read.xlsx(data_path,sheet = 'dbvsdbm') 

de_res = data_protein_DIA[,c(1,3:5)]
rel_data = data_protein_DIA[,c(1,3,17:22)]

# 处理差异结果
colnames(de_res)=c('protein','gene','FC','pvalue')
de_res$log2FC = log2(de_res$FC)
de_res$regulate = ifelse(
  de_res$log2FC>= log2(1.5) & de_res$pvalue<0.05,'up',
  ifelse(de_res$log2FC <= log2(2/3) & de_res$pvalue<0.05,'down','nosig')
)
de_res$ids = paste0(de_res$protein,'__',de_res$gene)
de_res=de_res[!is.na(de_res$pvalue),]
de_res$padjust = p.adjust(de_res$pvalue,method = 'fdr')
de_res = de_res[,c(1,2,7,3,5,4,8,6)]


# 处理表达数据
colnames(rel_data) = gsub('dbm','Control_',colnames(rel_data))
colnames(rel_data) = gsub('db','DCM_',colnames(rel_data))
rownames(rel_data) = paste0(rel_data$Protein.accession,'__',rel_data$Gene.name)
rel_data$Protein.accession=NULL
rel_data$Gene.name=NULL

# 保存处理结果
if(!dir.exists('04.protein_dia')){
  dir.create('04.protein_dia')
}

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
sheet1=openxlsx::addWorksheet(wb, sheetName = 'Protein Differential Results')
openxlsx::writeDataTable(wb=wb,sheet=sheet1,x=de_res,bandedRows = T)
nwid=nchar(colnames(de_res))
nwid=nwid+10
openxlsx::setColWidths(wb=wb,sheet=sheet1,cols=1:ncol(de_res),widths=nwid)
openxlsx::freezePane(wb,sheet=sheet1, firstRow = TRUE)
openxlsx::saveWorkbook(wb, '04.protein_dia/00.Differential_results_protein.xlsx', overwrite = TRUE)

write.table(de_res,file='04.protein_dia/00.de_results.txt',sep='\t',row.names = F,quote=F)
write.table(rel_data,file='04.protein_dia/00.rel_data.txt',sep='\t',quote=F)

################################################################################
#  获得脂肪酸β氧化相关基因
library(msigdbr)
library(dplyr)
gobp_set=msigdbr(species = 'Mus musculus',category = 'C5',subcategory = 'GO:BP')
pathway_set = msigdbr(species = 'Mus musculus',category = 'C2')
all_pathway_bp_process = rbind(gobp_set,pathway_set)
Kbhb_relation_set=all_pathway_bp_process[grepl("FATTY_ACID_BETA_OXIDATION",all_pathway_bp_process$gs_name),]
fatty_acid_genes = unique(Kbhb_relation_set$gene_symbol)

# de_res[de_res$gene%in% c('Cpt1a','Hmgcs2','Acsm5'),]
# 均在DCM上调
# intersect(c('Cpt1a','Hmgcs2','Acsm5'),fatty_acid_genes)
# Cpt1a

# 热图
library(pheatmap)
data__scale=apply(rel_data,1,scale)
rownames(data__scale)=colnames(rel_data)
data__scale=t(data__scale)

plot_genes = de_res[de_res$gene%in% unique(c(fatty_acid_genes,'Cpt1a','Hmgcs2','Acsm5')),'ids']
data_4_1 = data__scale[plot_genes,]
plot_genes_names = de_res[match(rownames(data_4_1),de_res$ids),'gene']
rownames(data_4_1)=plot_genes_names
data_4_1 = data_4_1[rownames(data_4_1) %in%de_res[de_res$regulate!='nosig','gene'] ,]

anno_colors=list(
  group=color_samp
)
mat_dist = dist(data_4_1)
hclust_1 = hclust(mat_dist)
test=rev(hclust_1$labels[hclust_1$order])

data_4_1 = data_4_1[,c(paste0('DCM_',1:3),paste0('Control_',1:3))]
anno_col=data.frame(
  row.names = colnames(data_4_1),
  group=rev(group)
)


pdf('04.protein_dia/fig_4_1_heatmap_sig_fatty_betta.pdf',width=4,height=4)
pheatmap(data_4_1[test,],scale='none',
         show_rownames = T,cluster_rows = T,
         show_colnames = F,
         cluster_cols = F,
         cellwidth = 15,
         treeheight_row = 10,
         annotation_col = anno_col,annotation_colors = anno_colors,
         annotation_names_col = F,
         color=colorRampPalette(c('#191970','white','#FF4500'))(100),
         border_color = NA
)
dev.off()
#-------------- 热图 所有脂肪酸beita氧化相关基因。
data_4_1.1 = data__scale[plot_genes,]
plot_genes_names = de_res[match(rownames(data_4_1.1),de_res$ids),'gene']
rownames(data_4_1.1)=plot_genes_names

pdf('04.protein_dia/fig_4_1_heatmap_full_fatty_betta.pdf',width=4,height=10)
pheatmap(data_4_1.1,scale='none',
         show_rownames = T,cluster_rows = T,
         show_colnames = F,
         cluster_cols = F,
         cellwidth = 15,
         treeheight_row = 10,
         annotation_col = anno_col,annotation_colors = anno_colors,
         annotation_names_col = F,
         color=colorRampPalette(c('#191970','white','#FF4500'))(100),
         border_color = NA
)
dev.off()
#  火山图-----------------------------------------------------------
library(ggrepel)
label_genes = rownames(data_4_1)
label_genes = label_genes[label_genes%in%unique(de_res[de_res$log2FC>1.3&de_res$regulate=='up','gene'])]
label_genes=c(label_genes,'Cpt1a',de_res[de_res$log2FC< -3&de_res$regulate=='down','gene'])


aa = table(de_res$regulate)
names(aa)=c('DCM_Down','No_Sig','DCM_Up')
aa = as.data.frame.table(aa)
labeldata = de_res[de_res$gene %in% label_genes,]


pdf("04.protein_dia/4_1_volcano_fatty_beta.pdf",width = 6,height = 5)
ggplot(de_res,aes(x=log2FC,y=-log10(pvalue)))+
  geom_point(size=2,aes(color=regulate))+
  labs(x='log2(Fold Change)',y="-log10(pvalue)")+
  geom_vline(xintercept = c(log2(2/3),log2(1.5)),lty=2,col ="grey30",lwd=0.5)+
  geom_hline(yintercept=-log10(0.05),lty=2,col = "grey30",lwd=0.5)+ 
  theme_test()+theme(
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
  geom_text_repel(data=labeldata,aes(x=log2FC,y=-log10(pvalue),label=gene),size=3.5,max.overlaps = 30,show.legend = F)+
  geom_point(data=labeldata,aes(x=log2FC,y=-log10(pvalue)),shape=21,size=2.5,col='black')

dev.off()

################################################################################
# 差异基因 GO  KEGG富集
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
#------------------------------------------------------
DCM_UP_genes = unique(de_res[de_res$regulate=='up','gene'])
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
write.table(ego_ALL_up,file='04.protein_dia/4_2_enrich_GO_DCM_UP.txt',sep='\t',quote=F,row.names = F)

# KEGG 富集 up
enrich_kegg_up <- enrichKEGG(gene = DCM_UP_genes_ids$ENTREZID,
                             organism = 'mmu',
                             pvalueCutoff = 0.05,
                             pAdjustMethod = 'none',
)
enrich_kegg_up = DOSE::setReadable(enrich_kegg_up,OrgDb='org.Mm.eg.db',keyType = 'ENTREZID')
enrich_kegg_up = as.data.frame(enrich_kegg_up)
write.table(enrich_kegg_up,file='04.protein_dia/4_2_enrich_KEGG_DCM_UP.txt',sep='\t',quote=F,row.names = F)

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
plot_top_go_up = rbind(plot_top_go_up,go_fatty_relation[1:5,])
plot_top_go_up = rbind(plot_top_go_up,go_fatty_relation[go_fatty_relation$ONTOLOGY=='MF',])
plot_top_go_up = plot_top_go_up[!duplicated(plot_top_go_up$ID),]
plot_top_go_up = plot_top_go_up[!plot_top_go_up$ID%in% c('GO:0043177','GO:0031406'),]


plot_top_go_up$Description2 = capitalize(plot_top_go_up$Description)
plot_top_go_up=plot_top_go_up %>%
  arrange(ONTOLOGY,desc(Count))
plot_top_go_up$Description2=factor(plot_top_go_up$Description2,levels=rev(plot_top_go_up$Description2))


pdf('04.protein_dia/4_2_enrich_GO_DCM_UP_top5.pdf',width = 5,height = 3.5)
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
plot_up_kegg$Description2=gsub('\\s-\\sMus\\smusculus\\s\\(house\\smouse\\)','',plot_up_kegg$Description)

plot_up_kegg = plot_up_kegg[order(-log10(plot_up_kegg$pvalue),decreasing = F),]
plot_up_kegg$Description2 = factor(plot_up_kegg$Description2,levels=plot_up_kegg$Description2)
label_col = rep('black',15)
label_col[grepl('fatty',plot_up_kegg$Description2,ignore.case = T)]='red'

pdf('04.protein_dia/4_2_enrich_KEGG_DCM_UP_top15.pdf',width = 5.5,height = 3.5)
ggplot(plot_up_kegg,aes(x=-log10(pvalue),y=Description2,color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+
  xlim(sort(-log10(plot_up_kegg$pvalue))[c(1,15)]+c(-1,1))+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.4,'cm'),
        axis.text.y=element_text(size=12,color=label_col))+
  labs(x=bquote(-log[10]~(pvalue)))+
  scale_y_discrete(position = 'right')+
  scale_color_gradient(low='blue',high='red')+
  guides(color='none')
dev.off()
#########################################################################
# DCM DOWN 基因富集分析
DCM_DOWN_genes = unique(de_res[de_res$regulate=='down','gene'])
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
write.table(ego_ALL_down,file='04.protein_dia/4_2_enrich_GO_DCM_DOWN.txt',sep='\t',quote=F,row.names = F)

# KEGG 富集 down
enrich_kegg_down <- enrichKEGG(gene = DCM_DOWN_genes_ids$ENTREZID,
                             organism = 'mmu',
                             pvalueCutoff = 0.05,
                             pAdjustMethod = 'none',
)
enrich_kegg_down = DOSE::setReadable(enrich_kegg_down,OrgDb='org.Mm.eg.db',keyType = 'ENTREZID')
enrich_kegg_down = as.data.frame(enrich_kegg_down)
write.table(enrich_kegg_down,file='04.protein_dia/4_2_enrich_KEGG_DCM_DOWN.txt',sep='\t',quote=F,row.names = F)

# GO可视化 down
library(dplyr)
library(stringr)
library(Hmisc)

plot_top_go_down = ego_ALL_down %>% 
  group_by(ONTOLOGY) %>% 
  slice_head(n=5)

plot_top_go_down$Description2 = capitalize(plot_top_go_down$Description)
plot_top_go_down=plot_top_go_down %>%
  arrange(ONTOLOGY,desc(Count))
plot_top_go_down$Description2=factor(plot_top_go_down$Description2,levels=rev(plot_top_go_down$Description2))


pdf('04.protein_dia/4_2_enrich_GO_DCM_DOWN_top5.pdf',width = 5.5,height = 3.5)
ggplot(plot_top_go_down,aes(x=Count,y=Description2,color=-log10(p.adjust),shape=ONTOLOGY))+
  geom_point(size=3)+xlim(c(min(plot_top_go_down$Count)-2,max(plot_top_go_down$Count)+2))+theme_bw()+
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        legend.key.size = unit(0.5,'cm'),
        axis.text.y=element_text(size=12))+
  labs(shape='Ontology',color=bquote(-log[10]~(FDR)),x='Count')+
  scale_y_discrete(position = 'right',label=function(x){stringr::str_wrap(x,width = 60)} )+
  scale_color_gradient(low='blue',high='red')+
  guides(shape=guide_legend(order=1))
dev.off() 

# KEGG可视化 down
library(tidyverse)
library(cowplot)

plot_down_kegg = head(enrich_kegg_down,15)
plot_down_kegg$Description2=gsub('\\s-\\sMus\\smusculus\\s\\(house\\smouse\\)','',plot_down_kegg$Description)

plot_down_kegg = plot_down_kegg[order(-log10(plot_down_kegg$pvalue),decreasing = F),]
plot_down_kegg$Description2 = factor(plot_down_kegg$Description2,levels=plot_down_kegg$Description2)

pdf('04.protein_dia/4_2_enrich_KEGG_DCM_DOWN_top15.pdf',width = 5.5,height = 3.5)
ggplot(plot_down_kegg,aes(x=-log10(pvalue),y=Description2,color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+
  xlim(sort(-log10(plot_down_kegg$pvalue))[c(1,15)]+c(-1,0.5))+
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
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)


fig_data_4_3 = rel_data[match(de_res[de_res$gene%in%c('Cpt1a','Hmgcs2','Acsm5'),'ids'],rownames(rel_data)),]
rownames(fig_data_4_3) = gsub('^.*__','',rownames(fig_data_4_3))
fig_data_4_3 = cbind(gene=rownames(fig_data_4_3),fig_data_4_3)
fig_data_4_3=melt(as.data.frame(fig_data_4_3),value.name ='relexp' ,id.vars = 'gene',variable.name = 'sample')
fig_data_4_3$gene=factor(fig_data_4_3$gene,levels=c('Acsm5','Cpt1a','Hmgcs2'))
fig_data_4_3$group = gsub('_.*$','',fig_data_4_3$sample)
fig_data_4_3$group = factor(fig_data_4_3$group,levels=c('DCM','Control'))

fig_data_4_3$relexp = as.numeric(fig_data_4_3$relexp)

pdf('04.protein_dia/4_3_gene_exp_compare.pdf',width = 4,height = 3.5)
ggboxplot(fig_data_4_3,x='gene',y='relexp',color = 'group',add='jitter')+
  scale_color_manual(values=color_samp)+
  theme(axis.title.x = element_blank())+
  stat_compare_means(aes(group=group),label.y=2.5,method='anova',label='p.signif')+
  labs(y='Relative Expression')
dev.off()
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


#----GSVA分析
library(GSVA)
data_4_4=rel_data[de_res$ids,]
identical(rownames(data_4_4),de_res$ids)# TRUE
rownames(data_4_4) = make.unique(de_res$gene)

gp=gsvaParam(as.matrix(data_4_4), kbhb_set,
             maxDiff =FALSE,
             kcdf='Gaussian',
             minSize = 3)
res_4_4_gsva <- gsva(gp,BPPARAM=BiocParallel::SnowParam(workers=6))
res_4_4_gsva = as.data.frame(res_4_4_gsva)

# 绘图展示GSVA score
library(pheatmap)
pdf('04.protein_dia/4_4_kbhb_gsva_score_heatmap.pdf',width = 10,height = 8)
pheatmap(res_4_4_gsva, show_colnames = T, 
         scale = "none",angle_col = "45",
         cluster_row = T,cluster_col = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

#找差异GSVA得分条目
library(limma)
group_list=factor(group,levels=c('DCM','Control'))

design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(res_4_4_gsva)
contrast.matrix<-makeContrasts("DCM-Control",levels = design)

fit <- lmFit(res_4_4_gsva,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯
DE_term <- topTable(fit2, coef = 1,n=Inf,adjust.method = 'none')

# 输出GSVA表
res_4_4_gsva1=cbind(pathway=rownames(res_4_4_gsva),res_4_4_gsva)
res_4_4_gsva1=res_4_4_gsva1[order(res_4_4_gsva1$pathway),]
res_4_4_gsva2=cbind(res_4_4_gsva1,DE_term[match(rownames(res_4_4_gsva1),rownames(DE_term)),])
res_4_4_gsva2$B=NULL
res_4_4_gsva2$t=NULL

wb = openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,'GSVA_Kbhb')
openxlsx::writeData(wb,'GSVA_Kbhb',x=res_4_4_gsva2)
openxlsx::addFilter(wb,1,1,cols=1:ncol(res_4_4_gsva2))
nwid=rep(15,ncol(res_4_4_gsva2))
nwid[1]=50
openxlsx::setColWidths(wb,1,cols=1:ncol(res_4_4_gsva2),widths = nwid)
openxlsx::freezePane(wb,1,firstRow = T)
openxlsx::saveWorkbook(wb,'04.protein_dia/4_4_GSVA_pathway_table_protein.xlsx',overwrite = T)

# 显著的过程通路
DE_set <- DE_term[DE_term$logFC > 1&DE_term$adj.P.Val<0.05 ,]
dat <- res_4_4_gsva[match(rownames(DE_set),rownames(res_4_4_gsva)),]
library(ComplexHeatmap)
dat = dat[,c(paste0('DCM_',1:3),paste0('Control_',1:3))]

top_anno = HeatmapAnnotation(
  group=rev(group),
  annotation_legend_param = list(
    title='group'
  ),
  col=list(group=color_samp),
  annotation_name_gp = gpar(col='black')
)
ht = Heatmap(as.matrix(dat),col=colorRampPalette(c("navy", "white", "firebrick3"))(50),row_dend_width = unit(5,'mm'),top_annotation = top_anno,show_column_names = F,cluster_columns = F,heatmap_legend_param = list(title=""),row_names_max_width =unit(20,'cm'))

pdf('04.protein_dia/4_4_kbhb_gsva_score_heatmap_sig.pdf',width = 11,height = 2.5)
draw(ht,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()
#------------------------------------ GSEA
library(clusterProfiler)
library(enrichplot)
enrich_data_4_4 = de_res$log2FC
names(enrich_data_4_4) = de_res$gene
enrich_data_4_4 = enrich_data_4_4[de_res$gene!='--']
enrich_data_4_4 = sort(enrich_data_4_4,decreasing = T)
enrich_data_4_4=enrich_data_4_4[enrich_data_4_4!=0]

#GSEA 富集分析
enrich_4_4_kbhb = GSEA(geneList = enrich_data_4_4,TERM2GENE = Kbhb_relation_set,seed=123,pvalueCutoff = 0.05,minGSSize = 3,maxGSSize = 200,pAdjustMethod = 'none')
#gseaplot2(enrich_2_4_kbhb,geneSetID = enrich_2_4_kbhb@result$Description)
write.table(enrich_4_4_kbhb@result,file='04.protein_dia/4_4_GSEA_kbhb_sig.txt',sep='\t',quote=F,row.names = F)
gsea_sigs = c(
  'GOBP_FATTY_ACID_BETA_OXIDATION',
  'KEGG_FATTY_ACID_METABOLISM',
  'KEGG_CITRATE_CYCLE_TCA_CYCLE',
  'REACTOME_GLUCONEOGENESIS'
)

for(iterm in gsea_sigs){
  pdf(paste0('04.protein_dia/4_4_GSEA_sig_term_',iterm,'.pdf'),width=6,height = 5)
  gpt = gseaplot2(enrich_4_4_kbhb,geneSetID = iterm,title=iterm)
  print(gpt)
  dev.off()
}


################################################################################
save.image(paste0('04.protein_dia.',Sys.Date(),'.RData'))
################################################################################
