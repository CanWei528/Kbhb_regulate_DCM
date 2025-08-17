# 八、联合代谢组分析X蛋白-KXX位点β羟基丁酰化修饰串的下游代谢
setwd('03.results')
if(!dir.exists('08.kbhb_and_meta')){
  dir.create('08.kbhb_and_meta')
}
color_samp = c('DCM'='#D96244','Control'='#5EB5C8')
######################################################################
# 使用数据
# 代谢数据
data_2 = read.table('01.data/代谢-小鼠(MJ20240401212) - 10 - 1/interaction_results/4.ExpDiff/01.TwoGroupExpDiff/ExpDiff_s1s2s3s4s5_d1d2d3d4d5/mix/DiffTest/db_vs_dbm.diff.exp.xls',header=T,sep='\t',quote='"',comment.char = '!')

exp_data = data_2[,c(1,7:16)]
exp_data = exp_data[!grepl('^metab_',exp_data$Metabolite),]
rownames(exp_data)=exp_data$Metabolite
exp_data$Metabolite=NULL
colnames(exp_data) = gsub('dbm_s','Control_',colnames(exp_data))
colnames(exp_data) = gsub('db_d','DCM_',colnames(exp_data))
exp_data = exp_data[,!grepl('4|5',colnames(exp_data))]

# 读取kbhb修饰数据
data_total = openxlsx::read.xlsx('01.data/kbhb报告/4-Differentially_expressed_protein/T-test_analysis/Differentially_annot.xlsx',sheet = 'dbvsdbm')
norm_exp = data_total[,c(1:3,8,17:22)]
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

# 处理相对表达值
norm_exp$kbhb_site = paste0(norm_exp$Protein.accession,'_',norm_exp$Amino.acid,norm_exp$Position,'_',norm_exp$Gene.name)
norm_exp[,1:4]=NULL
rownames(norm_exp)=norm_exp$kbhb_site
norm_exp$kbhb_site=NULL
colnames(norm_exp) = gsub('dbm','Control_',colnames(norm_exp))
colnames(norm_exp) = gsub('db','DCM_',colnames(norm_exp))
group=c(rep('Control',3),rep('DCM',3))

key_genes = 'Atp5f1a'
key_sites_de = de_res[de_res$Gene.name==key_genes,]
key_sites_de$sites2 = paste0(key_sites_de$Amino.acid,key_sites_de$Position)
key_sites_de = key_sites_de[!is.na(key_sites_de$regulate),]
key_sites_de = key_sites_de[key_sites_de$log2FC>0.2,]

key_sites_de = key_sites_de[,c(6,8,11)]
#####################################################################
# 1、分析X蛋白不同β羟基丁酰化修饰位点丁酰化水平与代谢产物（3-Hydroxybutanoic Acid（β-羟基丁酸（β-OHB））；Isocitric Acid）的相关性（图A）
library(ggplot2)
library(RColorBrewer)
data_8_1 = rbind(
  exp_data[c('3-Hydroxybutanoic Acid','Isocitric Acid'),],
  norm_exp[key_sites_de$kbhb_site,]
)
data_8_1 = cbind(name=rownames(data_8_1),data_8_1)
data_8_1 = reshape2::melt(data_8_1,id.vars='name',value.name = 'exp',variable.name='sample')

data_8_1$type=ifelse(
  data_8_1$name%in%c('3-Hydroxybutanoic Acid','Isocitric Acid'),'meta','kbhb_sites'
)
# y_data 修饰表达数据
y_data = data_8_1[data_8_1$type!='meta',]
y_data$sites = key_sites_de[match(y_data$name,key_sites_de$kbhb_site),'sites2']
y_data = y_data[,c(2,3,5)]
sample_order=c(paste0('Control_',1:3),paste0('DCM_',1:3))


plot_list_8_1 = list()
for(cur_meta in c('3-Hydroxybutanoic Acid','Isocitric Acid')){
  #cur_meta = 'Isocitric Acid'
  # x_data 代谢表达数据
  x_data = data_8_1[data_8_1$name==cur_meta,]
  x_data = x_data[match(sample_order,x_data$sample),]
  plot_data_m = NULL
  stat_m = NULL
  for(site in unique(y_data$sites)){
    tmp_1 = y_data[y_data$sites==site,]
    tmp_1 = tmp_1[match(sample_order,tmp_1$sample),]
    plot_data_m = rbind(plot_data_m,
                        data.frame(sample=sample_order,
                                   xpos=x_data$exp,
                                   ypos=tmp_1$exp,
                                   type=site)
                        )
    corr = cor.test(x_data$exp,tmp_1$exp)
    p_value = signif(corr$p.value,3)
    cor_value= signif(corr$estimate[['cor']],3)
    stat_m = rbind(stat_m,
                   data.frame(site=site,cor=cor_value,pvalue=p_value))
    
  }
  plot_data_m$group = gsub('_.*$','',plot_data_m$sample)
  plot_data_m$type = factor(plot_data_m$type,levels=sort(unique(plot_data_m$type)))
  stat_m = stat_m[match(levels(plot_data_m$type),stat_m$site),]
  labs_corr = paste0(stat_m$site,', R=',stat_m$cor)
  plot_data_m$group = factor(plot_data_m$group,levels=c('DCM','Control'))
  
  sites_col = brewer.pal(6,'Dark2')
  names(sites_col) = levels(plot_data_m$type)
  
  # x轴 代谢物表达，y轴kbhb修饰位点水平
  pp = ggplot(plot_data_m,aes(x=xpos,y=ypos))+
    geom_line(aes(color=type),linewidth=1)+
    geom_point(aes(shape=group),size=2,stroke=1)+
    theme_test()+
    labs(x=paste0('Expression of ',cur_meta),y='Expression of Kbhb on Atp5f1a')+
    scale_shape_manual(values=c('DCM'=20,'Control'=1))+
    scale_color_manual(values=sites_col,breaks = levels(plot_data_m$type),labels=labs_corr)+
    theme(
      legend.title = element_blank(),
      axis.text = element_text(size=10),
      axis.title = element_text(size=12)
    )+
    guides(color=guide_legend(order=1))
    #+geom_point(aes(color=type,shape=group),size=2,stroke=1,show.legend = F)
    
  plot_list_8_1[[cur_meta]] = pp
}

pdf('08.kbhb_and_meta/8_1_key_sites_and_meta_corr_line_black_point.pdf',width=6,height = 6)
cowplot::plot_grid(plotlist = plot_list_8_1,ncol=1,align='hv')
dev.off()
#####################################################################
# 2、分析X蛋白不同β羟基丁酰化修饰位点丁酰化水平与其它代谢产物的相关性，筛选可能串扰的下游代谢产物（图E）
# ciber 列是基因集 行是样本 改成 行是代谢物 列是样本
# exp_use_log 行是基因，列是样本，改成，行是位点，列是样本
# 类推
library(Hmisc)
library(dplyr)
# exp_data # 代谢物表达
# norm_exp[key_sites_de$kbhb_site,] # 位点修饰水平

nc = t(rbind(exp_data,norm_exp[key_sites_de$kbhb_site,]))
m = rcorr(nc)$r[1:nrow(exp_data),(ncol(nc)-length(key_sites_de$kbhb_site)+1):ncol(nc)]

##计算p值
p = rcorr(nc)$P[1:nrow(exp_data),(ncol(nc)-length(key_sites_de$kbhb_site)+1):ncol(nc)]

tmp <- matrix(case_when(as.vector(p) < 0.001 ~ "***",
                        as.vector(p) < 0.01 ~ "**",
                        as.vector(p) < 0.05 ~ "*",
                        TRUE ~ ""), nrow = nrow(p))
colnames(m) = key_sites_de[match(colnames(m),key_sites_de$kbhb_site),'sites2']
#------------------------------------------------------------------------------
# 找sub 部分绘图 m  tmp
pos_inx = rowSums(m>0.7)>2
#neg_inx = rowSums(m< -0.7)>3
index = pos_inx
m.sub = m[index,]
tmp.sub = tmp[index,]

tmp.sub = tmp.sub[nchar(rownames(m.sub))< 48,]
m.sub = m.sub[nchar(rownames(m.sub))< 48,]

dim(tmp.sub)
##绘制热图
pdf('08.kbhb_and_meta/8_2_kbhb_key_sites_corr_meta.pdf',width=7,height = 8)
pheatmap::pheatmap(m.sub,
                   border_color = '#F0F8FF',
                   display_numbers =tmp.sub,
                   fontsize_number = 10,
                   angle_col = 0,
                   color = colorRampPalette(c("#4682B4", "white", "#d71e22"))(50),
                   treeheight_col = 0,
                   cellwidth = 30,
                   #cellheight = 10,
                   treeheight_row = 15)
dev.off()
#"#8B0000"
#####################################################################
# 4、基于1、2结果筛选X蛋白最佳β羟基丁酰化修饰位点
# Atp5f1a_K239 最佳
# Atp5f1a_K161和Atp5f1a_K126次选
#####################################################################
save.image('08.kbhb_and_meta.0308.RData')
#######################################################################