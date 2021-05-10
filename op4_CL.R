#Clean workspace and memory ----
rm(list=ls())
gc()

#Set working directory----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
setwd(gsub("%s","",rootDir))
set.seed(8)

#Load libraries----
library("data.table")
library("matrixStats")
library("plyr")
library("rstatix")
library("MASS")
library("HiDimDA")
library("cola")
source("funs/fun_processing_LDA.R")
source("funs/fun_cl.R")
library("ComplexHeatmap")
library("xlsx")

#Set variables ----
#Input
op2_data_path="./data/op1_data_rdy.rda"
op3_data_path1="./output/op3_dga_analysis1.rda"
#Output
op4_output_path="./output/op4_CL1.rda"
op4_output_path2="./output/op4/ALL_data_heatmap.rda"
op4_mat_path="./output/resulting_excels//lpcm.txt"
op4_info_path="./output/resulting_excels//info.xlsx"
#Set variables
server=T;ncores=16

#Load data ----
load(op2_data_path)
load(op3_data_path1)

#High dimensional LDA ----
write.table(lcpm, op4_mat_path, quote = F, row.names = T, col.names = T)

#Extract data and compose the annotation df for the heatmap ----
comb_cl=paste(infoc$source,infoc$isto,sep="_")
infoc$comb_cl=comb_cl
group=as.factor(comb_cl)

x=unique(infoc$col_source);names(x)=unique(infoc$source)
isto_df=unique(infoc[,c(2,8)])
y=isto_df[,2];names(y)=isto_df[,1]
ann_df=infoc[,c(1,2)];ann_df$source=as.factor(ann_df$source)
ref=as.numeric(as.factor(paste(ann_df$source,ann_df$isto,sep="_")))
col_list=list(source=x,isto=y)

#Consensus clustering ----
min_intr2=min(setting_df$intr2)
setting_df=setting_df[min_intr2==setting_df$intr2,]
best_settings=as.numeric(rownames(setting_df))
top_genes2=unique(unlist(top_genes2_l[best_settings]))
rl=rl2_l[[best_settings[1]]]

#Check 
best_meth=get_best_meth(rl,as.character(group),10)

#Get the params of the best cons. clustering 
top_dist=best_meth$top_dist
top_cl=best_meth$top_cl
top_k=best_meth$top_k
res = rl[top_dist, top_cl]
#Prepare the path for the plot
plot_cl=paste("output/op4/op4_cl_",top_dist,"_",top_cl,"_k",top_k,".png",sep="")

#Produce the heatmap ----
#Extract the cluster membership labels
cls_lab=get_classes(res, k = top_k)$class
infoc$cluster=cls_lab
write.xlsx(infoc, op4_info_path, sheetName="info", append=F, row.names = F, col.names = T)

#Extract col annotation
ann_c = HeatmapAnnotation(df = ann_df,
                          col = col_list)

#Extract top genes
mat=res@.env[["data"]]

x=get_signatures(res, k = top_k, plot = F)
x=x[order(x$fdr,decreasing = F),]
top_genes=x$which_row

mat2=mat[top_genes,]

#Scale the matrix
mat3=t(scale(t(mat2)))

png(plot_cl, units="in", width=22, height=12, res=300)
Heatmap(
  mat3,
  top_annotation = ann_c,
  name="Consesus clustering",
  show_column_names = F,
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(colnames(mat3), rot = 60, offset = unit(1, "npc"), just = "right"),
    annotation_height = max_text_width(colnames(mat3))
  ),
  show_row_names = F,
  show_row_dend = F,
  cluster_columns = TRUE,
  #row_names_gp = gpar(fontsize = 12),
  #right_annotation = ann_r2,
  use_raster = TRUE, raster_quality = 2,
  column_split = cls_lab
)
dev.off()
save(rl,res,top_k,mat,infoc,lcpm,file=op4_output_path)
save(mat3,mat2,ann_c,cls_lab,file=op4_output_path2)
