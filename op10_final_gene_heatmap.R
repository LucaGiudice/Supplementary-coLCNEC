#Clean workspace and memory ----
rm(list=ls())
gc()

#Set working directory----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
options(java.parameters = "-Xmx600000m")
setwd(gsub("%s","",rootDir))
set.seed(8)

#Load libraries----

library("data.table")
library("matrixStats")
library("plyr")
library("xlsx")
library("ComplexHeatmap")
library("scales")
library("Rfast")

autri = function(v){
  diff=1-((abs(v[1]-v[2]))/2)
  C <- c(1-diff,1-diff)
  B1 <- c(1+v[1],0)
  B2 <- c(1,1+v[2])
  
  return(abs((C[1]*(B1[2]-B2[2]) + B1[1]*(B2[2]-C[2]) + B2[1]*(B1[2]-C[2]))/2))
}
see = function(m,nr=5,nc=5){
  if(nrow(m)<5){
    nr=nrow(m)
  }
  if(ncol(m)<5){
    nc=ncol(m)
  }
  
  m[1:nr,1:nc]
}
check_elems = function(x,y){
  indx_ord=match(x,y)
  corr_ord=seq(1,length(y))
  all.equal(indx_ord,corr_ord)
  err=0
  for(k in corr_ord){
    check=indx_ord[k]==k
    if(!check){
      err=c(err,k)
    }
  }
  err=err[-1]
  
  if(length(err)==0){cat("Matching \n")}
  if(length(err)!=0){cat("NO Matching, returning indexes \n")}
  return(err)
}

#Set variables ----
#Input
out_op4_path="output/op4_CL1.rda"
data4heat_path="output/op4/ALL_data_heatmap.rda"
out_op4x_path="output/op4x_artef_analysis.rda"
must_gs_path="data/must_genes.txt"
#output
out_op14_path_png="output/op10/op10_final_genes.png"
out_op14_path_lcpm="output/op10/op10_final_lcpm.txt"
out_op14_path_rda="output/op10_final_genes.rda"

#Load data ----
load(out_op4_path);rm(mat,res,rl);
load(out_op4x_path)
must_gs=fread(must_gs_path,data.table = F,header = F)$V1

#Remove outliers -----
lcpm=lcpm[-match(artefatti_gs,rownames(lcpm)),]

#Get the must in genes for the heatmap selected due the analysis ----
indxs=match(must_gs,rownames(lcpm))
indxs=indxs[!is.na(indxs)]
must_lcpm=lcpm[indxs,]

#Get the best genes of the clustering for the heatmap selected due to the cl -----
load(data4heat_path)
cl_gs=setdiff(rownames(mat2),must_gs)
cl_lcpm=lcpm[match(cl_gs,rownames(lcpm)),]

#Finalize the heatmap -----
hlcpm=rbind(must_lcpm,cl_lcpm)
mat3=t(scale(t(hlcpm)))

#Get the best genes of the clustering for the heatmap selected due to the cl -----
png(out_op14_path_png, units="in", width=22, height=12, res=600)
Heatmap(
  mat3,
  top_annotation = ann_c,
  name="Consesus clustering",
  show_column_names = F,
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(colnames(mat3), rot = 60, offset = unit(1, "npc"), just = "right"),
    annotation_height = max_text_width(colnames(mat3))
  ),
  show_row_names = T,
  show_row_dend = F,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 1),
  #right_annotation = ann_r2,
  use_raster = TRUE, raster_quality = 2,
  column_split = cls_lab
)
dev.off()

write.table(hlcpm, out_op14_path_lcpm, quote = F, row.names = T, col.names = T)
save(hlcpm,file=out_op14_path_rda)

