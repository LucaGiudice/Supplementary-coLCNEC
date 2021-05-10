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
library("GSVA")

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
out_op4x_path="output/op4x_artef_analysis.rda"
data4heat_path="output/op4/ALL_data_heatmap.rda"
gsea_pathway_path="data/msigdb_pathways.rds"
fast=T;pathway_type="msigdb_h";
#output
out_op13_path_png="output/op9/op9_pathways_overex_%s.png"
out_op13_path_table="output/op9/op9_pathways_overex_%s.xlsx"
out_op13_path_lcpm="output/op9/op9_pathways_overex_%s.txt"
out_op13_path_rda="output/op9_pathways_overex_%s.rda"

#Load data ----
load(out_op4_path)
rm(mat,res,rl)
load(out_op4x_path)
pathways=readRDS(gsea_pathway_path)

lcpm=lcpm[-match(artefatti_gs,rownames(lcpm)),]

#Writing report ----
groups=unique(infoc$cluster)
check_elems(infoc$case,colnames(lcpm))

keep=grep(pathway_type,names(pathways))
pathways=pathways[keep]

#Compute GSVA
if(!fast){
  gsva_res=gsva(lcpm,gset.idx.list=pathways,method="ssgsea",abs.ranking=TRUE,min.sz=1,
                max.sz=999,parallel.sz=5)
  
  #Format to get GSVA for cluster
  tmp=t(aggregate(t(gsva_res), list(infoc$cluster), mean))
  colnames(tmp)=paste("CL",tmp[1,],sep="")
  gsva_res=tmp[-1,]
  save(gsva_res,file="old/ssgsea_h.rda")
}else{
  load("old/ssgsea_h.rda")
}

#Select a cluster and the other groups
out_op13_path_rda=sprintf(out_op13_path_rda,pathway_type)
out_op13_path_lcpm=sprintf(out_op13_path_lcpm,pathway_type)
out_op13_path_png=sprintf(out_op13_path_png,pathway_type)
gsva_res_spe=gsva_res

#Row normalize
gsva_res_sc=apply(gsva_res_spe,1,function(x){
  y=rescale(x,to=c(-1,1))
  return(y)
})
gsva_res_sc=t(gsva_res_sc)

#Select a cluster and the other groups
best10=c("start")
opp_grs_l=list()
opp_grs_l[[1]]=c("CL3","CL6")
opp_grs_l[[2]]=c("CL4","CL8")
opp_grs_l[[3]]=c("CL1","CL5")
opp_grs_l[[4]]=c("CL2","CL8")
opp_grs_l[[5]]=c("CL4","CL9")
opp_grs_l[[6]]=c("CL4","CL1")
opp_grs_l[[7]]=c("CL1","CL2")
opp_grs_l[[8]]=c("CL3","CL6")
opp_grs_l[[9]]=c("CL1","CL5")
opp_grs_l[[10]]=c("CL4","CL7")

for(in_cl in 1:ncol(gsva_res_sc)){
  if(in_cl==6){next;}
  cat(in_cl,"on 10 \n")
  in_cl_name=colnames(gsva_res_sc)[in_cl]
  out_cl=opp_grs_l[[in_cl]]
  
  #Find the in and out scores
  in_gs=gsva_res_spe[,in_cl]
  out_gs=gsva_res_spe[,out_cl]
  
  #Compute the difference of the cluster enrichment score vs others
  diff_g=0
  for(g in 1:length(in_gs)){
    #diff_g=c(diff_g,dist(c(in_gs[g],mean(out_gs[g,])))[1])
    vs=c(in_gs[g],out_gs[g,])
    names(vs)[1]=in_cl_name
    dists=as.matrix(dist(vs))
    diff_g=c(diff_g,mean(as.matrix(dist(vs))[-1,1]))
  }
  diff_g=diff_g[-1]
  
  mat=cbind(rownames(gsva_res_sc),diff_g)
  mat=as.data.frame(mat,stringsAsFactors = F)
  mat$diff_g=as.numeric(mat$diff_g)
  colnames(mat)=c("pathway_enrich","diff")
  mat=mat[with(mat, order(-diff)), ]
  best10=c(best10,mat$pathway_enrich[1:20])
  
  if(in_cl==colnames(gsva_res_sc)[1]){
    app=F
  }else{
    app=T
  }
  #write.xlsx(mat,file = out_op5_path_table,sheetName = in_cl_name,row.names = F,col.names = T,
  #           append=app)
}
best10=best10[-1]
best10=unique(best10)
best10=best10[!is.na(best10)]

#Subset with the best pathways and row normalize
gsva_res_sc2plot=gsva_res_sc[best10,]

png(out_op13_path_png, units="in", width=22, height=12, res=600)
mat3=gsva_res_sc2plot
Heatmap(
  mat3,
  name="ssGSEA pathways",
  show_column_names = F,
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(colnames(mat3), rot = 60, offset = unit(1, "npc"), just = "right"),
    annotation_height = max_text_width(colnames(mat3))
  ),
  show_row_names = T,
  show_column_dend = T,
  show_row_dend = F,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),
  #right_annotation = ann_r2,
  use_raster = TRUE, raster_quality = 2,
)
dev.off()
write.table(mat3,file=out_op13_path_lcpm,quote = F,row.names = T,col.names = T)
save(gsva_res_spe,mat3,file=out_op13_path_rda)




