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
#output
out_op13_path_png="output/op9/op9_overexp_excl2_lcpm.png"
out_op13_path_table="output/op9/op9_overexp_excl2_table.xlsx"
out_op13_path_lcpm="output/op9/op9_overexp_excl2_lcpm.txt"
out_op13_path_rda="output/op9_overexp_excl2_analysis.rda"

#Load data ----
load(out_op4_path)
rm(mat,res,rl)
load(out_op4x_path)

#Remove outlier genes
lcpm=lcpm[-match(artefatti_gs,rownames(lcpm)),]

#Writing report ----
groups=unique(infoc$cluster)
check_elems(infoc$case,colnames(lcpm))
min=0.3;max=0.8

#Local
for(k_cl_in in 1:length(groups)){
  cat(k_cl_in,"/",length(groups),"\n")
  
  group_in=groups[k_cl_in]
  indxs_in=which(infoc$cluster==group_in)
  lcpm_in=lcpm[,indxs_in]
  qin=rowQuantiles(lcpm_in,probs = min)
  
  groups_out=setdiff(groups,group_in)
  
  for(k_cl_out in 1:length(groups_out)){
    
    group_out=groups_out[k_cl_out]
    indxs_out=which(infoc$cluster==group_out)
    lcpm_out=lcpm[,indxs_out]
    qout=rowQuantiles(lcpm_out,probs = max)
    
    bol_sat=qin>qout
    if(k_cl_out==1){
      bol_df=bol_sat
    }else{
      bol_df=cbind(bol_df,bol_sat)
    }
  }
  
  sum_score=rowSums(bol_df)
  if(k_cl_in==1){
    score_df=sum_score
  }else{
    score_df=cbind(score_df,sum_score)
  }
}
colnames(score_df)=paste("CL",groups,sep="")

r_min_scores=apply(score_df,1,function(x){
  y=x[order(x,decreasing = T)]
  z=min(y[1:2])
  return(z)
})

  print(quantile(r_min_scores,probs=c(0.6,0.7,0.8,0.9,1)))
th=8
bol_keep=r_min_scores>=th
  print(sum(bol_keep))

score_df1=score_df[bol_keep,]
for(k_g in 1:nrow(score_df1)){
  g=score_df1[k_g,]
  grs=colnames(score_df1)[which(g>=th)]
  stats_g=c(rownames(score_df1)[k_g],grs)
  if(k_g==1){
    stats_gs=stats_g
  }else{
    stats_gs=rbind(stats_gs,stats_g)
  }
}
colnames(stats_gs)=c("gene","excl1","excl2");rownames(stats_gs)=seq(1,nrow(stats_gs))
  
lcpm1=lcpm[bol_keep,]

#Heatmap
load(data4heat_path)

#Scale the matrix
mat3=t(scale(t(lcpm1)))

png(out_op13_path_png, units="in", width=22, height=12, res=600)
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

write.xlsx(stats_gs,file=out_op13_path_table,sheetName = "stats_genes",row.names = F)
write.table(lcpm1, out_op13_path_lcpm, quote = F, row.names = T, col.names = T)
save(lcpm1,stats_gs,file=out_op13_path_rda)

