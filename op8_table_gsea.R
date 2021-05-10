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

autri = function(v){
  diff=1-((abs(v[1]-v[2]))/2)
  C <- c(1-diff,1-diff)
  B1 <- c(1+v[1],0)
  B2 <- c(1,1+v[2])
  
  return(abs((C[1]*(B1[2]-B2[2]) + B1[1]*(B2[2]-C[2]) + B2[1]*(B1[2]-C[2]))/2))
}

#Set variables ----
#Input
out_op4_path="output/op4_data_prep.rda"
out_op5_path="output/op5_DE_results.rda"
out_op6_path="output/op6_GS_analysis.rda"
#Output
out_path5="output/resulting_excels/AAA_Contrasts.xlsx"
out_path7a="output/resulting_excels/op8_GS_int_SPIA_results.xlsx"
out_path7b="output/resulting_excels/op8_GS_int_GSEA_results.xlsx"
out_path7c="output/resulting_excels/op8_GS_int_camera_results.xlsx"

#Load data ----
load(out_op4_path)
load(out_op5_path)
rm(lcpm,mat,res,rl)

#Writing report ----
labels=names(de_contrast_l)
groups=paste("CL",unique(infoc$cluster),sep="")

#Load data ----
load(out_op6_path)

#Writing report ----
labels=names(de_contrast_l)

for(k_cl in 1:length(groups)){
  cat(k_cl,"/",length(groups),"\n")
  group=groups[k_cl]
  indxs_in=grep(group,labels)
  indxs_out=setdiff(seq(1,length(labels)),indxs_in)
  de_in_l=GS_res_l[indxs_in]
  de_out_l=GS_res_l[indxs_out]
  
  de_gs_in=sapply(de_in_l,function(x){
    spia_df=x$spia_df[x$spia_df$pGFdr<=0.05,]
    gs=spia_df$Name
    return(gs)
  })
  de_gs_in=unlist(de_gs_in)
  de_gs_in=as.data.frame(table(de_gs_in))
  de_gs_in$de_gs_in=as.character(de_gs_in$de_gs_in)
  de_gs_in=de_gs_in[order(de_gs_in$Freq,decreasing = T),]
  de_gs_in$Freq=(de_gs_in$Freq*100)/length(de_in_l)
  colnames(de_gs_in)[1]="gs_spia"
  
  de_gs_out=sapply(de_out_l,function(x){
    spia_df=x$spia_df[x$spia_df$pGFdr<=0.05,]
    gs=spia_df$Name
    return(gs)
  })
  de_gs_out=unlist(de_gs_out)
  de_gs_out=as.data.frame(table(de_gs_out))
  de_gs_out=de_gs_out[order(de_gs_out$Freq,decreasing = T),]
  de_gs_out$de_gs_out=as.character(de_gs_out$de_gs_out)
  de_gs_out$Freq=(de_gs_out$Freq*100)/length(de_out_l)
  colnames(de_gs_out)[1]="gs_spia"
  gs_df=merge(x=de_gs_in,y=de_gs_out,by="gs_spia",all.x = T, all.y = F)
  gs_df$Freq.y[is.na(gs_df$Freq.y)]=0
  
  scores_m=gs_df[,c(2,3)]
  scores_m[,1]=scales::rescale(scores_m[,1],to=c(0,1))
  scores_m[,2]=scales::rescale(scores_m[,2],to=c(1,0))
  scores_m=as.matrix(scores_m)
  score_fin=apply(scores_m,1,autri)
  gs_df$score_sum=score_fin
  gs_df=gs_df[order(gs_df$score_sum,decreasing = T),]
  colnames(gs_df)=c("gs_spia","score_in","score_out","score_sum")
  if(k_cl==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(gs_df, out_path7a, sheetName=group, append=app, row.names = F, col.names = T)
  gc()
}

for(k_cl in 1:length(groups)){
  cat(k_cl,"/",length(groups),"\n")
  group=groups[k_cl]
  indxs_in=grep(group,labels)
  indxs_out=setdiff(seq(1,length(labels)),indxs_in)
  de_in_l=GS_res_l[indxs_in]
  de_out_l=GS_res_l[indxs_out]
  
  de_gs_in=sapply(de_in_l,function(x){
    gsea_df=x$GSEA_df[x$GSEA_df[,4]<=0.05,]
    gs=rownames(gsea_df)
    return(gs)
  })
  de_gs_in=unlist(de_gs_in)
  de_gs_in=as.data.frame(table(de_gs_in))
  de_gs_in$de_gs_in=as.character(de_gs_in$de_gs_in)
  de_gs_in=de_gs_in[order(de_gs_in$Freq,decreasing = T),]
  de_gs_in$Freq=(de_gs_in$Freq*100)/length(de_in_l)
  colnames(de_gs_in)[1]="gs_GSEA"
  
  de_gs_out=sapply(de_out_l,function(x){
    gsea_df=x$GSEA_df[x$GSEA_df[,4]<=0.05,]
    gs=rownames(gsea_df)
    return(gs)
  })
  de_gs_out=unlist(de_gs_out)
  de_gs_out=as.data.frame(table(de_gs_out))
  de_gs_out=de_gs_out[order(de_gs_out$Freq,decreasing = T),]
  de_gs_out$de_gs_out=as.character(de_gs_out$de_gs_out)
  de_gs_out$Freq=(de_gs_out$Freq*100)/length(de_out_l)
  colnames(de_gs_out)[1]="gs_GSEA"
  gs_df=merge(x=de_gs_in,y=de_gs_out,by="gs_GSEA",all.x = T, all.y = F)
  gs_df$Freq.y[is.na(gs_df$Freq.y)]=0
  
  scores_m=gs_df[,c(2,3)]
  scores_m[,1]=scales::rescale(scores_m[,1],to=c(0,1))
  scores_m[,2]=scales::rescale(scores_m[,2],to=c(1,0))
  scores_m=as.matrix(scores_m)
  score_fin=apply(scores_m,1,autri)
  gs_df$score_sum=score_fin
  gs_df=gs_df[order(gs_df$score_sum,decreasing = T),]
  colnames(gs_df)=c("gs_GSEA","score_in","score_out","score_sum")
  if(k_cl==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(gs_df, out_path7b, sheetName=group, append=app, row.names = F, col.names = T)
  gc()
}