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
options(java.parameters = "-Xmx60000m")
library("xlsx")

#Set variables ----
#Input
out_op1_path="output/op4_data_prep.rda"
out_op2_path="output/op5_DE_results.rda"
out_op3_path="output/op6_GS_analysis.rda"
#Output
out_path1a="output/resulting_excels/AAA_Contrasts.xlsx"
out_path1b="output/resulting_excels/AAA_INFO.xlsx"
out_path2="output/resulting_excels/op4_DE_genes_results.xlsx"
out_path31="output/resulting_excels/op5_GS_SPIA_results.xlsx"
out_path32="output/resulting_excels/op5_GS_GSEA_results.xlsx"
out_path33="output/resulting_excels/op4_DE_gene_enrichment_results.xlsx"
out_path34="output/resulting_excels/op5_GS_camera_results.xlsx"

#Load data ----
load(out_op1_path)

#Writing report ----
write.xlsx(contrasts, out_path1a, sheetName="contrasts", 
           append=F, row.names = F, col.names = T)
write.xlsx(infoc, out_path1b, sheetName="info", 
           append=F, row.names = F, col.names = T)

labels=names(dataXcontrast_l)

#Load data ----
load(out_op2_path)

#Writing report ----
for(k in 1:length(de_contrast_l)){
  cat("DE \n")
  cat(k,"\n")
  df=de_contrast_l[[k]]$tT
  label=labels[k]
  if(nrow(df)==0){df[1,]=rep("NA",ncol(df))}
  if(k==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(df, out_path2, sheetName=label, append=app, row.names = T, col.names = T)
  gc()
}

for(k in 1:length(de_contrast_l)){
  cat("Camera \n")
  cat(k,"\n")
  df=de_contrast_l[[k]]$camera_res
  label=labels[k]
  if(nrow(df)==0){df[1,]=rep("NA",ncol(df))}
  if(k==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(df, out_path34, sheetName=label, append=app, row.names = T, col.names = T)
  gc()
}

#Load data ----
load(out_op3_path)

#Writing report ----
for(k in 1:length(de_contrast_l)){
  cat("SPIA \n")
  cat(k,"\n")
  df=GS_res_l[[k]]$spia_df
  df=df[df[,1]<=0.1,]
  label=labels[k]
  if(nrow(df)==0){df[1,]=rep("NA",ncol(df))}
  if(k==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(df, out_path31, sheetName=label, append=app, row.names = T, col.names = T)
  gc()
}

for(k in 1:length(GS_res_l)){
  cat(k,"/","7","\n")
  df=as.data.frame(GS_res_l[[k]]$GSEA_df)
  df=na.omit(df)
  if(sum(df[,4]<=0.05)<2){next;}
  df=df[df[,4]<=0.05,]
  df1=df[order(df$set.size,decreasing = F),]
  df1=df1[1:50,]
  df2=df[order(abs(df$stat.mean),decreasing = T),]
  df2=df2[1:50,]
  df3=rbind(df1,df2)
  df3=na.omit(df3)
  df3=unique(df3)
  if(nrow(df)==0){next;}
  label=labels[k]
  if(nrow(df)==0){next;}
  if(k==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(df3, out_path32, sheetName=label, append=app, row.names = T, col.names = T)
  gc()
}

for(k in 1:length(de_contrast_l)){
  cat("DE enr \n")
  cat(k,"\n")
  df=GS_res_l[[k]]$de_enr
  df=df[df[,2]>=0.80,]
  df=df[order(df[,1],decreasing = T),]
  df=df[(df[,1]>=1.5),]
  label=labels[k]
  if(nrow(df)==0){df=rbind(df,rep("NA",ncol(df)))}
  if(k==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(df, out_path33, sheetName=label, append=app, row.names = T, col.names = T)
  gc()
}