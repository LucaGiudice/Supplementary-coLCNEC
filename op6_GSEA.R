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
library("SPIA")
library("gage")

#Set variables ----
#Input
out_op4_path="output/op4_data_prep.rda"
out_op5_path="output/op5_DE_results.rda"
gsea_pathway_path="data/msigdb_pathways.rds"
#Output
out_op6_path="output/op6_GS_analysis.rda"

#Load data ----
load(out_op4_path)
load(out_op5_path)
gsea_pathways=pathways=readRDS(gsea_pathway_path)

#Data analysis ----
spia_col=c("Name","NDE","tA","pGFdr","Status","KEGGLINK")
GS_res_l=list()
res0=matrix(rep("no significative data",10),2,5);colnames(res0)=seq(1,5);rownames(res0)=seq(1,2)

for(k in 1:length(dataXcontrast_l)){
  contr_name=names(dataXcontrast_l)[k]
  info_df=dataXcontrast_l[[k]]$infoc
  lcpm=dataXcontrast_l[[k]]$lcpm
  tT_ov=de_contrast_l[[k]]$tT_ov
  
  #Convert the ensg rownames in lcpm with gene names ----
  #Match the order with rawdata and lcpm
  tT_ov=tT_ov[rownames(lcpm),]

  #GSEA ----
  #pathways=pathways[1:50]
  GSEA_res=gage(exprs = lcpm, gsets = gsea_pathways, 
                ref = which(info_df$group==0), 
                samp =  which(info_df$group==1), 
                compare='as.group', 
                same.dir=T)
  GSEA_res_g=GSEA_res$greater
  GSEA_res_s=GSEA_res$less
  GSEA_res=rbind(GSEA_res_g,GSEA_res_s)
  
  #Get the lFC of the de entrez genes and the source of experiment ----
  #Remove the rawdata having duplicated entrez genes
  entr_tT_ov=tT_ov[!duplicated(tT_ov$entrez),]
  entr_tT_ov=na.omit(entr_tT_ov)
  #Get the DE rawdata
  if(sum(entr_tT_ov$adj.P.Val<=0.05)<1){
    GS_res_l[[k]]=list(spia_df=res0,GSEA_df=GSEA_res,de_enr=res0);
    next;
  }
  de_tT=entr_tT_ov[entr_tT_ov$adj.P.Val<=0.05,]
  #Get the lFC of the de entrez genes
  de_entr2FC=de_tT$logFC
  names(de_entr2FC)=de_tT$entrez
  #Get the entrez source of experiment
  all_entr=entr_tT_ov$entrez
  spia_res=spia(de=de_entr2FC,all=all_entr,organism="hsa",nB=200,plots=F,
                verbose=T,beta=NULL,combine="fisher",data.dir="data/SPIA_db/")
  #Keep useful information and order the pathways
  spia_res=spia_res[,match(spia_col,colnames(spia_res))]
  spia_res=spia_res[order(spia_res$pGFdr,decreasing = F),]

  #Reconvert the entrez of the genes perturbating pathways to symbol
  gsXpath=substr(spia_res[,6],49,200)
  gsXpath_l=strsplit(gsXpath, "+", fixed=T)
  gsXpath_l=lapply(gsXpath_l,function(x){
    gsSYMB=as.data.frame(rownames(de_tT)[match(x,de_tT$entrez)])
    gsSYMB=as.data.frame(t(gsSYMB))
    rownames(gsSYMB)="gsSYMB"
    return(gsSYMB)
  })
  gsXpath_df=rbindlist(gsXpath_l,fill = T)
  #Update the spia result
  spia_res2=cbind(spia_res,gsXpath_df)
  
  #Enrichment analysis of DE genes -----
  prop_de=0
  de_pathways=pathways
  tT=de_contrast_l[[k]]$tT
  genes=rownames(tT)
  for(i in 1:length(pathways)){
    gene_set=pathways[[i]];n_set=length(gene_set)
    de_gene_set=intersect(gene_set,genes);n_de=length(de_gene_set)
    
    if(n_de!=0){
      prop=n_de/n_set
      prop_de=c(prop_de,prop)
      de_gene_set_v=tT$logFC[match(de_gene_set,genes)]
      names(de_gene_set_v)=de_gene_set
      de_pathways[[i]]=de_gene_set_v
    }else{
      prop_de=c(prop_de,0)
      de_pathways[[i]]=de_gene_set
    }
  }
  prop_de=prop_de[-1]
  keep_pathways=prop_de>=0.5
  de_pathways=de_pathways[keep_pathways]
  prop_de=prop_de[keep_pathways]
  
  keep_pathways=sapply(de_pathways,length)!=1
  de_pathways=de_pathways[keep_pathways]
  prop_de=prop_de[keep_pathways]
  
  mean_lFCxPATH=sapply(de_pathways,function(x){mean(abs(x))})
  pathway_data=cbind(mean_lFCxPATH, prop_de)
  
  if(nrow(pathway_data)==0){pathway_data=res0}
  
  #Build the enrichemnt object list containing everything ----
  GS_res_l[[k]]=list(spia_df=spia_res2,GSEA_df=GSEA_res,de_enr=pathway_data)
}
save(GS_res_l,file=out_op6_path)
  
  
  
  