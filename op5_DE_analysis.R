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
library("limma")
library("edgeR")
library("RColorBrewer")
source("funs/fun_processing_DE.R")

#Set variables ----
#input variable
op4_data_path="./output/op4_data_prep.rda"
hsa_db_path="data/human_db.rda"
pathways_path="data/msigdb_pathways.rds"
#output variables
op5_data_path="./output/op5_DE_results.rda"

#Load data ----
cat("Loading data \n")
load(op4_data_path)
load(hsa_db_path)
phts=readRDS(pathways_path)

#Interate over the contrasts ---
de_contrast_l=list()
for(k in 1:length(dataXcontrast_l)){
  cat(k,"on",length(dataXcontrast_l),"\n")
  #Extract the matrices related to the running contrast
  contr_name=names(dataXcontrast_l)[k]
  counts=dataXcontrast_l[[k]]$counts
  info_df=dataXcontrast_l[[k]]$infoc
  des=dataXcontrast_l[[k]]$design
  
  #DE analysis core ----
  #Define variable for the DE analysis
  group=as.factor(info_df$cls)
  samplenames=info_df$isi_names
  labels_name=info_df$isi_names
  
  #Create edgelist object and compute L and M variables
  x=DGEList(counts = counts, group = group, remove.zeros = TRUE)
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  
  #Remove genes which are low expressed and non informative of the difference between groups
  keep.exprs <- filterByExpr(x, group=group)
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  
  #Normalize data by library size
  x <- calcNormFactors(x, method = "TMM")  
  
  #DE analysis core ----
  #Select the column containing the result of the contrast of interest
  coeff_contr=ncol(des)
  #Analysis
  v <- voom(x, design=des, plot=FALSE);
  vfit <- lmFit(v, design=des)
  efit <- eBayes(vfit)
  tT_ov=topTable(efit,p.value=1,adjust.method="BH",number=Inf,coef=coeff_contr)
  #Composition of the resulting DE tables
  tT_ov$entrez=mapvalues(rownames(tT_ov), from=human_db$Symbol, to=human_db$entrez, warn_missing = FALSE)
  tT_ov$ens=mapvalues(rownames(tT_ov), from=human_db$Symbol, to=human_db$ensembl_gene_id, warn_missing = FALSE)
  tT=tT_ov[tT_ov$adj.P.Val<=0.05,]
  tT=tT[,-c(2,3,4,6)]
  tT=tT[order(tT$logFC,decreasing = T),]
  
  #Computation of camera analysis
  idx <- ids2indices(phts,id=rownames(v))
  cam.res = camera(v,idx,des,contrast=ncol(des))
  cam.res = cam.res[cam.res$FDR<=0.05,]
  
  #Compress results
  de_contrast_l[[contr_name]]=list(tT_ov=tT_ov,tT=tT,camera_res=cam.res)
  
}

save(de_contrast_l,file=op5_data_path)
