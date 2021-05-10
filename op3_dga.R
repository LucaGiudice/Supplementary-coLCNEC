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
library("cola")
library("xlsx")
source("funs/fun_cl.R")
options(java.parameters = "-Xmx6000000m")

#Set variables ----
#Input
op1_data_path="data//op1_data_rdy.rda"
gs2keep_path="data/genes2keep.txt"
#output
op3_dga_path1="output/op3_dga_analysis1.rda"
op3_dga_path2="output/op3_dga_analysis2.rda"
out_path="output/resulting_excels/op3_dga.xlsx"
out_path2="output/resulting_excels/op3_dgaXpat.xlsx"

#Data loading ----
#Load the batch corrected count matrix
load(op1_data_path)
gs2keep=fread(gs2keep_path,data.table = F,header = F)
gs2keep=unique(gs2keep)

#Combine source with isto
comb_cl=paste(infoc$source,infoc$isto,sep="_")
infoc$comb_cl=comb_cl

#Get group indexes ----
classes=unique(infoc$comb_cl)
ref_groups=infoc$comb_cl
gr_l=list();
for(k in 1:length(classes)){
  cl=classes[k]
  indxs=grep(paste("^",cl,"$", sep=""), infoc$comb_cl)
  gr_l[[cl]]=indxs
}

#Define set of tuning parameters ----
probs_tol=cbind(seq(0.15,0.30,0.025),
                seq(0.85,0.70,-0.025))
ths1=seq(60,90,5)
ths2=seq(40,70,5)
combs=expand.grid(1:nrow(probs_tol),1:length(ths1))
n_maxk=11
ncores=16
prev_var2=0
top_genes1_l=top_genes2_l=list();rl1_l=rl2_l=list();

x=unique(infoc$col_source);names(x)=unique(infoc$source)
isto_df=unique(infoc[,c(2,8)])
y=isto_df[,2];names(y)=isto_df[,1]
ann_df=infoc[,c(1,2)];ann_df$source=as.factor(ann_df$source)
ref=as.numeric(as.factor(paste(ann_df$source,ann_df$isto,sep="_")))
col_list=list(source=x,isto=y)

for(k in 1:nrow(combs)){
  cat(k,"on",nrow(combs),"\n")
  comb=unlist(combs[k,])
  var2=comb[2]
  prob_tol=probs_tol[comb[2],]
  th1=ths1[comb[1]]
  th2=ths2[comb[1]]
  
  #Discriminative analysis of genes and patients ----
  if(prev_var2!=var2){
    cat("  doing discriminative analysis \n")
    res=do_DiscAnalysis(lcpm,gr_l,prob_tol=prob_tol)
    #the table of gene x number of contrasts that it separates
    perc_discr_df=res$perc_discr_df
    #the list containing how much the patients of a class separate from the median for each gene
    dist_med_dfs_l=res$dist_med_dfs_l
  }
  prev_var2=var2
  #Get the genes which separate at least two classes by at least th in their contrasts
  top_genes1=unique(c(get_topDiscGenes1(perc_discr_df,th = th1),gs2keep$V1))
  #Get the genes which separate most of classes by at least th in their contrasts
  top_genes2=unique(c(get_topDiscGenes2(perc_discr_df,th = th2),gs2keep$V1))
  
  #Consensus clustering1 ----
  if(length(top_genes1)>=100000 & length(top_genes1)<=2000){
    cat("  doing cc1 \n")
    #Prepare matrix
    mat1=adjust_matrix(lcpm[match(top_genes1,rownames(lcpm)),],sd_quantile = 0.0001)
    #Set number of genes
    n_top_genes1=c(round(nrow(mat1)/4),round(nrow(mat1)/2),nrow(mat1))
    #Clustering
    rl1=run_all_consensus_partition_methods(
      mat1, max_k = n_maxk, top_n=n_top_genes1, verbose=F, mc.cores = ncores,
      anno_col=col_list,anno=ann_df,
      #partition_method=c("kmeans","pam","mclust")
      )
    intr1=get_best_meth(rl1,ref_groups,10)$min_intrs
  }else{
    intr1=1000
  }
  
  #Consensus clustering2 ----
  if(length(top_genes2)>=100 & length(top_genes2)<=2000){
    cat("  doing cc2 \n")
    #Prepare matrix
    mat2=adjust_matrix(lcpm[match(top_genes2,rownames(lcpm)),],sd_quantile = 0.0001)
    #Set number of genes
    n_top_genes2=c(round(nrow(mat2)/4),round(nrow(mat2)/2),nrow(mat2))
    #Clustering
    rl2=run_all_consensus_partition_methods(
      mat2, max_k = n_maxk, top_n=n_top_genes2, verbose=F, mc.cores = ncores,
      anno_col=col_list,anno=ann_df,
      #partition_method=c("kmeans","pam","mclust")
      )
    intr2=get_best_meth(rl2,ref_groups,10)$min_intrs
  }else{
    intr2=1000
  }
  if(k == 1){
    setting_df=data.frame(prob_tol=I(list(prob_tol)),th1=th1,th2=th2,intr1=intr1,intr2=intr2)
  }else{
    setting_df=rbind(setting_df,
                     data.frame(prob_tol=I(list(prob_tol)),th1=th1,th2=th2,intr1=intr1,intr2=intr2))
  }
  print(setting_df)
  
  top_genes1_l[[k]]=top_genes1
  top_genes2_l[[k]]=top_genes2
  rl2_l[[k]]=rl2
}
save(setting_df,top_genes1_l,top_genes2_l,rl1_l,rl2_l,file=op3_dga_path1)
#load(op3_dga_path1)

#Get the best settings and all their discriminative genes
min_intr2=min(setting_df$intr2)
setting_df=setting_df[min_intr2==setting_df$intr2,]
top_genes2=unique(unlist(top_genes2_l[as.numeric(rownames(setting_df))]))
#Get the parameters of the best and most restrictive setting
top_prob_tol=setting_df$prob_tol[[1]]
#Get the outlier and outsider analysis of the best and most restrictive setting
res=do_DiscAnalysis(lcpm,gr_l,prob_tol=top_prob_tol)
dist_med_dfs_l=res$dist_med_dfs_l

#Outliers analysis ----
res2=do_outlAnalysis(dist_med_dfs_l)
#the table of gene x how much separates the members of each class
#the genes are order in a column but not between columns
#the gene at top of one class could be at the bottom of another class
gs_outl_df=res2$gs_outl_df
#the table of gene x how much separates the members of each class
#the genes are ordered between columns
#the gene at the top of the matrix is relative high in all the classes
outlier_genes_rank_df=res2$outlier_genes_rank_df

#Outsiders analysis ----
do_outsAnalysisXpat=function(dist_med_dfs_l,top_genes){
  #Evaulate how members of a class separate from the median with the genes
  #that separate well all the classes
  #Take into account only the members and their behaviour on top genes
  dist_med_dfs_l2=lapply(dist_med_dfs_l,function(x){return(x[top_genes,])})
  dist_med_dfs_xPAT=list()
  #Iterate over each matrix of the list
  for(k in 1:length(dist_med_dfs_l2)){
    #Analysis on the members for the top genes
    dist_med2=dist_med_dfs_l2[[k]]
    
    x=dist_med2
    x=x[,order(colMedians(x),decreasing = T)]
    
    z=apply(x,2,function(x){
      paste(names(x)," ",round(x,2),sep="")
    })
    for(col in 1:ncol(z)){
      ord=order(x[,col],decreasing = T)
      z[,col]=z[ord,col]
    }
    rownames(z)=seq(1,nrow(z))
    dist_med_dfs_xPAT[[k]]=z
  }
  names(dist_med_dfs_xPAT)=names(dist_med_dfs_l2)
  return(dist_med_dfs_xPAT)
}
#Table which rank how much a member separates from its class using the best
#genes which separate the classes in the clustering
pats_outs_df=do_outsAnalysis(dist_med_dfs_l,top_genes2)
pats_outs_df_l=do_outsAnalysisXpat(dist_med_dfs_l,top_genes2)

#Writing 
write.xlsx(pats_outs_df, out_path, sheetName="outsiders",
           append=F, row.names = T, col.names = T)
write.xlsx(gs_outl_df, out_path, sheetName="outliers",
           append=T, row.names = T, col.names = T)
write.xlsx(outlier_genes_rank_df, out_path, sheetName="outlier_rank",
           append=T, row.names = T, col.names = T)

for(k in 1:length(pats_outs_df_l)){
  label_table=names(pats_outs_df_l)[k]
  x=pats_outs_df_l[[k]]
  x=x[1:300,]
  if(k==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(x, out_path2, sheetName=label_table,
             append=app, row.names = T, col.names = T)
}

save(top_genes2,pats_outs_df,gs_outl_df,outlier_genes_rank_df,pats_outs_df_l,file=op3_dga_path2)



