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
out_op14_path="output/op10_final_genes.rda"
out_op13_path1="output/op9_overexp_excl1_analysis.rda"
out_op13_path2="output/op9_overexp_excl2_analysis.rda"
out_op6_path="output/op5_DE_results.rda"
#output
out_op14_path_xlsx="output/resulting_excels/op11_final_gene_list.xlsx"

#Load data ----
load(out_op14_path)
load(out_op6_path)
load(out_op13_path1)
stats_gs1=stats_gs;colnames(stats_gs1)=c("gene","excl_cl")
stats_gs1=as.data.frame(stats_gs1,stringsAsFactors = F)
load(out_op13_path2)
stats_gs2=stats_gs;colnames(stats_gs2)=c("gene","excl_pair1","excl_pair2")
stats_gs2=as.data.frame(stats_gs2,stringsAsFactors = F)
stats_gs=merge(stats_gs1,stats_gs2,by="gene",all=T)
rownames(stats_gs)=stats_gs[,1];stats_gs=stats_gs[,-1]

#Get final genes of interest
gs=rownames(hlcpm)

for(kc in 1:length(de_contrast_l)){
  cat(kc,"in:",length(de_contrast_l),"\n")
  tT_ov=de_contrast_l[[kc]]$tT_ov[,c(1,2,5)]
  tT_ov$DE=tT_ov$adj.P.Val
  tT_ov$DE[tT_ov$adj.P.Val<=0.05]=1
  tT_ov$DE[tT_ov$adj.P.Val>0.05]=0
  
  missing_gs=setdiff(gs,rownames(tT_ov))
  add_m=matrix(NA,nrow=length(missing_gs),ncol=ncol(tT_ov))
  rownames(add_m)=missing_gs
  tT_ov=as.matrix(tT_ov)
  tT_ov=rbind(tT_ov,add_m)
  tT_ov=as.data.frame(tT_ov,stringsAsFactors = F)
  colnames(tT_ov)=paste(colnames(tT_ov),kc,sep="_n")
  
  if(kc==1){
    tT_ov_gs=tT_ov[gs,]
    df_gs_stats=tT_ov_gs
  }else{
    df_gs_stats=merge(df_gs_stats,tT_ov,by="row.names")
    rownames(df_gs_stats)=df_gs_stats[,1]
    df_gs_stats=df_gs_stats[,-1]
  }
}

m_genes=merge(df_gs_stats,stats_gs,by="row.names",all.x=T)
rownames(m_genes)=m_genes[,1];m_genes=m_genes[,-1]
m_genes_sub=m_genes[,c(181,182,183)]
m_genes=m_genes[,-c(181,182,183)];m_genes=cbind(m_genes_sub,m_genes)


header_df=c("exclusive_genes","exclusive_genes","exclusive_genes",
            rep(names(de_contrast_l),each=ncol(tT_ov)))
header_df=rbind(header_df,colnames(m_genes))
df_m_genes=as.matrix(m_genes)
full_m <- data.frame()
full_m <- rbind(header_df, df_m_genes)
rownames(full_m)[1:2]=c("header1","header2")

write.xlsx(full_m,file = out_op14_path_xlsx,row.names = T,
           col.names = F)



