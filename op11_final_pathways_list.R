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
out_op13_path1="output/op9_pathways_overex_msigdb_c6.rda"
out_op13_path2="output/op9_pathways_overex_msigdb_h.rda"
#output
out_op14_path_xlsx="output/resulting_excels/op11_final_pathways_list.xlsx"

#Load data ----
load(out_op13_path1)
gsva_res_sc1=gsva_res_spe
load(out_op13_path2)
gsva_res_sc2=gsva_res_spe

#Compose the final report of pathways
pathway_report=rbind(gsva_res_sc1,gsva_res_sc2)
write.xlsx(pathway_report,file=out_op14_path_xlsx,sheetName = "pathway_report",
           row.names = T, col.names = T)


