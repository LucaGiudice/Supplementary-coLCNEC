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
library('PCAtools')

check_elems = function(x,y){
  #Check if the vectors have the same length
  if(length(x)==length(y)){
    
    #Check if all elements in y are in x
    indx_ord=match(y,x)
    if(sum(!is.na(indx_ord))==0){
      cat("ERROR:Second vector has something that doesn't appear in the first \n")
      return(0)
    }
    
    #Check if all elements in x are in y
    indx_ord=match(x,y)
    if(sum(!is.na(indx_ord))==0){
      cat("ERROR:First vector has something that doesn't appear in the second \n")
      return(0)
    }
    
    #Check if all elements are also in the same position
    corr_ord=seq(1,length(y))
    err=0
    for(k in corr_ord){
      check=indx_ord[k]==k
      if(is.na(check)){
        err=c(err,k)
      }else{
        if(!check){
          err=c(err,k)
        }
      }
    }
    err=err[-1]
    
    if(length(err)==0){cat(">>Patient's genetic profiles and clinical data match \n");return(1)}
    if(length(err)!=0){cat("ERROR:Patient's genetic profiles and clinical data not match \n");return(0)}
  }else{
    cat("ERROR:Patient's genetic profiles and clinical data have different lengths \n");return(0)
  }
}
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
out_op12_path="output/op12/reduction_plot.png"

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

#Prepare the data for plot-----
hlcpm=rbind(must_lcpm,cl_lcpm)
check_elems(colnames(hlcpm),infoc$case)

infoc$isto[infoc$source=="CoLCNEC"]="CoLCNEC"
infoc$col_isto[infoc$source=="CoLCNEC"]="#008080"
infoc$cluster=paste("CL",infoc$cluster,sep="")

map_col_cl=cbind(paste("CL",seq(1,10),sep=""),
                 c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                   "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))
infoc$col_cl=mapvalues(infoc$cluster,from = map_col_cl[,1], to=map_col_cl[,2])
infoc=infoc[,c(3,2,8,11,12)]

#PCA ----
rownames(infoc)=infoc$case
p <- pca(mat=hlcpm, metadata = infoc, center=T, scale=F, rank=20)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

colkey_v=infoc$col_isto;names(colkey_v)=infoc$isto;colkey_v[!duplicated(colkey_v)]
shapekey_v=seq(1,10);names(shapekey_v)=unique(infoc$isto)
png(out_op12_path, units="in", width=22, height=12, res=300)
biplot(p,
       legendPosition = 'top', colLegendTitle = "isto groups",
       colby = "isto", colkey = colkey_v,
       shapeLegendTitle = "clusters", shape = "cluster", shapekey = shapekey_v,
       drawConnectors = FALSE,
       pointSize = 6, labSize = 4, labhjust = 1 
       )
dev.off()






