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
library("readxl")
library("cola")
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
adj.rand.index = function(cl1, cl2){
  tab <- table(cl1, cl2)
  f2 <- function(n) n * (n - 1)/2
  sum.f2 <- function(v) sum(f2(v))
  marg.1 <- apply(tab, 1, sum)
  marg.2 <- apply(tab, 2, sum)
  n <- sum(tab)
  prod <- sum.f2(marg.1) * sum.f2(marg.2)/f2(n)
  num <- (sum.f2(as.vector(tab)) - prod)
  den <- 0.5 * (sum.f2(marg.1) + sum.f2(marg.2)) - prod
  return(num/den)
}
#Set variables ----
#Input variables
out_op1_path="data/op1_data_rdy.rda"
op2_data_path1="./output/op4_CL1.rda"
op2_data_path2="./output/op4x_artef_analysis.rda"
#Output variables
out_op4_path="./output/op4_data_prep.rda"

#Load data ----
load(out_op1_path)
load(op2_data_path1)
load(op2_data_path2)
#Extract the cluster membership labels
cls_lab=paste("CL",get_classes(res, k = top_k)$class,sep="")
infoc$cls=cls_lab
infoc$comb_names=paste(substr(infoc$source,1,2),infoc$isto,sep="_")

#Order the info matrix and as consequence the counts matrix ----
infoc=infoc[order(infoc$isto,decreasing = F),]
lcpm=lcpm[,infoc$case]
counts=counts[,infoc$case]
check_elems(infoc$case,colnames(lcpm))
check_elems(infoc$case,colnames(counts))
lcpm=lcpm[-match(artefatti_gs,rownames(lcpm)),]
counts=counts[-match(artefatti_gs,rownames(counts)),]

#Simplify the name of the samples
isi_names=make.unique(infoc$comb_names,sep="")
infoc$isi_names=isi_names
colnames(lcpm)=isi_names
colnames(counts)=isi_names
check_elems(infoc$isi_names,colnames(lcpm))
check_elems(infoc$isi_names,colnames(counts))

#Get how the classes are organized by isto ----
indxs_byCL_l=list()
#For each condition
for(cl in unique(infoc$cls)){
  indexs=grep(paste("^",cl,"$", sep=""), infoc$cls)
  indxs_byCL_l[[cl]]=indexs
}

#Create contrasts -----
contrasts=data.frame(t(combn(unique(infoc$cls), 2)), stringsAsFactors = F)
colnames(contrasts)=c("A","B")

#Create the data for each contrast ----
dataXcontrast_l=list()
for(k in 1:nrow(contrasts)){
  cat(k,"on",nrow(contrasts),"\n")
  #Select the class to create the data for comparison
  Acl=contrasts$A[k]
  Bcl=contrasts$B[k]
  #Create bin contrast and submatrix of info and lcpm
  indxs1=indxs_byCL_l[[Acl]]
  indxs2=indxs_byCL_l[[Bcl]]
  indxs=c(indxs1,indxs2)
  sub_lcpm=lcpm[,indxs]
  sub_counts=counts[,indxs]
  sub_infoc=infoc[indxs,]
  bin_contrast=ifelse(sub_infoc$cls == Acl, 1, 0)
  #Create design matrix
  batch=droplevels(sub_infoc$batch)
  Condition=sub_infoc$cls
  if(length(unique(batch))!=1 & adj.rand.index(as.character(batch),Condition)<=0.7){
    design=model.matrix(~0+batch+Condition)
    design=1-design
  }else{
    cat(k,"here \n")
    design=model.matrix(~Condition)
  }
  #Contrast
  contrast_name=paste(Acl,"vs",Bcl,sep="")
  #Pack everything
  sub_infoc$group=bin_contrast
  dataXcontrast_l[[contrast_name]]=list(infoc=sub_infoc,lcpm=sub_lcpm,counts=sub_counts,
                                        design=design,
                                        contrast_name=contrast_name)
}

save(dataXcontrast_l,contrasts,infoc,lcpm,counts,file=out_op4_path)
