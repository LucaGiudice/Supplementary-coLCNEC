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
artef_path="data/artefatti_analysis/artefatti.txt"
fake_artef_path="data/artefatti_analysis/fake_artefatti.txt"
#output
out_op4x_path_png="output/op4x/op4x_artef_lcpm.png"
out_op4x_path_rda="output/op4x_artef_analysis.rda"
out_op4x_path_gs="output/op4x/artefatti_gs.txt"

#Load data ----
#Get norm and corr count matrix and info
load(out_op4_path)
rm(mat,res,rl)
#Get manually detected true and false outlier genes
artef_gs=fread(artef_path,data.table = F,header = F)
artef_gs=artef_gs$V1
fake_arte_gs=fread(fake_artef_path,data.table = F,header = F)
fake_arte_gs=fake_arte_gs$V1

#Compute stats about gene distributions in the clusters ----
groups=unique(infoc$cluster)
groups=setdiff(groups,c(5,6,10))
check_elems(infoc$case,colnames(lcpm))

genes=rownames(lcpm)
sel_genes="start"

#Local
for(k_cl in 1:length(groups)){
  cat(k_cl,"/",length(groups),"\n")
  group=groups[k_cl]
  indxs_in=which(infoc$cluster==group)
  lcpm_in=lcpm[,indxs_in]

  
  gr_pos_max_diff=apply(lcpm_in,1,function(x){
    x=x[order(x,decreasing = F)]
    y=diff(x)
    z=(which.max(y)*100)/length(x)
    return(z)
  })

  gr_max_diff=apply(lcpm_in,1,function(x){
    x=x[order(x,decreasing = F)]
    y=diff(x)
    z=max(y)
    return(z)
  })
  
  gr_top_mean=apply(lcpm_in,1,function(x){
    x=x[order(x,decreasing = T)]
    n_top=round((length(x)*3)/10)
    if(n_top<2){n_top=2}
    y=x[seq(1,n_top)]
    z=mean(y)-median(x)
    return(z)
  })
  
  
  gr_skew=apply(lcpm_in,1,function(x){
    z=skew(x)
    return(z)
  })
  
  gr_kurt=apply(lcpm_in,1,function(x){
    z=kurt(x)
    return(z)
  })
  
  gr_lin=apply(lcpm_in,1,function(x){
    x=x[order(x,decreasing = F)]
    y=seq(1,length(x))
    res_lm=lm(x ~ y)
    z=abs(res_lm$coefficients[1])
    return(z)
  })
  
  
  if(k_cl==1){
    pos_jump_df=gr_pos_max_diff
    val_jump_df=gr_max_diff
    diff_med_df=gr_top_mean
    skew_df=gr_skew
    lin_df=gr_lin
    kurt_df=gr_kurt
  }else{
    pos_jump_df=cbind(pos_jump_df,gr_pos_max_diff)
    val_jump_df=cbind(val_jump_df,gr_max_diff)
    diff_med_df=cbind(diff_med_df,gr_top_mean)
    skew_df=cbind(skew_df,gr_skew)
    lin_df=cbind(lin_df,gr_lin)
    kurt_df=cbind(kurt_df,gr_kurt)
  }
}
colnames(pos_jump_df)=colnames(val_jump_df)=colnames(diff_med_df)=
  colnames(skew_df)=colnames(lin_df)=colnames(kurt_df)=paste("CL",groups,sep="")

pos_jump_df=as.data.frame(pos_jump_df)
val_jump_df=as.data.frame(val_jump_df)
diff_med_df=as.data.frame(diff_med_df)
skew_df=as.data.frame(skew_df)
lin_df=as.data.frame(lin_df)
kurt_df=as.data.frame(kurt_df)

#Tune the parameters to maximize the reduction of true outliers and
#minimize the false outliers
bol_pos_jump_df=rowSums(pos_jump_df>=95)>=2
bol_skew_df=rowSums(skew_df>=2.9)<3
bol_val_jump_df=rowSums(val_jump_df>=2)>=3
bol_diff_med_df=rowSums(diff_med_df>=2)>=2
bol_lin_df=rowSums(lin_df>=1.9)<1
bol_keep=bol_lin_df & bol_pos_jump_df & bol_skew_df & bol_val_jump_df & bol_diff_med_df

#Filter
lcpm1=lcpm[bol_keep,]

#load the heatmap supporting data
load(data4heat_path)

#Scale the matrix
mat3=t(scale(t(lcpm1)))

png(out_op4x_path_png, units="in", width=22, height=12, res=600)
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

#Save everything
artefatti_gs=rownames(lcpm1)
write.table(artefatti_gs,file = out_op4x_path_gs, quote = F, row.names = F, col.names = F)
save(artefatti_gs,file=out_op4x_path_rda)

