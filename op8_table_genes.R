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
get_test_diff = function(lcpm,clusters,g_name="OTP",cl_int="CL2"){
  n_comp=0
  indxs=grep(paste("^",g_name,"$", sep=""), rownames(lcpm))
  if(length(indxs)!=0){
    g_row=lcpm[indxs,]
    g_df=data.frame(group=clusters,value=g_row,stringsAsFactors = T)
  
    PT = FSA::dunnTest(value ~ group, data=g_df, method="bh")
    cls_tests=PT[["res"]]
    cls_tests=cls_tests[cls_tests$P.adj<=0.05,]
  
    if(nrow(cls_tests)==0){
      n_comp=0
    }else{
      comp_names=strsplit(as.character(cls_tests$Comparison),split = " ",fixed = T)
      n_comp=sum(sapply(comp_names,function(x){
        pos_ch=grep(paste("^",cl_int,"$", sep=""), x)
        check=F
        if(length(pos_ch)!=0){
          check=T
        }else{
          check=F
        }
        return(check)
      }))
    }
  }
  return(n_comp)
}
get_r_exp = function(lcpm,clusters,g_name,cl_int){
  r_exp=0
  indxs=grep(paste("^",g_name,"$", sep=""), rownames(lcpm))
  if(length(indxs)!=0){
    g_row=lcpm[indxs,]
    g_df=data.frame(group=clusters,value=g_row,stringsAsFactors = T)
    
    exp_v=tapply(g_df$value, g_df$group, median)
    r_exp_v=rank(-exp_v)
    r_exp=r_exp_v[grep(paste("^",cl_int,"$", sep=""), names(r_exp_v))]
  }
  return(r_exp)
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
out_op4_path="output/op4_data_prep.rda"
out_op5_path="output/op5_DE_results.rda"
out_op6_path="output/op6_GS_analysis.rda"
#Output
out_path5="output/resulting_excels/AAA_Contrasts.xlsx"
out_path6a="output/resulting_excels/op8_DE_int_genes_results.xlsx"
out_path6aa="output/resulting_excels/op8_DE_int_top_genes_results.xlsx"
out_path6b="output/resulting_excels/op8_DE_int_gene_enrichment_results.xlsx"
out_path7a="output/resulting_excels/op8_GS_int_SPIA_results.xlsx"
out_path7b="output/resulting_excels/op8_GS_int_GSEA_results.xlsx"
out_path7c="output/resulting_excels/op8_GS_int_camera_results.xlsx"

#Load data ----
load(out_op4_path)
load(out_op5_path)

#Writing report ----
labels=names(de_contrast_l)
groups=paste("CL",unique(infoc$cluster),sep="")
check_elems(infoc$isi_names,colnames(lcpm))

cls_de_genes_stats=list()
for(k_cl in 1:length(groups)){
  cat(k_cl,"/",length(groups),"\n")
  group=groups[k_cl]
  if(k_cl==1){
    indxs_in=unique(c(1,setdiff(grep(group,labels,value=F),grep("CL10",labels,value=F))))
  }else{
    indxs_in=grep(group,labels,value=F)
  }
  indxs_out=setdiff(seq(1,length(labels)),indxs_in)
  de_in_l=de_contrast_l[indxs_in]
  de_out_l=de_contrast_l[indxs_out]
  
  de_gs_in=sapply(de_in_l,function(x){
    gs=rownames(x$tT)
    return(gs)
  })
  de_gs_in=unlist(de_gs_in)
  de_gs_in=as.data.frame(table(de_gs_in))
  de_gs_in$de_gs_in=as.character(de_gs_in$de_gs_in)
  de_gs_in=de_gs_in[order(de_gs_in$Freq,decreasing = T),]
  de_gs_in$Freq=(de_gs_in$Freq*100)/length(de_in_l)
  colnames(de_gs_in)[1]="gs"
  
  de_gs_out=sapply(de_out_l,function(x){
    gs=rownames(x$tT)
    return(gs)
  })
  de_gs_out=unlist(de_gs_out)
  de_gs_out=as.data.frame(table(de_gs_out))
  de_gs_out=de_gs_out[order(de_gs_out$Freq,decreasing = T),]
  de_gs_out$de_gs_out=as.character(de_gs_out$de_gs_out)
  de_gs_out$Freq=(de_gs_out$Freq*100)/length(de_out_l)
  colnames(de_gs_out)[1]="gs"
  gs_df=merge(x=de_gs_in,y=de_gs_out,by="gs",all.x = T, all.y = F)
  gs_df$Freq.y[is.na(gs_df$Freq.y)]=0
  colnames(gs_df)=c("genes","score_in","score_out")
  
  cl_int=group
  clusters=as.factor(paste("CL",infoc$cluster,sep=""))
  for(k_g in 1:nrow(gs_df)){
    g_name1=gs_df[k_g,1]
    test_score=(get_test_diff(lcpm,clusters=clusters,g_name=g_name1,cl_int=cl_int)*100)/(length(groups)-1)
    r_exp=get_r_exp(lcpm,clusters,g_name1,cl_int)
    if(k_g==1){
      test_score_v=test_score
      r_exp_v=r_exp
    }else{
      test_score_v=c(test_score_v,test_score)
      r_exp_v=c(r_exp_v,r_exp)
    }
  }
  gs_df$test_score=test_score_v
  gs_df$r_exp=r_exp_v
  cls_de_genes_stats[[cl_int]]=gs_df
}
save(cls_de_genes_stats,file="output/op8_de_gs_analysis1.rda")

for(k_el in 1:length(cls_de_genes_stats)){
  cat(k_el,"on",length(cls_de_genes_stats),"\n")
  group=names(cls_de_genes_stats)[k_el]
  gs_df1=cls_de_genes_stats[[k_el]]
  out_gs_dfs=setdiff(seq(1,length(cls_de_genes_stats)),k_el)
  
  for(k_g in 1:nrow(gs_df1)){
    g_int=gs_df1$genes[k_g]
    max_score=0
    max_name_cl=name_cl="NA"
    
    for(k_el_out in out_gs_dfs){
      name_cl=names(cls_de_genes_stats)[k_el_out]
      gs_df2=cls_de_genes_stats[[k_el_out]][,c(1,2)]
      indx=match(g_int,gs_df2$genes)
      if(!is.na(indx)){
        score=gs_df2$score_in[indx]
      }else{
        score=0
      }
      
      if(score>max_score){
        max_score=score
        max_name_cl=name_cl
      }
    }
  
    if(k_g==1){
      max_score_v=max_score
      max_name_cl_v=max_name_cl
    }else{
      max_score_v=c(max_score_v,max_score)
      max_name_cl_v=c(max_name_cl_v,max_name_cl)
    }
  }
  
  gs_df1$max_in_score_out=max_score_v
  gs_df1$max_in_cl_out=max_name_cl_v
  cls_de_genes_stats[[k_el]]=gs_df1
  
  if(k_el==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(gs_df1, out_path6a, sheetName=group, append=app, row.names = F, col.names = T)
  
  x=gs_df1
  x=x[x$score_in>=quantile(x$score_in,probs = 0.80) &
        x$score_out<=quantile(x$score_out,probs = 0.20) & 
        x$test_score>=quantile(x$test_score,probs = 0.80) &
        x$r_exp<=quantile(x$r_exp,probs = 0.20) &
        x$max_in_score_out<=quantile(x$max_in_score_out,probs = 0.20),]
  x=x[order(x$score_out,decreasing = F),]
  
  write.xlsx(x, out_path6aa, sheetName=group, append=app, row.names = F, col.names = T)
  
}
save(cls_de_genes_stats,file="output/op8_de_gs_analysis2.rda")


for(k_cl in 1:length(groups)){
  cat(k_cl,"/",length(groups),"\n")
  group=groups[k_cl]
  if(k_cl==1){
    indxs_in=unique(c(1,setdiff(grep(group,labels,value=F),grep("CL10",labels,value=F))))
  }else{
    indxs_in=grep(group,labels,value=F)
  }
  indxs_out=setdiff(seq(1,length(labels)),indxs_in)
  
  de_in_l=de_contrast_l[indxs_in]
  de_out_l=de_contrast_l[indxs_out]

  de_gs_in=sapply(de_in_l,function(x){
    gs=rownames(x$camera_res)
    return(gs)
  })
  de_gs_in=unlist(de_gs_in)
  de_gs_in=as.data.frame(table(de_gs_in))
  de_gs_in$de_gs_in=as.character(de_gs_in$de_gs_in)
  de_gs_in=de_gs_in[order(de_gs_in$Freq,decreasing = T),]
  de_gs_in$Freq=(de_gs_in$Freq*100)/length(de_in_l)
  colnames(de_gs_in)[1]="gs_camera"
  
  de_gs_out=sapply(de_out_l,function(x){
    gs=rownames(x$camera_res)
    return(gs)
  })
  de_gs_out=unlist(de_gs_out)
  de_gs_out=as.data.frame(table(de_gs_out))
  de_gs_out=de_gs_out[order(de_gs_out$Freq,decreasing = T),]
  de_gs_out$de_gs_out=as.character(de_gs_out$de_gs_out)
  de_gs_out$Freq=(de_gs_out$Freq*100)/length(de_out_l)
  colnames(de_gs_out)[1]="gs_camera"
  gs_df=merge(x=de_gs_in,y=de_gs_out,by="gs_camera",all.x = T, all.y = F)
  gs_df$Freq.y[is.na(gs_df$Freq.y)]=0
  
  scores_m=gs_df[,c(2,3)]
  scores_m[,1]=scales::rescale(scores_m[,1],to=c(0,1))
  scores_m[,2]=scales::rescale(scores_m[,2],to=c(1,0))
  scores_m=as.matrix(scores_m)
  score_fin=apply(scores_m,1,autri)
  gs_df$score_sum=score_fin
  gs_df=gs_df[order(gs_df$score_sum,decreasing = T),]
  colnames(gs_df)=c("gs_camera","score_in","score_out","score_sum")
  if(k_cl==1){
    app=F
  }else{
    app=T
  }
  write.xlsx(gs_df, out_path6b, sheetName=group, append=app, row.names = F, col.names = T)
  gc()
}