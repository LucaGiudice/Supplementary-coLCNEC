see = function(m,nr=5,nc=5,view=FALSE){
  if(nrow(m)<5){
    nr=nrow(m)
  }
  if(ncol(m)<5){
    nc=ncol(m)
  }
  if(view==TRUE){
    View(m[1:nr,1:nc])
  }else{
    m[1:nr,1:nc]
  }
}

#' Function to compute the adj.rand.index between two partitions|clusters
#' Roxygen Documentation:
#' @param cl1 vector of integers, first partition
#' @param cl2 vector of integers, second partition
#' @return measure of similarity between the two partitions
#' @export
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

get_lab_clusters=function(res_list,k_cl){
  m_stats=get_stats(res_list, k = k_cl)
  best_method=names(which.max(m_stats[,4]))
  best_method=substr(best_method,5,20)
  res = res_list["ATC", best_method]
  clusters=paste0(best_method,get_classes(res, k = 2)$class)
  return(clusters)
}

vec2col=function(ov_genes,excl_genes,lab1="OV",lab2="EXCL"){
  gene_lab=rep(lab1,length(ov_genes))
  indxs=match(excl_genes,ov_genes)
  ov_genes=c(ov_genes[-indxs],excl_genes)
  gene_lab=c(gene_lab[-indxs],rep(lab2,length(excl_genes)))
  res=list(genes_ord=ov_genes,gene_lab=gene_lab)
  return(res)
}

vec2col_up=function(ov_genes,excl_genes,lab2="EXCL",gene_lab){
  indxs=match(excl_genes,ov_genes)
  ov_genes=c(ov_genes[-indxs],excl_genes)
  gene_lab=c(gene_lab[-indxs],rep(lab2,length(excl_genes)))
  res=list(ov_genes=ov_genes,gene_lab=gene_lab)
  return(res)
}


color_list = function(ann_df){
  
  ch2col=function(var,var_name,colors,l){
    map2color<-function(x,pal,limits=NULL){
      if(is.null(limits)) limits=range(x)
      pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
    }
    
    na_indxs=var!="NA"
    v=as.numeric(var[na_indxs])
    mypal <- colorRampPalette( colors )(500)
    col_v=map2color(v,mypal)
    col_v=c(rep("#999999",sum(!na_indxs)),col_v)
    names(col_v)=var
    col_v=col_v[unique(names(col_v))]
    l[[var_name]]=col_v
    return(l)
  }
  
  col_vars_l=list()
  #col_vars_l=ch2col(ann_df[,1],colnames(ann_df)[1],c("blue","orange","red"),col_vars_l)
  #col_vars_l=ch2col(ann_df[,3],colnames(ann_df)[3],c("purple","yellow","green"),col_vars_l)
  
  bin2col=function(var,var_name,colors,l){
    names(colors)=unique(var)
    l[[var_name]]=colors
    return(l)
  }
  
  colors=c("#0099FF","#FF6666")
  col_vars_l=bin2col(ann_df[,2],colnames(ann_df)[2],colors,col_vars_l)
  # colors=c("#999999","#993333","#FF33FF")
  # col_vars_l=bin2col(ann_df[,3],colnames(ann_df)[3],colors,col_vars_l)
  # colors=c("#999999","#CCFF33","#00CCCC")
  # col_vars_l=bin2col(ann_df[,6],colnames(ann_df)[6],colors,col_vars_l)
  colors=c("#00CC99","#993333","#FF0033","#FFF000","#0033FF","#9900FF","#FF00CC")
  col_vars_l=bin2col(ann_df[,1],colnames(ann_df)[1],colors,col_vars_l)
  
  return(col_vars_l)
  
}