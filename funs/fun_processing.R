Quality_check=function(df){
  m=as.matrix(df)
  m[m == 0] <- NA
  m <- log2(m)
  boxplot(m,
          ylab = "Intensity Distribution",
          xlab = "Samples")
}

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

make_venn=function(set_A,set_B,name_A,name_B,title_name,cex_int=2){
  require(gridExtra)
  require(VennDiagram)
  common_set=intersect(set_A,set_B)
  n_common_set=length(common_set)
  venn.plot <- draw.pairwise.venn(length(set_A), length(set_B), 
                                  n_common_set, c(name_A, name_B), scaled = TRUE,
                                  fill = c("blue", "green"), euler.d=TRUE, 
                                  ind = FALSE, cex=cex_int, cat.cex=cex_int)
  grid.arrange(gTree(children=venn.plot), top=textGrob(title_name, gp=gpar(fontsize=15)))
}

rfilt_fun = function(m){
  df=m
  m=as.matrix(m[,c(-1,-2)])
  rsd=rowSds(m)>1
  rme=rowMedians(m)>1
  rfilt=rme & rsd
  df=df[rfilt,]
  return(df)
}

rem_dups = function(m){
  m=m[,-1]
  gs=m$gene_name
  dups=unique(gs[duplicated(gs)])
  for(k in 1:length(dups)){
    gs=m$gene_name
    g=dups[k]
    indxs=grep(g,gs,fixed = T)
    m1=as.matrix(m[indxs,-1])
    m=m[-indxs,]
    ex1=rowMedians(m1)>=1
    if(sum(ex1)!=0){
      if(sum(ex1)==1){
        prof=m1[ex1,]
      }else{
        prof=colMedians(m1[ex1,])
      }
    }else{
      prof=colMedians(m1)
    }
    
    prof=c(g,prof)
    if(k==1){
      prof_df=prof
    }else{
      prof_df=rbind(prof_df,prof)
    }
  }
  colnames(prof_df)=colnames(m)
  m_new=rbind(m,prof_df)
  return(m_new)
}

check_elems = function(x,y){
  indx_ord=match(x,y)
  if(sum(!is.na(indx_ord))==0){cat("NO MATCHES");return(0)}
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
  
  if(length(err)==0){cat("Matching \n");return(0);}
  if(length(err)!=0){cat("NO Matching, returning indexes \n");return(err)}
}

check_dist = function(RNA){
  m=as.matrix(RNA);m=apply(m,2,as.numeric)
  print(colQuantiles(m[,1:20]))
}