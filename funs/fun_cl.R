get_best_meth = function(rl,ref,topn=5){
  top_methods=rownames(suggest_best_k(rl))[1:topn]
  top_dists=sub("\\:.*", "", top_methods)
  top_cls=sub('.*:', '', top_methods)
  
  for(i in 1:topn){
    meth=top_cls[i]
    dist=top_dists[i]
    obj1 = rl[dist, meth]
    df_cl_classes=get_classes(obj1)
    for(k in 1:ncol(df_cl_classes)){
      cl_v1=df_cl_classes[,k]
      mono_intrusi=as.data.frame(table(cl_v1))$Freq
      mono_intrusi=sum(mono_intrusi[mono_intrusi<=2])
      cls=unique(cl_v1)
      for(cl1 in cls){
        indxs=which(cl_v1==cl1)
        groups=ref[indxs]
        groups_freq=as.data.frame(table(groups))
        intrusi=sum(groups_freq$Freq)-groups_freq$Freq[which.max(groups_freq$Freq)]
        if(cl1==cls[1]){
          intrusi_tot=intrusi+mono_intrusi
        }else{
          intrusi_tot=intrusi_tot+intrusi
        }
      }
      if(k==1){
        adj_ri_v=intrusi_tot
      }else{
        adj_ri_v=c(adj_ri_v,intrusi_tot)
      }
    }
    if(i==1){
      df_adj_ri_v=adj_ri_v
    }else{
      df_adj_ri_v=cbind(df_adj_ri_v,adj_ri_v)
    }
  }
  rownames(df_adj_ri_v)=seq(1:ncol(df_cl_classes))
  colnames(df_adj_ri_v)=top_methods
  top_method=colnames(df_adj_ri_v)[colSums(df_adj_ri_v==min(df_adj_ri_v))!=0][1]
  top_dist=sub("\\:.*", "", top_method)
  top_cl=sub('.*:', '', top_method)
  top_k=rownames(df_adj_ri_v)[rowSums(df_adj_ri_v==min(df_adj_ri_v))!=0][1]
  top_k=as.numeric(top_k)+1
  min_intrs=min(df_adj_ri_v)
  cat("dist:",top_dist," top_cl:",top_cl," top_k:",top_k," min_intr:",min_intrs,"\n",sep="")
  res=list(top_dist=top_dist,top_cl=top_cl,top_k=top_k,min_intrs=min_intrs)
  return(res)
}
do_DiscAnalysis=function(m,gr_l,prob_tol=c(0.2,0.8)){
  cont_tol=0
  cont_tol_v=0
  combs_l=dist_med_dfs_l=list()
  #Create contrasts for each class 
  # e.g 1vs2, 1vs3, etc
  # e.g 2vs1, 2vs3, etc
  for(k in 1:length(gr_l)){
    combs_m=cbind(
      rep(k,length(gr_l)-1),
      setdiff(seq(1,length(gr_l)),k)
    )
    combs_l[[k]]=combs_m
  }
  
  #For every set of contrasts: 
  # get how much each gene separates a class from the others 
  # get how much a patient goes far away from the group for each gene
  for(comb in 1:length(combs_l)){
    #cat(comb,"on",7,"\n")
    #Get the contrasts of a specific class
    combs_m=combs_l[[comb]]
    
    for(k in 1:nrow(m)){
      #Get the expression of a gene
      g_v=m[k,]
      
      for(i in 1:length(gr_l)){
        #Get the distribution of the gene in a specific class
        g1_v=g_v[gr_l[[i]]]
        q1_v=quantile(g1_v,prob=prob_tol)
        #Create the table of distributions of the gene
        #for all the classes
        if(i==1){
          q1_df=q1_v
        }else{
          q1_df=rbind(q1_df,q1_v)
        }
        #Get the median of the gene inside a class
        med_g1_v=median(g1_v)
        #Compute how much the patients tend to separate from the med.
        dist_med_g1_v=abs(g1_v-med_g1_v)
        #Create the df of how the patient of each class separate with
        #respect the same gene
        if(k==1 && comb==1){
          dist_med_dfs_l[[i]]=dist_med_g1_v
        }else{
          if(comb==1){
            dist_med_dfs_l[[i]]=rbind(dist_med_dfs_l[[i]],dist_med_g1_v)
          }
        }
      }
      rownames(q1_df)=names(gr_l)
      
      #For each gene df of distributions x class
      #Compare due to the contrast and detect how many times
      #the gene doesn't overlap
      for(j in 1:nrow(combs_m)){
        gr1=combs_m[j,1];gr2=combs_m[j,2]
        eval=(q1_df[gr1,1]>q1_df[gr2,2]) || (q1_df[gr2,1]>q1_df[gr1,2])
        if(eval==T){
          cont_tol=cont_tol+1
        }
      }
      cont_tol_v=c(cont_tol_v,cont_tol)
      cont_tol=0
    }
    cont_tol_v=cont_tol_v[-1]
    if(comb==1){
      cont_tol_df=cont_tol_v
    }else{
      cont_tol_df=cbind(cont_tol_df,cont_tol_v)
    }
    cont_tol_v=0
  }
  
  #Format and finish the table of gene x number of contrasts that
  #it separates
  colnames(cont_tol_df)=names(gr_l)
  rownames(cont_tol_df)=rownames(m)
  n_comps=(length(gr_l)-1)
  perc_tol_df=t(apply(cont_tol_df,1,function(x){
    y=(100*x)/n_comps
    return(y)
  }))
  
  #Format and finish the list containing how much the patients of a
  #class separate from the median for each gene
  names(dist_med_dfs_l)=names(gr_l)
  dist_med_dfs_l=lapply(dist_med_dfs_l,function(x){
    rownames(x)=rownames(m)
    return(x)
  })
  
  res=list(perc_discr_df=perc_tol_df,dist_med_dfs_l=dist_med_dfs_l)
  return(res)
}
get_topDiscGenes1=function(perc_discr_df,th=80){
  bol_df=t(apply(perc_discr_df,1,function(x){
    y=x>th
    return(y)
  }))
  rowS=rowSums2(bol_df)
  bol_top_genes=rowS>1
  top_genes=rownames(bol_df)[bol_top_genes]
  return(top_genes)
}
get_topDiscGenes2=function(perc_discr_df,th=60){
  rowM=rowMeans(perc_discr_df)
  bol_top_genes=rowM>th
  top_genes=rownames(perc_discr_df)[bol_top_genes]
  return(top_genes)
}
do_outlAnalysis=function(dist_med_dfs_l){
  nrow_df=nrow(dist_med_dfs_l[[1]])
  for(k in 1:length(dist_med_dfs_l)){
    gs_v=rep("-",nrow_df)
    #Now process how all genes behave in seperating the memembers of each class
    dist_med1=dist_med_dfs_l[[k]]
    #Compute how much the members separate for each gene
    outlier_genes=rowMedians(dist_med1)
    #Compute the rank of the genes in the original matrix order s.t. I can compare
    #their rank between their behaviour in the classes
    outlier_genes_rank=rank(outlier_genes,ties.method = "min")/length(outlier_genes)
    #Process the vector with the median value of gene behaviour
    names(outlier_genes)=rownames(dist_med1)
    #Order
    outlier_genes=outlier_genes[order(outlier_genes,decreasing = T)]
    #Format
    outlier_genes_v=paste(names(outlier_genes),round(outlier_genes,2))
    gs_v[1:length(rownames(dist_med1))]=outlier_genes_v
    
    if(k==1){
      gs_outl_df=gs_v
      outlier_genes_rank_df=outlier_genes_rank
    }else{
      gs_outl_df=cbind(gs_outl_df,gs_v)
      outlier_genes_rank_df=cbind(outlier_genes_rank_df,outlier_genes_rank)
    }
  }
  colnames(gs_outl_df)=colnames(outlier_genes_rank_df)=names(dist_med_dfs_l)
  rownames(outlier_genes_rank_df)=rownames(dist_med_dfs_l[[1]])
  rownames(gs_outl_df)=1:nrow(gs_outl_df)
  gs_outl_df=as.data.frame(gs_outl_df)
  rowMR=rowMeans2(outlier_genes_rank_df)
  outlier_genes_rank_df=outlier_genes_rank_df[order(rowMR,decreasing = T),]
  res=list(gs_outl_df=gs_outl_df,outlier_genes_rank_df=outlier_genes_rank_df)
  return(res)
}


do_outsAnalysis=function(dist_med_dfs_l,top_genes){
  #Evaulate how members of a class separate from the median with the genes
  #that separate well all the classes
  #Take into account only the members and their behaviour on top genes
  dist_med_dfs_l2=lapply(dist_med_dfs_l,function(x){return(x[top_genes,])})
  #Set variables
  ncol_df=max(sapply(dist_med_dfs_l2,function(x){ncol(x)}))
  #Iterate over each matrix of the list
  for(k in 1:length(dist_med_dfs_l2)){
    pats_v=rep("-",ncol_df)
    #Analysis on the members for the top genes
    dist_med2=dist_med_dfs_l2[[k]]
    #Compute how much each member separates from the medians of the top genes
    outsider_patients=colMedians(dist_med2)
    names(outsider_patients)=colnames(dist_med2)
    #Order from the worst memeber to the most representative one
    outsider_patients=outsider_patients[order(outsider_patients,decreasing = T)]
    #Format the vector in order to keep the order, name and how in average he seps
    outsider_patients_v=paste(names(outsider_patients),round(outsider_patients,2))
    pats_v[1:length(colnames(dist_med2))]=outsider_patients_v
    
    if(k==1){
      pats_outs_df=pats_v
    }else{
      pats_outs_df=cbind(pats_outs_df,pats_v)
    }
  }
  colnames(pats_outs_df)=names(dist_med_dfs_l2)
  rownames(pats_outs_df)=1:nrow(pats_outs_df)
  pats_outs_df=as.data.frame(pats_outs_df)
  return(pats_outs_df)
}