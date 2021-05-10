#This function:
#Find all the interaction involving genes_interest
#It returns their neighbours
#It returns their interactions with only the genes_2select

find_inter = function(genes_interest,string_db,genes_2select,no_cores){
  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  
  res.l=list()
  res.l=foreach(k = 1:length(genes_interest)) %dopar% {
    #Find the interactions of one selected gene
    g=genes_interest[k]
    #In the first column
    indxs1=grep(g, string_db[,1]);
    #In the second column
    indxs2=grep(g, string_db[,2]);
    
    #Extract the subnetworks
    net1=string_db[indxs1,c(1,2)];net2=string_db[indxs2,c(2,1)];
    colnames(net2)=c("main_protein","found_protein")
    colnames(net1)=c("main_protein","found_protein")
    #Assemble them
    net=rbind(net1,net2);rm(net1,net2);
    
    #Find the neighbours
    neigh=net[,2];
    #Find if some mirna targets are neighbours
    inters=match(genes_2select,neigh)
    inters=unique(inters[!is.na(inters)])
    #Extract the subnetwork
    inters_net=net[inters,]
    if(nrow(inters_net)==0){
      return(NA)
    }else{
      return(list(inters_net=inters_net,neigh=neigh))
    }
  }
  
  #Close cores ----
  stopCluster(cl)
  
  return(res.l)
}

#This is function is combined with find_inter()
#It returns the list without empty elements
keep_exist=function(l.res.1D){
  exists_net=sapply(l.res.1D,function(x){
    class(x[[1]])=="data.frame"
  }
  )
  l.res.1D=l.res.1D[exists_net]
  return(list(l=l.res.1D,bol_exists=exists_net))
}

#This function remove:
#Remove rows with node1 - node1
#Remove duplicated rows
remove_inter = function(edge_list,no_cores){
  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  
  #Remove rows with node1 - node1
  bol_keep=parallel::parApply(cl,edge_list,1,function(x){
    x[1]!=x[2]
  })
  edge_list=edge_list[bol_keep,]
  
  #Remove duplicated rows
  dat.sort = t(parallel::parApply(cl,edge_list, 1, sort))
  edge_list=edge_list[!duplicated(dat.sort),]
  
  stopCluster(cl)
  return(edge_list)
}


#This function:
#Find all the interaction involving genes_interest
find_neigh = function(genes_interest,string_db,no_cores){
  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  
  res.l=list()
  res.l=foreach(k = 1:length(genes_interest)) %dopar% {
    #Find the interactions of one selected gene
    g=genes_interest[k]
    #In the first column
    indxs1=grep(g, string_db[,1]);
    #In the second column
    indxs2=grep(g, string_db[,2]);
    
    #Extract the subnetworks
    net1=string_db[indxs1,c(1,2)];net2=string_db[indxs2,c(2,1)];
    colnames(net2)=c("main_protein","neigh_protein")
    colnames(net1)=c("main_protein","neigh_protein")
    #Assemble them
    net=rbind(net1,net2);rm(net1,net2);
    return(net)
  }
  
  #Close cores ----
  stopCluster(cl)
  
  return(res.l)
}