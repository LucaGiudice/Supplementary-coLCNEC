#Watch the difference between raw and filtered data
plot_comp_rawVSfilt = function(x,lcpm,M,L,samplenames){
  lcpm.cutoff <- log2(10/M + 2/L)
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", samplenames, text.col=col, bty="n")
  lcpm <- cpm(x, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", samplenames, text.col=col, bty="n")
}

compl_m = function(input_m,th=0.70){
  cat(dim(input_m),"\n")
  #Remove rows with 70% of zeros ----
  th=round(ncol(input_m)*th)
  keep_row=0
  for(k in 1:nrow(input_m)){
    x=input_m[k,]
    n_zeros=sum(x!=0)
    if(n_zeros>=th){
      keep_row=c(keep_row,k)
    }
  }
  keep_row=keep_row[-1]
  input_m=input_m[keep_row,]
  input_m[input_m==0]=NA
  cat(dim(input_m),"\n")
  return(input_m)
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

unspe_gsK_detection=function(gs_ints_l,k_comb=2){
  require("Rfast")
  gs2_ints=c("START")
  combs=Rfast::comb_n(n=seq(1,4), k=k_comb,simplify=TRUE)
  for(comb in 1:ncol(combs)){
    indxs=combs[,comb]
    gs2_ints=c(gs2_ints,Reduce(intersect, gs_ints_l[indxs]))
  }
  gs2_ints=unique(gs2_ints[-1])
  return(gs2_ints)
}

prep_gs_ints_list=function(DE_genes_subt1,GSEA_enr_genes_l_subt1,jump=4){
  gs_ints_l=list()
  
  for(k in 1:length(DE_genes_subt1)){
    n_genes=length(DE_genes_subt1[[k]])
    if(n_genes<50){
      DE_genes_subt1[[k]]=GSEA_enr_genes_l_subt1[[k]]
    }
  }
  
  for(i in 1:round(length(DE_genes_subt1)/2)){
    j=i+jump
    gs_ints=intersect(DE_genes_subt1[[i]],DE_genes_subt1[[j]])
    gs_ints_l[[i]]=gs_ints
  }
  
  return(gs_ints_l)
}