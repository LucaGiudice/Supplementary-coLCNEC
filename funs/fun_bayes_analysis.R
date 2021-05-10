read.peptides <- function(dat, cha){
  output <- NULL
  
  dat$Sequence <- as.character(dat$Sequence)
  dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
  dat$Quan.Usage <- as.character(dat$Quan.Usage)
  dat$Quan.Info <- as.character(dat$Quan.Info)
  dat$Isolation.Interference <- as.numeric(as.character(dat$Isolation.Interference))
  
  dat <- subset(dat, Isolation.Interference<=30)  
  dat <- subset(dat, Quan.Usage=="Used")
  dat <- subset(dat, Protein.Group.Accessions!="")
  dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))  
}

quantify.proteins <- function(dat, cha){
  e.function <- function(x, seq) tapply(x, seq, median)
  output <- NULL
  
  dat$Sequence <- toupper(dat$Sequence) # Capital letters
  accessions <- as.character(unique(dat$Protein.Group.Accessions))
  n.proteins <- length(accessions)
  n.cha <- length(cha)
  
  for(k in 1:n.proteins){
    id <- accessions[k]
    sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
    sdat[cha] <- log2(sdat[cha])
    sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
    pdat <- sdat[, -1]
    n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
    temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])          
    n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)    
    if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
    pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
    output <- rbind(output, pdat)
  }
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- round(output[,1:n.cha],3)
  row.names(output) <- accessions
  output <- as.data.frame(output)
  return(output)
}

eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

eb.fit.mult <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coef[, "tr2"]
  df.0 <- rep(fit.eb$df.prior, n)
  df.r <- fit.eb$df.residual
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
  t.mod <- fit.eb$t[, "tr2"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "tr2"]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb.mult <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.0, df.r, s2.0, s2, s2.post)
  results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
  return(results.eb.mult)
}