#########################################################################################################
# Author: Roula Tsonaka (s.tsonaka@lumc.nl)
#         Leiden University Medical Center (LUMC)
#########################################################################################################
# Title: SupportFunctions.R
#
# Purpose: several additional supporting functions to implement the 2-stage approach of Tsonaka et al. (2012) on a
#          longitudinal outcome.
# 
#
# Reference: Tsonaka, R., van der Helm - Mil, A., Houwing-Duistermaat, J. (2012). 
#             A two-stage mixed-effects model approach for gene-set analyses in candidate gene studies. 
#             Statistics in Medicine, 13, 1190 – 1202.
#
# Date: 05SEPTEMBER2016
###########################################################################################################




# Metropolis-Hastings
mh.b <- function(i, n.sims, n.gen, start.id, mu.id, sd.id, mu.gen, var.gen, df., beta.snp, start.gen, scale.id, scale.gen){
  dat.i <- data.snps[which(data.snps$"id." %in% i), ]
  Z2 <- model.matrix(~ gen.ind - 1, data = dat.i)
  dimnames(Z2) <- NULL
  bid.cur <- start.id
  bgen.cur <- start.gen
  draws. <- matrix(NA, nrow = n.sims, ncol = n.gen + 1)
  
  b.update <- function(bid.cur, bgen.cur) {
    argen <- rep(0, n.gen)
    arid <- 0  
    if(n.gen == 1){ offid <- beta.snp} else{ offid <- beta.snp + rowSums(Z2 * bgen.cur)}  
    if(sd.id > 0){
      bid.can <- rnorm(1, mean = mu.id, sd = scale.id * sd.id)/sqrt(rchisq(1, df.)/df.)                    
      eta.n <- offid + bid.can            
      eta.dn <- offid + bid.cur
      ex.n1 <- sum(dbinom(dat.i$snps, size = dat.i$size, prob = plogis(eta.n), log = TRUE)) + dnorm(bid.can, 0, sd.b.snps[2], log = TRUE) + dmvnorm(bgen.cur, rep(0, n.gen), diag(sd.b.snps[1]^2, 3), log = TRUE)
      ex.n2 <- dgt(bid.cur, mu = mu.id, sigma = scale.id * sd.id, df = df., log = TRUE)
      ex.dn1 <- sum(dbinom(dat.i$snps, size = dat.i$size, prob = plogis(eta.dn), log = TRUE)) + dnorm(bid.cur, 0, sd.b.snps[2], log = TRUE) + dmvnorm(bgen.cur, rep(0, n.gen), diag(sd.b.snps[1]^2, 3), log = TRUE)
      ex.dn2 <- dgt(bid.can, mu = mu.id, sigma = scale.id * sd.id, df = df., log = TRUE)
      accept.prob <- (ex.n1 - ex.n2 - ex.dn1 + ex.dn2)
      #cat("\nID", exp(accept.prob) )
      if (!is.na(accept.prob) & log(runif(1)) <= accept.prob){ 
        bid.cur <- bid.can; arid <- 1#;cat("\nAccept!")
        } else{ bid.cur <- bid.cur}
    }else{
      bid.can <- bid.cur
      
    }                        
    offgen <- beta.snp + bid.cur
    for(g in 1:n.gen){
      bgen.can <- bgen.cur
      bgen.can[g] <- rnorm(1, mean = mu.gen[g], sd = scale.gen * sqrt(var.gen[g]))/sqrt(rchisq(1, df.)/df.)
      #rmvt(n = 1, mu = mu.gen, Sigma = scale.gen * diag(var.gen), df = df.)
      eta.n <- offgen + rowSums(Z2 * bgen.can)
      eta.dn <- offgen + rowSums(Z2 * bgen.cur)
      ex.n1 <- sum(dbinom(dat.i$snps, size = dat.i$size, prob = plogis(eta.n), log = TRUE)) + dmvnorm(bgen.can, rep(0, n.gen), diag(sd.b.snps[1]^2, 3), log = TRUE) + dnorm(bid.cur, 0, sd.b.snps[2], log = TRUE)
      ex.n2 <- dgt(bgen.cur[g], mu = mu.gen[g], sigma = scale.gen * sqrt(var.gen[g]), df = df., log = TRUE)
      #ex.n2 <- dmvt(bgen.cur, mu = mu.gen, Sigma = scale.gen * diag(var.gen), df = df., log = TRUE)
      ex.dn1 <- sum(dbinom(dat.i$snps, size = dat.i$size, prob = plogis(eta.dn), log = TRUE)) + dmvnorm(bgen.cur, rep(0, n.gen), diag(sd.b.snps[1]^2, 3), log = TRUE) + dnorm(bid.cur, 0, sd.b.snps[2], log = TRUE)
      ex.dn2 <- dgt(bgen.can[g], mu = mu.gen[g], sigma = scale.gen * sqrt(var.gen[g]), df = df., log = TRUE)
      #ex.dn2 <- dmvt(bgen.can, mu = mu.gen, Sigma = scale.gen * diag(var.gen), df = df., log = TRUE)
      accept.prob <- (ex.n1 - ex.n2 - ex.dn1 + ex.dn2)
      #cat("\nGENE", exp(accept.prob) )
      if (!is.na(accept.prob) & log(runif(1)) <= accept.prob){
        bgen.cur[g] <- bgen.can[g]; argen[g] <- 1#;cat("\nAccept!")
        } else {bgen.cur[g] <- bgen.cur[g]}                     
    }
    #cat("\nID", i, "\tSample", s.)  
    c(bid.cur, bgen.cur, arid, argen)
  }
  
  for (s. in 1:n.sims) {
    dr <- b.update(start.id, start.gen)
    draws.[s., ] <- dr[1:4]
    if(s. == 1){ ar.id <- dr[5] }else{ ar.id <- ar.id + dr[5]}
    if(s. == 1){ ar.gen <- dr[6:8] }else{ ar.gen <- ar.gen + dr[6:8]}
  }
  #print(c(ar.id, ar.gen)/n.sims)
  return(draws.)
}


dmvt <- function (x, mu, Sigma, df, log = FALSE) {
    if (!is.numeric(x)) 
        stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
        x <- rbind(x)
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p)) || ncol(x) != p) 
        stop("incompatible arguments")
    ed <- eigen(Sigma, symmetric = TRUE)
    ev <- ed$values
    if (!all(ev >= -1e-06 * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
    ss <- x - rep(mu, each = nrow(x))
    inv.Sigma <- ed$vectors %*% (t(ed$vectors)/ev)
    quad <- rowSums((ss %*% inv.Sigma) * ss)/df
    fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + 
        log(df)) + sum(log(ev)))
    if (log) 
        fact - 0.5 * (df + p) * log(1 + quad)
    else exp(fact) * ((1 + quad)^(-(df + p)/2))
}

rmvt <- function (n, mu, Sigma, df){
    mvrnorm(n, mu = mu, Sigma = Sigma)/rep(sqrt(rchisq(n, df)/df), length(mu))
}

dgt <- function(x, mu = 0, sigma = 1, df = stop("no df arg"), log = FALSE){
    if (log)
        dt((x - mu)/sigma, df = df, log = TRUE) - log(sigma)
    else
         dt((x - mu)/sigma, df = df) / sigma
}




matFun <- function(lis, FUN){
  if(!is.list(lis) || !all(sapply(lis, is.matrix)))
    stop("'lis' must be a list containing 2-dimensional arrays")
  dims <- sapply(lis, dim)
  n <- dims[1, 1]
  p <- dims[2, 1]
  if(!all(n == dims[1, ]) || !all(p == dims[2, ]))
    stop("the matrices must have the same dimensions")
  mat <- matrix(unlist(lis), n * p, length(lis))
  matrix(apply(mat, 1, FUN), n, p)
}

matSums <- function(lis){
  if(!is.list(lis) || !all(sapply(lis, is.matrix)))
    stop("'lis' must be a list containing 2-dimensional arrays")
  dims <- sapply(lis, dim)
  n <- dims[1, 1]
  p <- dims[2, 1]
  if(!all(n == dims[1, ]) || !all(p == dims[2, ]))
    stop("the matrices must have the same dimensions")
  out <- array(data = 0, dim = c(n, p))
  for(i in seq(along = lis))
    out <- out + lis[[i]]
  out
}

matMeans <- function(lis)
  matSums(lis) / length(lis)

