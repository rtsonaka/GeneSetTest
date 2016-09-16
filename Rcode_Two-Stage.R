#########################################################################################################
# Author: Roula Tsonaka (s.tsonaka@lumc.nl)
#         Leiden University Medical Center (LUMC)
#########################################################################################################
# Title: Rcode_TwoStage.R
#
# Aim: R code to run the 2-stage approach of Tsonaka et al. (2012).
#      
# Reference: Tsonaka, R., van der Helm - Mil, A., Houwing-Duistermaat, J. (2012). 
#             A two-stage mixed-effects model approach for gene-set analyses in candidate gene studies. 
#             Statistics in Medicine, 13, 1190 – 1202.
#
# Date: 05SEPTEMBER2016
###########################################################################################################


# load data: snps and longitudinal responses
data.snps <- read.table("C:/RoulaTsonaka/LeidenRA/Paper2/LaTex Files/Revision1/GIT/snps.txt")
# 150 snps have been simulated in 3 genes: 20 in gene1, 40 in gene2 and 90 in gene3.

data.y <- read.table("C:/RoulaTsonaka/LeidenRA/Paper2/LaTex Files/Revision1/GIT/data.txt")
# longitudinal data are simulated with SNPs c(10, 40, 105) directly associated with progression



# 1st stage: mixed-effects logistic regression
library(lme4)
library(aod)
library(mvtnorm)

model.snps0 <- try(glmer(cbind(snps, size - snps) ~ 1 + (1 | id./gen.ind), family = binomial, data = data.snps))
n.gen <- length(unique(data.snps$gen.ind))

if(!inherits(model.snps0, "try-error")){    
  beta.snp <- fixef(model.snps0)
  sd.b.snps <- unlist(lapply(VarCorr(model.snps0), function(x) sqrt(x[1])))
  vcov.s <-  matrix(vcov(model.snps0)@ "x", byrow = TRUE, ncol = length(fixef(model.snps0)))
  rr <- ranef(model.snps0, condVar = TRUE)
  eb.hat <- matrix((rr$"gen.ind:id.")[, "(Intercept)"], byrow = FALSE, ncol = n.gen) + rr$"id."[, "(Intercept)"]
  post.varb.id <- unlist(lapply(attr(rr$"id.", "postVar"), function(x) x))
  post.varb.gen <- matrix(unlist(lapply(attr(rr$"gen.ind:id.", "postVar"), function(x) x)), byrow = FALSE, ncol = n.gen)
}else{ 
  cat("\nSNP Model failed")
  }

# Store EB estimates per gene
ind. <- match(data.y$"idy", unique(data.y$"idy"))
data.y$"eb1" <- eb.hat[ind., 1] # EB estimates of RE for gene1
data.y$"eb2" <- eb.hat[ind., 2] # EB estimates of RE for gene2
data.y$"eb3" <- eb.hat[ind., 3] # EB estimates of RE for gene3

# 2nd stage: test gene effect on progression
fm1 <- as.formula(paste("Resp ~ Time.y * (eb1 + eb2 + eb3) + (Time.y|idy)"))
model.1 <- try(lmer(fm1, data = data.y, REML = FALSE))
fm11 <- as.formula(paste("Resp ~ Time.y + (eb1 + eb2 + eb3) + (Time.y|idy)"))
model.11 <- try(lmer(fm11, data = data.y, REML = FALSE))

if(!inherits(model.1, "try-error") && !inherits(model.11, "try-error")){ 
  beta.y <- by. <- fixef(model.1)
  vcov.by <-  vb.. <- matrix(vcov(model.1)@ "x", byrow = TRUE, ncol = length(fixef(model.1)))
  vby <- as.data.frame(VarCorr(model.1))

  if(ncol(ranef(model.1)$idy) == 1){
    sd.b.y <- vby[1, "sdcor"]
    }else{
      sd.b.y <- vby[, "sdcor"]
      }             

  # Wald test per gene and as gene-set    
  nam.t1 <- names(fixef(model.1))   
  pb <- (anova(model.1, model.11)$"Pr(>Chisq)")[2]
  terms.1 <- which(nam.t1 %in% c("Time.y:eb1"))
  p1w <- wald.test(b = by., Sigma = vb.., Terms = terms.1)$result$chi2["P"]
  terms.2 <- which(nam.t1 %in% c("Time.y:eb2"))
  p2w <- wald.test(b = by., Sigma = vb.., Terms = terms.2)$result$chi2["P"]
  terms.3 <- which(nam.t1 %in% c("Time.y:eb3"))
  p3w <- wald.test(b = by., Sigma = vb.., Terms = terms.3)$result$chi2["P"]
  p4w <- wald.test(b = by., Sigma = vb.., Terms = c(terms.1, terms.2, terms.3))$result$chi2["P"]
      

    
  # Metropolis-Hastings to account variability of EB estimates
  B <- 10
  df.t <- 4
  thetas <- matrix(NA, nrow = B, ncol = length(fixef(model.1)))
  vthetas <- vector(length = B, mode = "list")
  id <- data.snps$id.
  idn <- match(id, unique(id))
  n <- length(unique(id))
  data.new <- data.y
    
  res1 <- res2 <- res3 <- res4 <- matrix(NA, ncol = B, nrow = n)
  for(j in 1:n){
    id.indxx <- j                       
    start.id <- rr$"id."[id.indxx, "(Intercept)"] 
    start.gen <- matrix((rr$"gen.ind:id.")[, "(Intercept)"], byrow = FALSE, ncol = n.gen)[id.indxx, ]
    var.id <- post.varb.id[id.indxx] 
    var.gen <- post.varb.gen[id.indxx, ]  
    mh.draws <- mh.b(id.indxx, n.sims = B, n.gen = n.gen, start.id = start.id, mu.id = start.id, sd.id = sqrt(var.id), 
                       mu.gen = start.gen, var.gen = var.gen, df. = df.t, 
                       beta.snp = beta.snp, start.gen = start.gen, scale.id = 1, scale.gen = 1)
    res1[j, ] <- mh.draws[, 1]
    res2[j, ] <- mh.draws[, 2]
    res3[j, ] <- mh.draws[, 3]
    res4[j, ] <- mh.draws[, 4]
    #print(j)
  }
    
  for (b in seq_len(B)){
    data.new <- data.y    
    ebb.hat <- res1[, b] + cbind(res2[, b], res3[, b], res4[, b])
    data.new$"ebb1" <- ebb.hat[ind., 1]
    data.new$"ebb2" <- ebb.hat[ind., 2]
    data.new$"ebb3" <- ebb.hat[ind., 3]
    fm1 <- as.formula(paste("Resp ~ Time.y * (ebb1 + ebb2 + ebb3) + (Time.y|idy)"))
    model.1.n <- try(lmer(fm1, data = data.new, REML = FALSE))
    if(!inherits(model.1.n, "try-error")){
      thetas[b, ] <- fixef(model.1.n)
      vthetas[[b]] <- matrix(vcov(model.1.n)@ "x", byrow = TRUE, ncol = length(fixef(model.1.n)))
      }
      cat("\nSampling:", b)
  }
    
  V1 <- ((B + 1)/B) * var(thetas, na.rm =  TRUE)
  V2 <- matMeans(vthetas)
  V <- V1 + V2
  vcov.boot <- V                    
  terms.1 <- which(names(fixef(model.1.n)) %in% c("Time.y:ebb1"))
  p1b. <- wald.test(b = colMeans(thetas, na.rm = TRUE), Sigma = V, Terms = terms.1)$result$chi2["P"]
  p1b <- wald.test(b = by., Sigma = V, Terms = terms.1)$result$chi2["P"]
  terms.2 <- which(names(fixef(model.1.n)) %in% c("Time.y:ebb2"))
  p2b. <- wald.test(b = colMeans(thetas, na.rm = TRUE), Sigma = V, Terms = terms.2)$result$chi2["P"]
  p2b <- wald.test(b = by., Sigma = V, Terms = terms.2)$result$chi2["P"]
  terms.3 <- which(names(fixef(model.1.n)) %in% c("Time.y:ebb3"))
  p3b. <- wald.test(b = colMeans(thetas, na.rm = TRUE), Sigma = V, Terms = terms.3)$result$chi2["P"]
  p4b. <- wald.test(b = colMeans(thetas, na.rm = TRUE), Sigma = V, Terms = c(terms.1, terms.2, terms.3))$result$chi2["P"]         
  p3b <- wald.test(b = by., Sigma = V, Terms = terms.3)$result$chi2["P"]
  p4b <- wald.test(b = by., Sigma = V, Terms = c(terms.1, terms.2, terms.3))$result$chi2["P"]         
  
  
  pvals <- c(pb, p1w, p2w, p3w, p4w, p1b, p2b, p3b, p4b)
  names(pvals) <- c("LRT-plugin", "W1-plugin", "W2-plugin", "W3-plugin", "MW-plugin", 
                    "W1-MH", "W2-MH", "W3-MH", "MW-MH")
  #pvalues: 
  #"LRT-plugin", "W1-plugin", "W2-plugin", "W3-plugin", "MW-plugin": 
  #LRT and wald tests first per gene and then multivariate wald. We do not account for variability of EB.
  #"W1-MH", "W2-MH", "W3-MH", "W4-MH", "MW-MH":
  #Wald tests first per gene and then multivariate wald. We account for variability of EB.
  
}else{ 
  cat("\nModel failed")
  }# model fit error


