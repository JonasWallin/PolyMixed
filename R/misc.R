##
#  misc files
#
# D: 2019-03-27
##

#
#' @export
findDuplicate <- function(X) {
  ## find which columns are duplicated
  dupl.ind <- which(duplicated(X, MARGIN = 2)) # these columns are duplicated...
  dupl.repl <- which(duplicated(X, MARGIN = 2, fromLast = TRUE)) # ...with these columns
  dupl <- cbind(dupl.ind, dupl.repl)
  return(dupl)
}

calculateCorrelations <- function(t, markers, k = 5) {
  ## calculate correlations between t statistics
  # t -- vector od t statistics
  # markers -- markers on each chromosome
  # k -- correlations are calculated between (t_i, t_(i-1)), ..., (t_i, t_(i-k))
  if(length(markers)==1){
    cor.t = acf(t,lag=k ,plot = F,na.action = na.omit)$acf[2:(k+1)]
  }else{
    M <- sum(markers) - length(markers)
    cor.t <- numeric(k)
    cumm <- cumsum(markers)
    for (i in 1:k) {
      t1 <- t[-cumm][1:(M-i)]
      t2 <- t[-(cumm + i)][(i+1):M]
      cor.t[i] <- cor(t1, t2)
    }
  }
  return(cor.t)
}

calculateCrit <- function(t, markers, d = 1, alpha = 0.05, beta=NULL) {
  ## calculate critical value according to Biometrics correcion
  # t -- vector od t statistics
  # markers -- markers on each chromosome
  # d -- distance between markers (cM)
  # alpha -- significance levele

  C <- length(markers)
  L <- d*sum(markers)
  cor.t <- calculateCorrelations(t, markers)
  x <-  d* (1:max(sum(cor.t > 0.3), 2))
  cor.t <- cor.t[x]
  if(cor.t[1] <= 0)
    return(qnorm((1-alpha/2)^(1/L)))
  cor.t[cor.t<= 0] = 10^-14
  if(is.null(beta))
    (beta <- -coef(lm(log(cor.t) ~ x - 1)))

  v <- function(y) (2 / y * (pnorm(y / 2) - 0.5)) / (y / 2 * pnorm(y / 2) + dnorm(y / 2))
  Pr <- function(z) {
    1 - exp(-2*C*(1 - pnorm(z)) - 2 * beta * L * z * dnorm(z) * v(z * sqrt(2 * beta * d))) -
      alpha
  }
  crit <- uniroot(Pr, c(1, 10))$root
  return(crit)
}


#' simulateDesignMatrix
#' @export
simulateDesignMatrix <- function(chr, m, n, d, seed = NA) {
  ## generate a design matrix
  # chr -- number of chromosomes
  # m -- number of markers on each chromosome
  # n -- number of observations
  # d -- distance between markers (cM)
  # seed -- seed for a random number generator
  if (length(m) == 1) m <- rep(m, chr)
  m_cum <- c(0, cumsum(m))
  if(is.na(seed)==F)
    set.seed(seed)
  r <- (1 - exp(-2*d/100)) / 2
  X <- matrix(0, n, sum(m))
  for(i in 1:chr) {
    x <- 2 * sample(2, n, rep = TRUE) - 3 # szybsze od c(-1, 1)
    x <- matrix(rep(x, m[i]), ncol = m[i])
    z <- 2 * sample(2, n*m[i], rep = TRUE, prob = c(r, 1 - r)) - 3
    z <- cumprod(z)
    z <- t(matrix(z, ncol = n))
    X[, (m_cum[i]+1):m_cum[i+1]] <- (x*z + 1)/2
  }
  return(2 * X - 1) # -1/1 coding
}

#' calculate power, FP, FDR, FWER, distance and beta error
#' @param find      - (q x 1) list of position discovered
#' @param beta      - (q x 1) list of the coeffient of the discoveries
#' @param beta.true - (m x 1) vector of the true coeffients
#' @param qtl.pos   - (m x 1) vecotr of the true positions
#' @param maxCM     - (int) minimum abs dist coutning as a find
#' @export
summarizeFind.new <- function(find, beta, beta.true, qtl.pos, maxCM = 15){


  m <- length(qtl.pos)
  q <- length(find)
  find.matrix <- matrix(0,  q, m + 1)
  if(m > 0){
    dist.matrix <- matrix(NA, q, m)
    beta.matrix <- matrix(NA, q, m)
  }
  for (i in 1:q) {
    if (length(find[[i]]) == 0)
      next
    if(m == 0){
      find.matrix[i] <- length(find[[i]])
      next
    }
    ord <- order(find[[i]])
    find[[i]] <- find[[i]][ord]
    beta[[i]] <- beta[[i]][ord]
    dist <- sapply(find[[i]], function(x) min(abs(x - qtl.pos)))
    k <- length(dist)
    b <- length(beta[[i]])
    if (k) beta.qtls <- beta[[i]][1:k]
    close.qtls <- sapply(find[[i]], function(x) which.min(abs(x - qtl.pos)))
    close.qtls[dist > maxCM] <- m + 1

    if (length(close.qtls)) {
      tab <- table(close.qtls)
      names <- as.numeric(names(tab))
      find.matrix[i, names] <- tab
      if (sum(close.qtls != m + 1)) {
        names <- names[names != m + 1]
        tab <- tab[names(tab) != m + 1]
        dist.qtls <- dist[dist <= maxCM]
        dist.qtls <- tapply(dist.qtls, rep(names, tab), mean)
        beta.qtls <- beta.qtls[dist <= maxCM]
        beta.qtls <- tapply(beta.qtls, rep(names, tab), mean)
        dist.matrix[i, names] <- dist.qtls
        beta.matrix[i, names] <- beta.qtls
      }
    }


  }
  fp       <- mean(find.matrix[, (m+1), drop=FALSE])
  fdr      <- mean(find.matrix[, (m+1), drop=FALSE] / apply(find.matrix, 1, sum), na.rm = T)
  fwer     <- mean(find.matrix[, (m+1)]>0)

  if(m>0){
    dist     <- apply(dist.matrix, 2, mean, na.rm = TRUE)
    power    <- apply(find.matrix[, -(m+1), drop=FALSE], 2, function(x) mean(x > 0))
    beta.dif <- abs(sweep(beta.matrix, 2, beta.true))
    beta.err <- apply(beta.dif, 2, mean, na.rm = TRUE)# / abs(beta.true)
  }

  if(m > 0){
    results <- matrix(NA, 3, m + 3)
    colnames(results) <- c(paste('QTL',c(1:m),sep=''), "FP", "FDR", "FWER")
    rownames(results) <- c("Power", "Dist", "Beta_err")
    results[1, ] <- c(power, fp, fdr, fwer)
    results[2, 1:m] <- dist
    results[3, 1:m] <- beta.err
  }else{
    results <- matrix(NA, 1, 3)
    colnames(results) <- c("FP", "FDR", "FWER")
    rownames(results) <- c("Power")
    results[1, ] <- c(fp, fdr, fwer)
  }

  return(results)
}

