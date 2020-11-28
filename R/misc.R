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

#' @export
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
#' @export
simulateTraits <- function(X,
                           q,
                           mu,
                           tau,
                           cis = NULL,
                           qtl = c(0, 1, 1),
                           qtl2 = c(0, 1, 1),
                           qtl3 = c(0, 1, 1),
                           mixedd = 'normal') {
  ## generate vector of q traits for a given design matrix (X) with polygenci effects
  ## and up to 3 true QTLs
  # mu, tau -- polygenic effect (mean and sd)
  # cis -- c(mean, sd)
  # qtl -- c(position: marker, trait; beta)
  M <- ncol(X)
  if(mixedd == "normal"){
    beta.matrix <- matrix(rnorm(q * M, mu, tau), ncol = q, byrow = TRUE)
  }else if(mixedd  == "Laplace"){
    beta.matrix <- matrix(rep(mu,M) + sqrt(rgamma(q*M, 1)) * rnorm(q * M, 0, tau), ncol = q, byrow = TRUE)
  }
  if (!is.null(cis)) {
    cis.matrix <- diag(rnorm(M, cis[1], cis[2]))
  } else cis.matrix <- 0
  qtl.matrix <- matrix(0, M, q)
  qtl.matrix[qtl[2], qtl[3:length(qtl)]] <- qtl[1]
  qtl2.matrix <- matrix(0, M, q)
  qtl2.matrix[qtl2[2], qtl2[3:length(qtl2)]] <- qtl2[1]
  qtl3.matrix <- matrix(0, M, q)
  qtl3.matrix[qtl3[2], qtl3[3:length(qtl3)]] <- qtl3[1]
  e.matrix <- replicate(q, rnorm(n, 0, sigma))
  y <- X %*% (beta.matrix + cis.matrix + qtl.matrix + qtl2.matrix + qtl3.matrix) + e.matrix
  return(y)
}

#' @export
plotTrait <- function(t, markers, crit, smooth = FALSE, addline = NULL) {
  trait <- data.frame(Marker = 1:length(t), t = t)
  if (smooth) trait$t <- smooth.spline(trait$t)$y
  p <- ggplot(trait) +
    geom_line(aes(Marker, t), size = 0.3) +
    geom_vline(xintercept = cumsum(markers), linetype = 2) +
    geom_hline(yintercept = crit) +
    scale_x_continuous(limits = c(0, sum(markers)), expand = c(0, 0), breaks = cumsum(markers)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),axis.text.x = element_text(hjust = 1))
  if (!is.null(addline)) p <- p + geom_line(aes(Marker, addline), col = "red", size = 0.3)
  p
}

#' @export
plotTraits <- function(t, markers, crit, addline = NULL) {
  library(reshape2)
  traits <- melt(as.data.frame(cbind(Marker = 1:ncol(t), t)),
                 id = "Marker", value.name = "t")
  p <- ggplot(traits) +
    geom_line(aes(x = Marker, y = t, group = variable), size = 0.3) +
    geom_vline(xintercept = cumsum(markers), linetype = 2) +
    geom_hline(yintercept = crit) +
    scale_x_continuous(expand = c(0, 0), breaks = cumsum(markers)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),axis.text.x = element_text(hjust = 1))
  if (!is.null(addline)){
    line <- data.frame(Marker = 1:ncol(t), addline = addline)
    p <- p + geom_line(data = line, aes(Marker, addline), col = "red", size = 0.3)
  }
  p
}

#' @export
plotTraits2 <- function(t,start_marker, crit, addline = NULL) {

  traits <- melt(as.data.frame(cbind(Marker = start_marker +(1:length(t)), t)),
                 id = "Marker", value.name = "t")
  p <- ggplot(traits) +
    geom_line(aes(x = Marker, y = t, group = variable), size = 0.3) +
    #geom_vline(xintercept = cumsum(markers), linetype = 2) +
    geom_hline(yintercept = crit) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),axis.text.x = element_text(hjust = 1))
  if (!is.null(addline)){
    line <- data.frame(Marker = 1:ncol(t), addline = addline)
    p <- p + geom_line(data = line, aes(Marker, addline), col = "red", size = 0.3)
  }
  p
}

### Figures
#' @export
plotHot <- function(t, crit, markers) {
  ind <- as.data.frame(arrayInd(which(abs(t) > crit), dim(t)))
  colnames(ind) <- c("Trait", "Marker")
  ggplot(ind) +
    geom_point(aes(Marker, Trait), col = "red", size = 0.3) +
    geom_vline(xintercept = cumsum(markers), linetype = 2) +
    scale_x_continuous(limits = c(0, sum(markers)), expand = c(0, 0), breaks = cumsum(markers)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),axis.text.x = element_text(hjust = 1))
}

#' @export
replaceNA <- function(X) {
  ## replace missing values with the mean
  ## (function uses replaceOneNA)
  na_pos <- apply(X, 1, function(x) which(is.na(x)))
  for (i in seq_along(na_pos)) {
    if (length(na_pos[[i]]) == 0) next
    X[i, na_pos[[i]]] <- sapply(na_pos[[i]], replaceOneNA, x = X[i, ])
  }
  return(X)
}

replaceOneNA <- function(na, x) {
  ## replace one missing value
  if (!is.na(x[na])) stop("Value is not NA.")
  n <- length(x)
  left <- which(!is.na(x[na:1]))[1] - 1
  right <- which(!is.na(x[na:n]))[1] - 1
  if (is.na(left) & is.na(right)) stop("All values are NA.")
  if (is.na(left)) left <- Inf else if (is.na(right)) right <- Inf
  if (left < right) {
    replace <- x[na - left]
  } else if (left > right) {
    replace <- x[na + right]
  } else {
    replace <- (x[na - left] + x[na + right])/2
  }
  return(replace)
}

#' @export
pseudovalue <- function(m1, m2, d1, d2){
  ## calculate expected value of marker m between m1 and m2 (conditional)
  ## for m in {-1, 0, 1}
  r1 <- (1 - exp(-2*d1/100)) / 2
  r2 <- (1 - exp(-2*d2/100)) / 2
  r <- r1 + r2 - 2*r1*r2
  num <- (1 - r1)^(1 + m1) * r1^(1 - m1) * (1 - r2)^(1 + m2) * r2^(1 - m2) -
    (1 - r1)^(1 - m1) * r1^(1 + m1) * (1 - r2)^(1 - m2) * r2^(1 + m2)
  den <- r^abs(m1 - m2) * (1 - r)^(2 - abs(m1 - m2))
  return(num/den)
}

#' @export
pseudomarkers <- function(X, markers, len_cum, by = 2){
  ## put pseudomarkers between the first and last markers
  # markers -- vector, number of markers on each chromosome
  # len_cum -- distance from the begining of chromosome
  n <- nrow(X)
  chr <- length(markers)
  markers_cum <- c(0, cumsum(markers))
  # first marker on a given chromosome must be at the distance of 0:
  len_cum <- len_cum - rep(len_cum[markers_cum[-(chr + 1)] + 1], markers)
  chr_len <- len_cum[cumsum(markers)] # lengths of chromosomes
  markers <- c(0, cumsum(markers))
  new_markers <- numeric(chr)
  XX <- NULL
  l <- 1
  for(i in 1:chr) {
    new_markers[i] <- ceiling(chr_len[i] / by)
    XX <- cbind(XX, matrix(0, n, new_markers[i]))
    for(j in 1:new_markers[i]) {
      m <- sum(len_cum[(markers[i] + 1) : markers[i + 1]] <= (j - 1)*by) + markers[i]
      m1 <- X[, m]
      m2 <- X[, m + 1]
      d1 <- (j - 1)*by - len_cum[m]
      d2 <- len_cum[m + 1] - (j - 1)*by
      XX[, l] <- pseudovalue(m1, m2, d1, d2)
      l <- l + 1
    }
  }
  return(list(X = XX, markers = new_markers))
}
#' @export
CIMModel <- function(X, y, window = 10) {
  ## composite interval mapping
  M <- ncol(X)
  n <- nrow(X)
  t <- numeric(M)
  for (j in 1:M) {
    Xm <- matrix(1, n, 2)
    Xm[, 2] <- X[, j]
    add <- setdiff(1:M, (j - window):(j + window))
    Xm <- cbind(Xm, X[, add])
    fit <- fastLmPure(Xm, y)
    t[j] <- fit$coef[2] / fit$se[2]
  }
  t[!is.finite(t)] <- 0
  return(t)
}

