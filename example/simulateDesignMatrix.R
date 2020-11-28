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

calculateCorrelations <- function(t, markers, k = 5) {
  ## calculate correlations between t statistics
  # t -- vector od t statistics
  # markers -- markers on each chromosome
  # k -- correlations are calculated between (t_i, t_(i-1)), ..., (t_i, t_(i-k))
  M <- sum(markers) - length(markers)
  cor.t <- numeric(k)
  cumm <- cumsum(markers)
  for (i in 1:k) {
    t1 <- t[-cumm][1:(M-i)]
    t2 <- t[-(cumm + i)][(i+1):M]
    cor.t[i] <- cor(t1, t2)
  }
  return(cor.t)
}

calculateCrit <- function(t, markers, d = 1, alpha = 0.05) {
  ## calculate critical value according to Biometrics correcion
  # t -- vector od t statistics
  # markers -- markers on each chromosome
  # d -- distance between markers (cM)
  # alpha -- significance levele

  C <- length(markers)
  L <- sum(markers)
  cor.t <- calculateCorrelations(t, markers)
  x <- 1:max(sum(cor.t > 0.3), 2)
  cor.t <- cor.t[x]
  (beta <- -coef(lm(log(cor.t) ~ x - 1)))

  v <- function(y) (2 / y * (pnorm(y / 2) - 0.5)) / (y / 2 * pnorm(y / 2) + dnorm(y / 2))
  Pr <- function(z) {
    1 - exp(-2*C*(1 - pnorm(z)) - 2 * beta * L * z * dnorm(z) * v(z * sqrt(2 * beta * d))) -
      alpha
  }
  crit <- uniroot(Pr, c(1, 10))$root
  return(crit)
}
