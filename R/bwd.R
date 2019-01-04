#' Backward procedure for the change point detection
#' @description
#' Implements backward procedure for detecting single or multiple change points.
#'
#' @param y observed data
#' @param alpha target level that detemines stopping criterion. Default is 0.05
#' @param kmin minimum length of segements for checking possible change points
#' @param lastkgroup We can abvoid chekcing possible change points when we have less groups than "lastkgroup" to improve computational efficiency. Default is 0.01 * n
#' @param mu0 Baseline mean value whe detecting epidemic chang points. Defalut is \code{NULL}
#' @param normal if \code{TRUE} normal cutoff values are used, and if \code{FALSE} residual permuted cutoff values are used. Default is \code{TRUE}
#' @param n.permute number of permutation when computing the permuted cutoff. Defalut is 1000
#' @param h bandwidth size for variance esitimator
#' @return bwd object that contains information of detected segments and significance levels
#' @author Seung Jun Shin, Yicaho Wu, Ning Hao
#' @references Shin, Wu, and Hao (2018+) A backward procedure for change-point detection with applications to copy number variation detection, arXiv:1812.10107.
#' @seealso \code{\link{plot.bwd}}
#' @examples
#' # simulated data
#' set.seed(1)
#' n <- 1000
#' L <- 10
#'
#' mu0 <- -0.5
#'
#' mu <- rep(mu0, n)
#' mu[(n/2 + 1):(n/2 + L)] <- mu0 + 1.6
#' mu[(n/4 + 1):(n/4 + L)] <- mu0 - 1.6
#' y <- mu + rnorm(n)
#' alpha <- c(0.01, 0.05)
#'
#' # BWD
#'obj1 <- bwd(y, alpha = alpha)
#'
#' # Modified for epidemic changes with a known basline mean, mu0.
#'obj2 <- bwd(y, alpha = alpha, mu0 = 0)
#'
#'par(mfrow = c(2,1))
#'plot(obj1, y)
#'plot(obj2, y)
#'
#' @export
#' @importFrom stats qnorm quantile
#'
bwd <- function(y, alpha = 0.05, kmin = 3, lastkgroup = floor(0.01*n), mu0 = NULL, normal = T, n.permute = 1000, h = 10) {
  n <- length(y)

  if (normal) {
    cutoff <-  get.cutoff(n, alpha)
  } else
  {
    if (n > 50000) warning("sample size is large. It may takes long to get cut off values without normality assumption!")
    r <- estimateResidual(y, h = h)
    cutoffs <- NULL
    for (iter in 1:n.permute)
    {
      set.seed(iter)
      index <- sample(n)
      pre.obj <- BackwardStepOne(r[index], lastkgroup = lastkgroup, kmin = kmin)
      obj <- BackwardStepTwo(pre.obj, r[index], kmin, min.group = 0, h, mu0)
      cutoffs[iter] <- max(obj$DI)
    }
    cutoff <-  quantile(cutoffs, 1 - alpha)
  }

  obj1 <- BackwardStepOne(y, lastkgroup, kmin)
  obj2 <- BackwardStepTwo(obj1, y, kmin, min.group = 0, h, mu0)
  hat.sigma <- estimateSigma(y, h)

  n.alpha <- length(alpha)
  segment <- as.list(1:n.alpha)
  for (k in 1:n.alpha)
  {
    sel <- which(obj2$DI > cutoff[k])

    if (length(sel) > 0)
    {
      obj3 <- obj2$result[[min(sel)]]

      id1 <- c(1, obj3[,1])
      id2 <- c(obj3[,1] - 1, n)
      m1 <- mean(y[id1[1]:id2[1]])
      means <- c(m1, obj3[,2])

      if (is.numeric(mu0) & length(mu0) == 1) {
        z <- sqrt(id2[1]) * abs(m1 - mu0)/hat.sigma
        if (z < qnorm(1-0.05/2)) means <- c(mu0, obj3[,2])
      }
      size <- id2 - id1 + 1
      obj <- data.frame(from = id1, to = id2, mean = means, size = size)
    } else {
      obj <- data.frame(from = 1, to = n, mean = mean(y), size = n)
    }

    segment[[k]] <- obj
  }
  result <- list(segment = segment, alpha = alpha)
  class(result) <- "bwd"
  return(result)
}
