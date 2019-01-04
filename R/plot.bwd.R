#' plot for the backward procedure for the change point detection
#' @description
#' A plot of segments estimated by the backward procedure.
#'
#' @param x bwd object
#' @param y observed data
#' @param ... graphical parameters
#' @return plot of estimated segments
#' @author Seung Jun Shin, Yicaho Wu, Ning Hao
#' @references Shin, Wu, and Hao (2018+) A backward procedure for change-point detection with applications to copy number variation detection, arXiv:1812.10107.
#' @seealso \code{\link{bwd}}
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
#' @export plot.bwd
#' @importFrom graphics lines plot legend
plot.bwd <- function(x, y, ...){
  Segment <- x$segment

  plot(y, ...)
  K <- length(Segment)
  mu <- matrix(0, length(y), K)
  for (k in 1:K){
    segment <- Segment[[k]]
    for (j in 1:nrow(segment))
    {
      mu[segment[j,1]:segment[j,2], k] <- segment[j,3]
    }
    lines(mu[,k], col = k + 1, lty = k)
  }
  legend("bottomright", paste(x$alpha), col = 1:K + 1, lty = 1:K)
}

