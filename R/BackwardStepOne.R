#' @useDynLib bwd backward
"BackwardStepOne" <- function(y, lastkgroup = floor(0.01*n), kmin = 5) {
n <- length(y)
m <- y
z <- rep(1, n)
rss <- z[1:(n-1)] * z[2:n] / (z[1:(n-1)] + z[2:n]) * (m[1:(n-1)] - m[2:n])^2
l <- n
value <- rep(0, n)
index <- rep(0, n)
ngroup <- lastkgroup
sigma <- estimateSigma(y, h = 10)
cutoff <- 1.0e+12

temp <- .C("backward", as.double(m),
                       as.integer(z),
                       as.double(rss),
                       as.integer(n),
                       as.integer(l),
                       as.double(value),
                       as.integer(index),
                       as.double(sigma),
                       as.double(cutoff),
                       as.integer(ngroup),
                       as.integer(kmin))
        l <- temp[[5]]
        m <- temp[[1]][1:l]
        z <- temp[[2]][1:l]
        rss <- c(temp[[3]][1:(l-1)], NA)
        index <- c(1, cumsum(z)[-l]+1)
obj <- cbind(index, m, z, rss)
colnames(obj) <- c("index", "mean", "size", "rss")
return(obj)
}
