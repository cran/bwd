"BackwardStepTwo" <- function(obj, y = y, kmin = 2, min.group = 3, h = 10, mu0 = 0) {
  group <- obj[,1] # starting index
  means <- or.means <- obj[,2] # group means
  gSize <- obj[,3] # group size
  rss <- obj[,4]   # rss

  nn <- nrow(obj)
  hat.sigma <- estimateSigma(y, h)

  # epidemic change
  if (is.numeric(mu0) & length(mu0) == 1) {
    z <- sqrt(gSize) * (means - mu0)/hat.sigma
    mu0.id <- which(abs(z) < qnorm(1 - 0.05/2))
    means[mu0.id] <- mu0
    rss <- rss + c(diff(gSize * (or.means - means)^2), NA)
  }

  result <- as.list(1:(nn - min.group))
  result[[1]] <- data.frame(group = group, means = means, size = gSize)[-1,]
  di <- rep(0, nn-1)
  two.index <- two.size <- matrix(0, 2, nn-1)
  for (m in 1:(nn-min.group-1))
  {
    i <- which.min(rss)
    DI <- 0
    two.size[,m] <- both.size <- gSize[i:(i+1)]
    two.index[,m] <- both.index <- group[i:(i+1)]
    if (min(both.size) > kmin) {
      temp.group <- group[i:(i+1)]
      temp.means <- means[i:(i+1)]
      temp.gSize <- gSize[i:(i+1)]
      DI <- abs(temp.means[1] - temp.means[2])/(sqrt(1/temp.gSize[1] + 1/temp.gSize[2]) * hat.sigma)
    }
    di[m] <- DI
    group <- group[-(i+1)]

    merge.gSize <- gSize[c(i,i+1)]
    gSize[i] <- sum(merge.gSize)
    gSize <- gSize[-(i+1)]

    means[i] <- (merge.gSize[1] * means[i] + merge.gSize[2] * means[i+1])/gSize[i]

    # Epidemic change
    if (is.numeric(mu0) & length(mu0) == 1) {
      z <- sqrt(gSize[i]) * (means[i] - mu0)/hat.sigma
      shrinkage  <- abs(z) < qnorm(1 - 0.05/2)
      if (shrinkage) means[i] <- mu0
    }

    means <- means[-(i+1)]

    if (i == (length(rss)-1)) {
      temp.id <- i:(i+1)
      temp.gSize <- gSize[temp.id]
      temp.means <- means[temp.id]
      temp1 <- temp.gSize * temp.means
      temp2 <- temp1 * temp.means
      temp.rss <- sum(temp2) - (sum(temp1))^2/sum(temp.gSize)
      rss[i] <- temp.rss
    } else if (i == 1){
      temp.id <- 1:2
      temp.gSize <- gSize[temp.id]
      temp.means <- means[temp.id]
      temp1 <- temp.gSize * temp.means
      temp2 <- temp1 * temp.means
      temp.rss <- sum(temp2) - (sum(temp1))^2/sum(temp.gSize)
      rss[2] <- temp.rss
    } else {
      temp.id <- (i-1):(i+1)
      temp.gSize <- gSize[temp.id]
      temp.means <- means[temp.id]
      temp1 <- temp.gSize * temp.means
      temp2 <- temp1 * temp.means

      temp.rss <- (temp2[-1] + temp2[-3]) - (temp1[-1] + temp1[-3])^2/(temp.gSize[-1] + temp.gSize[-3])

      rss[i-1] <- temp.rss[1]
      rss[i+1] <- temp.rss[2]
    }
    rss <- rss[-i]
    result[[m+1]] <- data.frame(group = group, means = means, size = gSize)[-1,]
  }

  obj <- list(result = result, DI = di, index = two.index, size = two.size, sigma = hat.sigma)
  return(obj)
}
