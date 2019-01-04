get.cutoff <- function(n, alpha = 0.05)
{
  b <- 0.2411
  alpha.list <- c(0.1, 0.05, 0.025, 0.01, 0.005)
  a.list <- c(2.2135, 2.5049, 2.7434, 3.0121, 3.3869)
  
  id <- match(alpha, alpha.list)
  a <- a.list[id]
  
  if (all(is.na(a))) {
    stop("alpha should be 0.1, 0.05, 0.025, 0.01 or 0.005!")
  } else if (all(is.na(a))) {
    warning("alpha should be 0.1, 0.05, 0.025, 0.01 or 0.005!")
    a <- a[!is.na(a)]
  } else a <- a
  value <- a + b * log(n)
  return(value)
}

