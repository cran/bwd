#' @importFrom stats qt sd
estimateSigma <-
function(Y,h=10){
  n     = length(Y)
  YBar  = rep(0,n)
  for (i in 1:n) {
     a       = min(n,i+h)
     b       = max(1,i-h)
     YBar[i] = mean(Y[b:a])
  }
  return(sd(Y-YBar))
}
