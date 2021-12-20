###Log fit - be sure to use quotes around the variable names in the call
log.fit2 <- function(dep, ind, yourdata){
  #Self-starting...
  
  y <- yourdata[, dep]
  x <- yourdata[, ind]
  library(minpack.lm)
  log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
  #log.ss <- nlsLM(y ~ SSlogis(x, phi1, phi2, phi3) ,alg="plinear")
  
  
  phi1 <- summary(log.ss)$coef[1]
  
  phi2 <- summary(log.ss)$coef[2] #exp((summary(log.ss)$coef[2]) * (1/summary(log.ss)$coef[3]))
  
  phi3 <- summary(log.ss)$coef[3]
  
  
  out <- data.frame(cbind(c(phi1=phi1, phi2=phi2, phi3=phi3)))
  names(out)[1] <- "Logistic Curve"
  
  return(out)
}