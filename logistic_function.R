#Code found here: https://stats.stackexchange.com/questions/47802/whats-the-most-pain-free-way-to-fit-logistic-growth-curves-in-r
###Log fit - use quotes around the variable names 
log.fit <- function(dep, ind, data){
  #Self-starting...
  
  y <- data[, dep]
  x <- data[, ind]
  
  
  
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