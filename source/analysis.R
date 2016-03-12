## Analysis.R
#.
#.
#.
#.
#.

source("source/functions.R", encoding = "UTF-8")
load("data/precedent.RData") 

##- identify the shape of sample distribution through l-moments ratio diagram
###################################################################################################
#.  Papalexiou et al. [2012]
#.  round 10 stations each, J-shape and Bell-shape 

Mmat %>% filter(., t_2 >= 0.7 & t_3 >= 0.7) -> JShape
Mmat %>% filter(., t_2 <= 0.43 & t_3 <= 0.43) -> BellShape

##- AIC BIC table
###################################################################################################

info.tab <- sapply(1:49, function(x) {
  
  y <- by_station(x) 
  n <- length(y)
  GGP.info <- general.qqplot(y, model = "GGP", without.plot = TRUE)
  EGP.info <- general.qqplot(y, model = "EGP", without.plot = TRUE)
  ph_fit <- EMpht(y, type = 4, phases = 3, iter = 8000, curve.fit = FALSE)
  loglik <- sum(log(dphtype(y, prob = ph_fit$alpha, rates = ph_fit$T)))
  p <- length(ph_fit$alpha) ; k <- p * 2 - 1
  AIC <- 2*k - 2*loglik
  BIC <- k * log(n) - 2*loglik
  c(n, GGP.info, EGP.info, AIC, BIC)
  
})

info.tab %<>% t()

colnames(info.tab) <- c("recordLen",
  apply(expand.grid(c("AIC", "BIC"), c("GGP", "EGP","PH")), 1, paste, collapse="."))

info.tab <- data.frame(ID = paste("ID", stationIDs, sep = ""), info.tab)
setcolorder(info.tab, 
  c("ID", "recordLen", "AIC.GGP", "AIC.EGP", "AIC.PH", "BIC.GGP", "BIC.EGP", "BIC.PH"))
###################################################################################################



##- MSEs 
###################################################################################################



###################################################################################################


