## Analysis.R
#.
#.
#.
#.
#.

source("source/functions.R", encoding = "UTF-8")
load("data/ushcn.RData")
load("data/precedent.RData") 

##- identify the shape of sample distribution through l-moments ratio diagram
###################################################################################################
#.  Papalexiou et al. [2012]
#.  round 10 stations each, J-shape and Bell-shape 

Mmat %>% filter(., t_2 >= 0.7 & t_3 >= 0.7) -> JShape
Mmat %>% filter(., t_2 <= 0.43 & t_3 <= 0.43) -> BellShape

## sample histograms of each shapes of distribution
par(mfrow = c(2,2))
# a case of J-shaped sample
get.station.data(JShape$Id[8]) %>% 
  hist(., breaks = 60, main = "J-shaped sample distribution")

# log 
get.station.data(JShape$Id[8]) %>% 
  log %>% hist(., breaks = 60, main = "J-shaped sample distribution(log)")

# a case of bell-shaped sample
get.station.data(BellShape$Id[2]) %>% 
  hist(., breaks = 60, main = "Bell-shaped sample distribution")

get.station.data(BellShape$Id[1]) %>% 
  log %>% hist(., breaks = 60, main = "Bell-shaped sample distribution(log)")

par(mfrow = c(1,1))


## ids of ushcn texas data set 
IDs <- prep.tx %>% select(COOP_ID) %>% unique() %>% unlist(., use.names = FALSE)
## 43 stations
subset.ghcn.tx <- selectedInv %>%
  filter(substr(Id, 1, 2) == 'US') %>% filter(substr(Id, 6, 12) %in% IDs) 

## excluding these stations of tx 
which(!(IDs %in% substr(subset.ghcn.tx$Id, 6, 12)))

### USHCN data seems J-shape data according to L-moments ratio plot
### EMPHT fit well both Bell- and J-shape data


# calculate l-moments of prep.tx
Lmat.tx <- sapply(1:49, function(x) {
  
  vec = by_station(x)
  lmom::samlmu(vec, nmom = 4, ratios = TRUE)
  
})

Lmat.tx <- t(Lmat.tx) %>% data.table()
Lmat.tx <- cbind(COOP_ID = IDs, Lmat.tx)
Lmat.tx[, t_2 := l_2 / l_1]
setcolorder(Lmat.tx, c("COOP_ID", "l_1", "l_2", "t_2", "t_3", "t_4"))


plot(x = Lmat$t_2, y = Lmat$t_3, pch = 16,
  ylab = expression(italic('L') * "-skewness"), xlab = expression(italic('L') * "-variation"), 
  ylim = c(0, 1), xlim = c(0, 1))

points(x = Lmat.tx$t_2, y = Lmat.tx$t_3, pch = 15, col = 'red')
###################################################################################################




##- loop to save ph-fitting results in plot 
###################################################################################################
#. 
#.
#.
## ----------------- ## 
#. loop set up
selected.id = c(7, 11, 24, 44)
n.iter = 8000
dist.type = 4 # coxian 
is.logph = FALSE # default
dir = "../plot/"
## ----------------- ## 
#.
#.
for(id in selected.id) {
  
  imagename <- paste("ID", id, ".pdf", sep = "") 
  fit <- phase.select(x = by_station(ID = id), iter = n.iter, type = dist.type, max.phase = 7)
  dev.copy2pdf(file = paste(dir, imagename, sep = ""), width = 8) ; dev.off()
  
  ph.qqplot(fit = fit)
  dev.copy2pdf(file = paste(dir, "QQplot", imagename, sep = ""), width = 8) ; dev.off()
  
}
###################################################################################################





##- reproduce some plots in precedent studies 
###################################################################################################
#.
#.
#.

y = by_station(11)

#. 1) Exponential

par(mfrow = c(2,3))

general.qqplot(y, model = "exponential")

#. 2) Gamma

general.qqplot(y, model = "gamma")

#. 3) FK08, Furrer & Katz [2008]

general.qqplot(y, model = "FK08", theta = 32)

#. 4) HEG, Li et al. [2013]

general.qqplot(y, model = "HEG")

#. 5) GG & GB2, Papaelxiou et al. [2012]

general.qqplot(y, model = "GG")
general.qqplot(y, model = "GB2")


###################################################################################################



##- Goodness of fit 
###################################################################################################

## AIC BIC



###################################################################################################




##- Evaluate estimator (efficiency, asymptotic property, etc., ...)
###################################################################################################



###################################################################################################


