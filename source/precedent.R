## precedent.R
#.
#. Reproducing some results from precedent studies 
#.  Papalexiou et al. 2012
#.

suppressPackageStartupMessages({
  
  library(data.table)
  library(dplyr)
  library(evir)
  library(lmom)

})

###################################################################################################
### 0. Preprocessing 
###################################################################################################

# station inventory
stationInv <- fread("data/stationInv.txt")

# file directory
FILE.DIR = "data/DailyFiles/"
RAW.DATA.DIR = "data/DailyRaw/"

## data selection criteria
#. missing values less than or equal to 20%
#. suspicious flags less than or equal to 0.1%
#. excluding inappropriately coded stations 
selectedInv <- stationInv %>% filter(nYear >= 50 & flag <= 0.01 & miss <= 0.2 & codingErr == FALSE)

stationlist = selectedInv$Id
filelist = paste(stationlist, ".dly", sep = "")

###################################################################################################
### 1. L-variation vs L-skewness from Paplexitou et al. [2012]
###################################################################################################

# matrix of l-moment ratios (~ 4-th moment)
Lmat <- sapply(filelist, function(x) {
  
  vec <- fread(paste(FILE.DIR, x, sep = "")) %>% unlist(., use.names = FALSE)
  samlmu(vec, nmom = 4, ratios = TRUE)
  
  }, USE.NAMES = FALSE)

Lmat <- t(Lmat) %>% data.table()
Lmat <- cbind(Id = stationlist, Lmat)
Lmat[, t_2 := l_2 / l_1]
setcolorder(Lmat, c("Id", "l_1", "l_2", "t_2", "t_3", "t_4"))

# LMRD, L-variation vs L-skewness
plot(x = Lmat$t_2, y = Lmat$t_3, pch = 16,
  ylab = expression(italic('L') * "-skewness"), xlab = expression(italic('L') * "-variation"), 
  ylim = c(0, 1), xlim = c(0, 1))

# average point 
points(x = mean(Lmat$t_2), y = mean(Lmat$t_3), col = 'red', pch = 15)

if(!dir.exists("plot/")) dir.create("plot/")
dev.copy2pdf(file = "plot/lmrd.pdf", width = 8) ; dev.off()

##### add parametric line on LMRD - GPD, GG ,BurrXII




###################################################################################################
### 2. Peaks Over Threshold - estimating GPD shape parameter, Papalexiou et al.[2013 HESS]
###################################################################################################

## defining tail, Papalexiou et al.[2013 HESS]
#. number of years := number of extreme values 
#. then estimate shape parameter 

## matrix of gpd parameter estimates 
#. fit gpd with same criteria in Papalexiou et al. [2013 HESS]
#. if error, trying 0.99 threshold 

GPmat <- sapply(filelist, function(x) {
  
  vec <- fread(paste(FILE.DIR, x, sep = "")) %>% unlist(., use.names = FALSE)
  n_extrm <- round(selectedInv[Id == substr(x, 1, nchar(x)-4)]$nYear)
  fit <- tryCatch(gpd(data = vec, nextremes = n_extrm, method = "ml"),
    error = function(e) gpd(data = vec, threshold = 0.99, method = "ml"))
  fit$par.ests
  
}, USE.NAMES = FALSE)

GPmat <- t(GPmat) %>% data.table()
GPmat <- cbind(Id = stationlist, GPmat)

# empirical distribution of estimates of gpd shape parameter 
hist(GPmat$xi, breaks = 70, probability = TRUE, xlim = c(-0.8,1),
  xlab = expression( hat(xi) ), main = expression(paste("Empirical distribution of ", hat(xi))), 
  col = 'grey')

# drawing normal density curve 
x <- seq(-2, 2,length=1000)
y <- dnorm(x, mean = mean(GPmat$xi), sd = sd(GPmat$xi))
lines(x, y, type="l", lwd = 2) ; rm(list = c('x', 'y'))

dev.copy2pdf(file = "plot/dist_xi.pdf", width = 8) ; dev.off()

Mmat <- merge(Lmat, GPmat, by = "Id") 
rm(list = c('Lmat', 'GPmat'))

# save image to shorten computing time 
save.image(file = "data/precedent.RData")
