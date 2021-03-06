## figure.R
#.
#.  lines for saving various figures
#.
#.

source("source/functions.R", encoding = "UTF-8")
load("data/precedent.RData") 

if(!dir.exists("plot/")) dir.create("plot/")

##- loop to save fitting results of each selected station
###################################################################################################
#. 
#.
#.
## ----------------- ## 
#. loop and ph-fitting set up
selected.id = c(2, 7, 11, 24, 44)
n.iter = 8000
dist.type = 4 # coxian 
is.logph = FALSE # default
## ----------------- ## 
#.
#.

for(id in selected.id) {
  
  y <- by_station(ID = id)
  
  models <- c("exp", "gamma", "kappa", "GGP", "EGP", "PH")
  par(mfrow = c(2,3))
  
  if (file.exists(paste("results/qmat.", id, ".csv", sep = ""))) {
    
    qmat <- fread(paste("results/qmat.", id, ".csv", sep = ""))
    for(m in models[1:5]) generalQQplot(y, model = m)
    
    qqplot(qmat$PH, y, main = "Model : PH", xlim = c(0, max(y)),
      xlab = "Theoretical Quantiles", ylab = "Observed Sample")  
    abline(0, 1, col = 'red', lty = 2)
  
  } else {
    
    for (m in models) generalQQplot(y, model = m)
    
  }
  
  dev.copy2pdf(file = paste("plot/QQplots_", "ID", id, ".pdf", sep = ""), width = 8)
  dev.off()
  
}
###################################################################################################


##- Comparing multi componenet models 
###################################################################################################
#. 
#. this part must be runned after analysis.R 
#.
## ----------------- ## 
#. loop and ph-fitting set up
selected.id = c(2, 11, 44)

## ----------------- ## 
#.
#.

layout(matrix(1:9, ncol = 3, byrow = TRUE))
for(id in selected.id) {
  
  qmatName <- paste("results/qmat.", id, ".csv", sep = "")
  qmat <- fread(qmatName)
  y <- qmat$y
  
  qqplot(qmat$GGP, y, main = "Model : GGP", xlim = c(0, max(y)),
    xlab = "Theoretical Quantiles", ylab = "Ordered Sample")
  text(x = max(y) * 0.9, y = max(y) * 0.2, paste("ID", id, sep = ""), cex = 1.5)
  abline(0, 1, col = 'red', lty = 2)
  
  qqplot(qmat$EGP, y, main = "Model : EGP", xlim = c(0, max(y)),
    xlab = "Theoretical Quantiles", ylab = "Ordered Sample")
  text(x = max(y) * 0.9, y = max(y) * 0.2, paste("ID", id, sep = ""), cex = 1.5)
  abline(0, 1, col = 'red', lty = 2)
  
  qqplot(qmat$PH, y, main = "Model : PH", xlim = c(0, max(y)),
    xlab = "Theoretical Quantiles", ylab = "Ordered Sample")
  text(x = max(y) * 0.9, y = max(y) * 0.2, paste("ID", id, sep = ""), cex = 1.5)
  abline(0, 1, col = 'red', lty = 2)

}

dev.copy2pdf(file = "plot/QQplot_multicomp.pdf", width = 8, height = 8) 
dev.off()
###################################################################################################


##- empirical distribution of estimates of gpd shape parameter 
###################################################################################################
hist(Mmat$xi, breaks = 70, probability = TRUE, xlim = c(-0.8,1),
  xlab = expression( hat(xi) ), main = expression(paste("Empirical distribution of ", hat(xi))), 
  col = 'grey')

# drawing normal density curve 
x <- seq(-2, 2,length=1000)
y <- dnorm(x, mean = mean(Mmat$xi), sd = sd(Mmat$xi))
lines(x, y, type="l", lwd = 2) ; rm(list = c('x', 'y'))

dev.copy2pdf(file = "plot/dist_xi.pdf") ;
dev.off()
###################################################################################################


##- lmrd with USHCN and GHCN
###################################################################################################
## ids of ushcn texas data set 
IDs <- ushcn.tx %>% select(COOP_ID) %>% unique() %>% unlist(., use.names = FALSE)
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

plot(x = Mmat$t_2, y = Mmat$t_3, pch = 20, cex = 0.7,
  ylab = expression(italic('L') * "-skewness"), xlab = expression(italic('L') * "-variation"), 
  ylim = c(0, 1), xlim = c(0.2, 1))

# add regression line
reg <- lm(Mmat$t_3 ~ Mmat$t_2)
abline(reg, lty = 2)
legend(x = 0.2, y = 0.9, legend = "regression line", lty = 2, box.col = 'transparent')

# point texas data 
points(x = Lmat.tx$t_2, y = Lmat.tx$t_3, pch = 15, col = 'red', cex = 0.7)
legend("bottomright", legend = "USHCN-Texas", col = 'red', pch = 15, cex = 0.9,
  box.col = 'transparent' )

dev.copy2pdf(file = "plot/lmrd.pdf", width = 8) ; dev.off()
###################################################################################################


