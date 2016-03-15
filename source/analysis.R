## Analysis.R
#.
#.
#.
#.
#.

source("source/functions.R", encoding = "UTF-8")
load("data/precedent.RData") 


##- AIC BIC table of USCHN Texas
###################################################################################################

info.tab <- sapply(1:49, function(x) {
  
  y <- by_station(x)
  n <- length(y)
  c(n, generalInfo(y, model = "GGP"),
    generalInfo(y, model = "EGP"),
    generalInfo(y, model = "PH"))
  
})

info.tab %<>% t()

colnames(info.tab) <- c("recordLen",
  apply(expand.grid(c("AIC", "BIC"), c("GGP", "EGP","PH")), 1, paste, collapse="."))

info.tab <- data.frame(ID = paste("ID", 1:49, sep = ""), info.tab)
setcolorder(info.tab, 
  c("ID", "recordLen", "AIC.GGP", "AIC.EGP", "AIC.PH", "BIC.GGP", "BIC.EGP", "BIC.PH"))

write.table(info.tab, file = "results/infotable.txt", row.names = FALSE)
###################################################################################################



##- MSEs 
###################################################################################################

save.qmat <- function(id, model = "GGP EGP PH") {
  
  y = by_station(id)
  qmat <- generalQQplot(y, model = model, as.quantile = TRUE)
  qmat <- data.table(y = sort(y), qmat)
  qmat[, c("GGP.sq", "EGP.sq", "PH.sq") := list((y-GGP)^2, (y-EGP)^2, (y-PH)^2)]
  write.csv(qmat, file = paste("results/qmat.", id, ".csv", sep = ""), row.names = FALSE)

}


ids <- c(2, 7, 9, 27, 31, 48)
sapply(ids, save.qmat)
###################################################################################################


##- identify the shape of sample distribution through l-moments ratio diagram
###################################################################################################
#.  Papalexiou et al. [2012]
#.  round 10 stations each, J-shape and Bell-shape 

Mmat %>% filter(., t_2 >= 0.65 & t_3 >= 0.65) -> JShape
Mmat %>% filter(., t_2 <= 0.43 & t_3 <= 0.43) -> BellShape


#. take 1 bell-shaped sample and do same analysi.
y <- get.station.data(BellShape$Id[1])
hist(y, breaks = 60)
generalQQplot(y, model = "PH")

qmat.bell <- generalQQplot(y, model = "GGP EGP PH", as.quantile = TRUE)
qmat.bell <- data.table(y = sort(y), qmat.bell)
write.csv(qmat.bell, file = "results/qmat.bell.csv", row.names = FALSE)

generalQQviaQmat(qmat.bell, layout.matrix = matrix(1:3, ncol = 3))
dev.copy2pdf(file = "plot/QQplot_bell.pdf", width = 9, height = 4) ; dev.off()
generalInfo(y, model = "GGP EGP PH")
###################################################################################################



qmat <- fread("results/qmat.24.csv")
qmat %>% filter(., y >= quantile(y, 0)) %>% select(., contains(".sq")) %>% colSums()



y = get.station.data(BellShape$Id[3])


