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
  c(n, generalInfo(y, model = "GGP", with.loglik = T),
    generalInfo(y, model = "EGP", with.loglik = T),
    generalInfo(y, model = "PH", with.loglik = T))
  
})

info.tab %<>% t()

colnames(info.tab) <- c("recordLen",
  apply(expand.grid(c("loglik", "AIC", "BIC"), c("GGP", "EGP","PH")), 1, paste, collapse="."))

info.tab <- data.frame(ID = paste("ID", 1:49, sep = ""), info.tab)
setcolorder(info.tab, c("ID", "recordLen", 
  apply(arrange(expand.grid(c("loglik", "AIC", "BIC"), c("GGP", "EGP","PH")), Var1),
    1, paste, collapse=".")))

write.table(info.tab, file = "results/infotable.csv", row.names = FALSE)
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


quantmatch <- function(q = 0.9, ids = c(2, 5, 9, 11, 15, 24, 31, 44)) {

  mat <- t(sapply(ids, function(x) {
    
    qmat <- fread(paste("results/qmat.", x, ".csv", sep = ""))
    
    qmat %>% filter(., y >= quantile(y, q)) %>% select(., contains(".sq")) %>%
      colSums() / nrow(qmat) 
    
  })) %>% data.frame()
  
  
}

###################################################################################################


### phase type including 0s 



