## precedent.R
#.
#. Reproducing some results from precedent studies 
#.  Papalexiou et al. 2012, Li et al. 2012
#.

suppressPackageStartupMessages({
  
  library(data.table)
  library(dplyr)
  library(evir)
  library(lmom)
  library(ffbase)
  library(LaF)
  library(magrittr)
  library(R.utils)
  
})

##- GHCN

###################################################################################################
### 1-1. GHCN : Preprocessing 
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
### 1-2. GHCN : L-variation vs L-skewness from Paplexitou et al. [2012]
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


###################################################################################################
### 1-3. GHCN : Peaks Over Threshold - estimating GPD shape parameter, Papalexiou et al.[2013 HESS]
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

Mmat <- merge(Lmat, GPmat, by = "Id") 
rm(list = c('Lmat', 'GPmat'))


##- USHCN Texas


#. 49 weather stations in Texas 
#. 10 climate divisions divided by National Weather Service 
#. probably : Amarillo, Austin/San Antonio, Brownsville, Corpus Christi, 
#.   Dallas/Fort Worth, El Paso, Houston/Galveston, Lubbock, Midland/Odessa, San Angelo

###################################################################################################
### 1. USHCN : data import and tidy up 
###################################################################################################

url <- "http://cdiac.ornl.gov/ftp/ushcn_daily/state41_TX.txt.gz"
ushcn.zip <- "data/data/state41_TX.txt.gz"
ushcn.raw <- "data/state41_TX.txt"

if(!file.exists(ushcn.raw) & !file.exists(ushcn.zip))
  download.file(url, ushcn.zip)

if(!file.exists(ushcn.raw) & file.exists(ushcn.zip))
  R.utils::gunzip(ushcn.zip, ushcn.raw)

ushcn.tx.laf <- laf_open_fwf(ushcn.raw,
  column_widths = c(6, 4, 2, 4,  rep(c(5, 1, 1, 1), 31) ),
  column_types = c('character', 'integer', 'integer', 'character',
    rep(c('integer', 'character', 'character', 'character'), 31)))


ushcn.tx <- laf_to_ffdf(ushcn.tx.laf)

ushcn.tx <- as.data.table(ushcn.tx)
# missing values  
ushcn.tx[ushcn.tx == -9999] = NA

# assign column names 
colnames(ushcn.tx) <- c("COOP_ID", "YEAR", "MONTH", "ELEMNT",
  paste(c("VALUE", "MFLAG", "QFLAG", "SFLAG"), rep(1:31, each = 4), sep = "" ))

ushcn.tx %<>% filter(ELEMNT == "PRCP") %>% select(-ELEMNT)

tidy.ushcn.tx <- ushcn.tx %>% filter(., YEAR >= 1940, YEAR <= 2009) %>% 
  select_(.dots = c("COOP_ID", "YEAR", "MONTH", paste("VALUE", 1:31, sep = "")))

# 
rm(list = c('url', 'ushcn.raw', 'ushcn.zip', 'ushcn.tx.laf'))


# save image to shorten computing time 
save.image(file = "data/precedent.RData")
