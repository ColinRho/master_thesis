## USHCN.R

suppressPackageStartupMessages({
  
  library(data.table)
  library(dplyr)
  library(ffbase)
  library(LaF)
  library(magrittr)
  library(R.utils)
  
})

####
###


# 49 weather stations in Texas 
# 10 climate divisions divided by National Weather Service 
# probably : Amarillo, Austin/San Antonio, Brownsville, Corpus Christi, 
#   Dallas/Fort Worth, El Paso, Houston/Galveston, Lubbock, Midland/Odessa, San Angelo

###################################################################################################
### 1. data import and tidy up 
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


prep.tx <- laf_to_ffdf(ushcn.tx.laf)

prep.tx <- as.data.table(prep.tx)
# missing values  
prep.tx[prep.tx == -9999] = NA

# assign column names 
colnames(prep.tx) <- c("COOP_ID", "YEAR", "MONTH", "ELEMNT",
  paste(c("VALUE", "MFLAG", "QFLAG", "SFLAG"), rep(1:31, each = 4), sep = "" ))

prep.tx %<>% filter(ELEMNT == "PRCP") %>% select(-ELEMNT)

tidy.prep.tx <- prep.tx %>% filter(., YEAR >= 1940, YEAR <= 2009) %>% 
  select_(.dots = c("COOP_ID", "YEAR", "MONTH", paste("VALUE", 1:31, sep = "")))

# 
rm(list = c('url', 'ushcn.raw', 'ushcn.zip', 'ushcn.tx.laf'))

###################################################################################################
### 2. some functions
###################################################################################################


# flattening function
flat <- 
  
  # input vector y should be integer 
  
  function(y) {
    
    # flattening data
    z <- table(y)
    zz <- c()
    for ( i in 1:length(z) ) {
      
      interval <- c(as.numeric(names(z[i])) - 0.5, as.numeric(names(z[i])) + 0.5)
      v <- seq(from = interval[1], to = interval[2], length.out = z[i] + 2)
      zz <- c(zz, v[-c(1, length(v))])
      
    }
    
    zz
    
  }


# a function extracting daily data by each station without missing value and 0's 
by_station <-
  
  # insert a positive integer from 1 to 49 corresponding to station ID number needed 
  # flattening 
  
  function(ID, in.mm = TRUE, flat = FALSE) {
    
    id.value <- 
      tidy.prep.tx %>% select(COOP_ID) %>% unique() %>% slice(ID) %>% as.integer
    
    # extract a vector
    y <- 
      tidy.prep.tx %>% filter(as.integer(COOP_ID) == id.value) %>% 
      select(-c(COOP_ID, YEAR, MONTH)) %>% c() %>% unlist(., use.names = F) %>%
      subset(!(is.na(.)) & . != 0)
    
    if (flat) v <- flat(y) 
    else v <- y 
    
    # change unit(1/100 inch to mm)
    if (in.mm) v * 0.254
    else v
    
  }

save.image(file = "data/ushcn.RData")
