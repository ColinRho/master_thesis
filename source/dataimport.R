## dataimport.R
#.
#. This process may take more than a whole day 
#. over 8 GB 
#.
#. Criteria
#.  1. record length greater than or equal to 50 yrs
#.  2. % of missing and of suspicious flags are recorded
#.
#. <A measurement of precipitation is a tenths of mm >
#. 
#. This script should be runned only once before the others 
#. 
#.  8.3GB raw sets
#.  653M + 1.4M minimal sets
#. 

suppressPackageStartupMessages({
  
  library(data.table)
  library(dplyr)
  library(GhcnDaily)
  library(magrittr)
  library(rnoaa)
  
})

###################################################################################################

#- IMPORTING

## download daily inventory 
downloadDailyInventory(directory = "data/")

## call daily inventory which contains precipitation data
PRCP <- readDailyInventory(elements = "PRCP", filename = "data/DailyInv.txt")

## record length greater than or equal to 50 years 
## this set contatins information of each station
stationInv <- PRCP %>% mutate(YearLength = LastYear - FirstYear) %>% 
  filter(YearLength >= 50) %>% select(-YearLength, -Element) %>% data.table()

m = nrow(stationInv)

## directories to store data sets(raw sets, converted sets)
RAW.DATA.DIR = "data/DailyRaw/"
FILE.DIR = "data/DailyFiles/"

## only precipitation data
extractPRCP <- function(id) {
  
  path = paste(RAW.DATA.DIR, stationInv$Id[id], ".dly", sep = "")
  
  if (file.exists(path)) {
    
    dat = fread(path)
    dat %<>% filter(., element == 'PRCP') %>% select(., -element)
    write.table(dat, path, row.names = FALSE)
    
  } else {
    
    print("no data")
  
  }
  
}

## raw data download loop
err.list <- c()
for (id in 1:m) {
  
  fpath <- paste(RAW.DATA.DIR, stationInv$Id[id], ".dly", sep = "")
  
  # file is not downloaded yet
  if (!file.exists(fpath)) {
    
    tryCatch(invisible(ghcnd(stationid = stationInv$Id[id], path = RAW.DATA.DIR)),
      error = function(err) { err.list = append(err.list, stationInv$Id[id]) ; cat("error", "\n")})
    
    # only PRCP data
    extractPRCP(id = id)
    cat(".")
    
  } else { # file exists
    
    cat(".")
    
  }
  
  # print processing 
  if (id %% 50 == 0) 
    cat("\n",id/m * 100, "% done", "\n")
  
  if (id == m)
    cat("\n", "100% done", "\n")
  
}

###################################################################################################

#- REARRANGEMENT

## a directory to record minimal data sets
if(!dir.exists(FILE.DIR)) dir.create(FILE.DIR)

filelist <- list.files(path = RAW.DATA.DIR)

## extract only strictly positive values, add some summary statistics on 'stationInv' 
record_arrange <- function(filename, PATH = RAW.DATA.DIR, DIR = FILE.DIR) {
  
  if(!file.exists(paste(DIR, filename, sep = ""))) {
  
    # import raw data
    dt = fread(input = paste(PATH, filename, sep = ""))
    stationid = substr(filename, 1, nchar(filename) - 4)
    
    VALUE <- paste("VALUE", 1:31, sep = "")
    QFLAG <- paste("QFLAG", 1:31, sep = "")
    
    # arrange vectors
    values = dt %>% select_(., .dots = VALUE) %>% unlist(., use.names = F)
    values[values == -9999] = NA ; n =length(values)
    flags = dt %>% select_(., .dots = QFLAG) %>% unlist(., use.names = F)
    
    # some summary stats 
    dryDays = length(values[values == 0 & !is.na(values)])
    missingRatio = round(sum(is.na(values)) / n , 3)
    badFlagRatio = round(sum(flags == "G" | flags == "X", na.rm = T) / n, 3)
    
    # strictly positive values 
    pos = values[values != 0 & !is.na(values)]
    
    # add summaries on 'stationInv'
    stationInv[Id == stationid, nYear := nrow(dt)/12]  
    stationInv[Id == stationid, recordLen := n]
    stationInv[Id == stationid, dry := dryDays]
    stationInv[Id == stationid, miss := missingRatio]
    stationInv[Id == stationid, flag := badFlagRatio]
    
    # coding error & no records
    # stations contaion negative values 
    if(length(pos) != 0)
      stationInv[Id == stationid, codingErr := (sum(pos < 0) != 0)]
    
    else  # no records 
      stationInv[Id == stationid, codingErr := NA]
    
    # write files
    write.table(pos, paste(DIR, filename, sep = ""),
      row.names = FALSE)  
    
  } 
  
}

invisible(sapply(filelist, record_arrange))

## update stationInv  
write.table(stationInv, "data/stationInv.txt", row.names = FALSE)
    
