# ***************************************************************************
# Functions to download and read data from the FERC repository
# ***************************************************************************

get.url.714 <- function() {
  # gets string with url of form 714 in FERC website
  url.714 <- 'http://www.ferc.gov/docs-filing/forms/form-714/data/'
  return(url.714)
}

get.name.zip.714 <- function() {
  # gets string with name of zip file with data of form 714
  name.zip.714 <- 'form714-database.zip'
  return(name.zip.714)
}

downloadFERCdata <- function(ferc.zip){
  
  url.ferc <- get.url.714()
  zip.name <- get.name.zip.714()
  
  #ferc.zip <- tempfile()
  download.file(paste0(url.ferc, zip.name), ferc.zip)

  return(ferc.zip)
  
}

readListFERC <- function(get.online = TRUE,
                         ferc.zip = 
                           "/Users/Kiko/GoogleDrive/CMU/RIPS/FERC/after2006/form714-database.zip"){
  
  if (get.online) ferc.zip <- downloadFERCdata(ferc.zip=ferc.zip)
  
  # reads lists with names of utilities and respective codes
  respondent.id.file <- "form714-database/Respondent IDs.csv"
  list.names <- read.csv(file = unz(ferc.zip, filename = respondent.id.file),
                         stringsAsFactors = FALSE)
  
  list.names[ ,2] <- tolower(list.names[ ,2])
  
  return(list.names)
}


readFercNew <- function(get.online = TRUE,
                        ferc.zip = 
                          "/Users/Kiko/GoogleDrive/CMU/RIPS/FERC/after2006/form714-database.zip",
                        name.util="tennessee") {
  # Reads data from the FERC repository's zip file, either online or previoulsy 
  # downaloaded. This function only reads data more recent than 2006
  # source of data: http://www.ferc.gov/docs-filing/forms/form-714/data.asp
  #
  # Args:
  #   get.online: (boolean) TRUE if data should be read directly from online
  #               repository. FALSE if data was downloaded previously
  #   path.data: (string) path where the FERC data was saved
  #   name.util: (string) name (of part of name) of utility being searched
  #
  # Returns:
  #   data frame with load data after 2006
  
  # reads lists with names of utilities and respective codes
  list.names <- readListFERC(get.online, ferc.zip)
  
  # gets position and code for utility
  idx.tva <- grep(tolower(name.util), list.names$respondent_name)
  code.tva <- list.names$respondent_id[idx.tva]
  
  # reads whole data base
  name.file <- paste0("form714-database/",
                      "Part 3 Schedule 2 - Planning Area Hourly Demand.csv")
  data.ferc <- read.csv(file = unz(ferc.zip, filename =  name.file), 
                        stringsAsFactors = FALSE)
  
  # filters data from TVA
  data.tva <- data.ferc[data.ferc$respondent_id==code.tva, ]
  
  # remove data.ferc and call garbage collector
  rm(data.ferc)
  gc()
  
  # get demand data from data frame and convert to 1-D vector
  
  dates.load <- data.tva[ ,c("plan_date")]
  ls.hours <- lapply(X = dates.load, 
                     FUN = {function(x) {
                       return(seq(from=as.POSIXlt(x, 
                                                  format="%m/%d/%Y %H:%M:%S",
                                                  tz="UTC"), 
                                  by="hour", length.out = 24))}})
  
  hours.vector <- ldply(.data=ls.hours, .fun={ function(x){data.frame(time=x)}})
  
  a <- as.matrix(data.tva[,8:31])
  a <- t(a)
  nr <- dim(a)[1]
  nc <- dim(a)[2]
  dim(a) <- nr*nc
  a <- as.vector(unlist(a))
  
  df.ferc.1 <- data.frame(time = hours.vector, load = a)
  
  #if (get.online) {
  #  unlink(ferc.zip)
  #}
  
  return(df.ferc.1)
}

# ******** unzip files from FERC (pre 2006) ******
readFercOld <- function(read.zip = FALSE, flag.download = FALSE) {
  # Reads data from the FERC repository saved at 
  # path.data. This function only reads data older than 2006
  # source of data: http://www.ferc.gov/docs-filing/forms/form-714/data.asp
  #
  # Args:
  #   read.zip: (boolean) TRUE if zip files should be read
  #   flag.download: (boolean) TRUE if zip files should be downloaded. 
  #                   If FALSE, files will be read from local folder
  #
  # Returns:
  #   data frame with load data older than 2006
  
  #library(xlsx)
  
  if (read.zip) {
    url.ferc <- get.url.714()
    
    # paths of SERC zip files
    df.paths.zip <- data.frame(path=c("2004/SERC.zip", "2003/SERC.zip", 
                                      "2002/SERC3.zip", "2001/01SERC.zip", 
                                      "2000/2000SERC.zip", "99SERC.zip",
                                      "98SERC.zip", "97serc1.zip", "97serc2.zip", 
                                      "96SERC1.zip", "96SERC2.zip", "95SERC1.zip",
                                      "95SERC2.zip", "94SERC1.zip", "94SERC2.zip",
                                      "93SERC1.zip", "93SERC2.zip"),
                               year=c(2004, 2003, 2002, 2001, 2000, 1999, 1998,
                                      1997, 1997, 1996, 1996, 1995, 1995, 1994, 
                                      1994, 1993, 1993))
    
    files.tva.total <- c()
    
    if (flag.download) {
      n.files <- nrow(df.paths.zip)
    } else {
      a <- list.files()
      n.files <- length(a)
    }
    
    for (i in 1:n.files) {
      if (flag.download) {
        temp <- tempfile()
        download.file(paste0(url.ferc, df.paths.zip$path[i]),
                      temp)
      } else {
        temp <- a[i]
      }
      
      list.of.files.zip <- unzip(zipfile = temp, list = TRUE)
      list.of.files.zip <- list.of.files.zip[list.of.files.zip$Length > 0, ] 
      irow <- grep("Tennessee|tva", list.of.files.zip$Name, 
                   ignore.case=TRUE)
      if(length(irow) > 0) {
        files.tva <- list.of.files.zip$Name[irow]
        if (flag.download) {
          path.save <- paste0("./FERC/pre2006/", df.paths.zip$year[i]) 
        } else {
          path.save <- paste0("./FERC/pre2006/", substring(a[i], 1, nchar(a[i])-4))
        }
        
        dir.create(path=path.save, recursive = TRUE)
        unzip(zipfile = temp, files=files.tva, junkpaths = TRUE, 
              exdir = path.save)
        
        files.tva.total <- c(files.tva.total, files.tva)
      }
      if (flag.download) { 
        unlink(temp) 
      }
      #i <- i + 1
    }
    files.tva.total
  }
  
  path.data <- paste0(getwd(), "/FERC/pre2006")
  
  dir.years <- list.dirs(path = path.data, full.names = FALSE, 
                         recursive = FALSE)
  exclude <- match(c(".", "zip"), dir.years)
  if (sum(is.na(exclude)) > 0) {
    exclude <- exclude[!is.na(exclude)]
  }
  if (length(exclude) > 0) {
    dir.years <- dir.years[-exclude]
  }
  
  name.files <- c("TVALOAD.PRN", "LOADS.PRN", "TVALOAD.PRN", "TVALOAD.PRN", 
                  "output.prn", "TVALOADS.PRN", "CY99.csv", "CY2000.csv", 
                  "TVA2001Load&Cost.csv", "cy02pscosts.csv", 
                  "CY03714hourlydata.xls")
  
  path.files <- paste("./FERC/pre2006", dir.years, name.files, sep="/")
  
  df.ferc.2 <- data.frame(date=as.Date(character()), 
                          hour=numeric(), 
                          load=numeric())
  for (i in 1:length(path.files)) {
    ext <- substring(path.files[i], first = nchar(path.files[i])-2)
    if (tolower(ext) == "prn") {
      dat <- read.table(file=path.files[i], stringsAsFactors = FALSE)
      dat <- dat[-which(!complete.cases(dat[,c(1:3)])), ]
      if (nchar(dat[1,1]) == 8) {
        date.format <- "%Y%m%d"
      } else {
        date.format <- "%y%m%d"
      }
      dat[,1] <- as.Date(dat[,1], format=date.format, tz='UTC')
      names(dat) <- names(df.ferc.2)
      df.ferc.2 <- rbind(df.ferc.2, dat)
    } else {
      # check cases where the file is xls and convert to csv
      if (tolower(ext) == 'xls' || tolower(ext) == 'xlsx') {
        csv.file <- convert.xls2csv(xls.file = path.files[i], 
                                    csv.file = "./temp.csv")
        out <- system(paste0("cp ", csv.file, " ", dirname(path.files[i])))
        csv.file <- paste0(csv.file, '.0')
        path.files[i] <- paste0(dirname(path.files[i]), '/', basename(csv.file))
      }
      
      dat <- read.csv(file=path.files[i], stringsAsFactors = FALSE, 
                      strip.white = TRUE)
      if (sum(!complete.cases(dat[,c(1:3)])) > 0 ) {
        dat <- dat[-which(!complete.cases(dat[,c(1:3)])), ]
      }
      dat <- dat[ ,c(1:3)]
      if (class(dat[ ,3]) == "character") {
        dat[ ,3] <- gsub(",", "", dat[ ,3])
        dat[ ,3] <- as.numeric(dat[ ,3])
      }
      dat[,1] <- as.Date(as.character(dat[ ,1]), format="%m/%d/%y", tz='UTC')
      names(dat) <- names(df.ferc.2)
      df.ferc.2  <- rbind(df.ferc.2, dat)
    }
  }
  df.ferc.2$hour <- (df.ferc.2$hour - 1)
  
  time <- as.POSIXct(paste0(as.character(df.ferc.2$date), " ", 
                            sprintf("%d:00:00", df.ferc.2$hour)),
                     "%Y-%m-%d %H:%M:%S", tz = "UTC")
  df.ferc.2 <- data.frame(time=time, df.ferc.2)
  df.ferc.2 <- subset(df.ferc.2, select = -c(date, hour))
  
  #plot(x=df.ferc$time, y=df.ferc$load, type="l")
  #lines(x=df.ferc$time, 
  #      y=filter(x=df.ferc$load, filter=rep(1/(365*24), 365*24)), 
  #      col="red")
  
  return(df.ferc.2)  
}