recessions <- function(){
  # list of recessions
  folder.csv <- "~/GoogleDrive/CMU/RIPS/socioeconomic/"
  df.recessions <- read.csv(file=paste0(folder.csv, "NBER.csv"), 
                            stringsAsFactors = FALSE)
  
  # clean
  df.recessions <- df.recessions[ ,c(1:4)]
  df.recessions <- df.recessions[complete.cases(df.recessions), ]
  df.recessions$Peak.month <- paste("01", df.recessions$Peak.month)
  df.recessions$Trough.month <- paste("01", df.recessions$Trough.month)
  df.recessions$Peak.month <- as.Date(df.recessions$Peak.month, 
                                      format = "%d %B %Y")
  df.recessions$Trough.month <- as.Date(df.recessions$Trough.month, 
                                        format = "%d %B %Y")
  
  return(df.recessions)
}