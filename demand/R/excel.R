convert.xls2csv <- function(xls.file, csv.file = "./temp.csv") {
  # converts XLS file to CSV
  # needs command line tool Gnumeric installed (only in mac or linux)
  # http://gnumeric.org/
  # to install in mac: brew install gnumeric
  
  xls.file <- paste0("\'",xls.file,"\'")
  
  cat("\nConverting xls file to csv... Please wait...\n")
  
  out <- system(paste0("ssconvert -S ",xls.file, " ",csv.file))
  
  if (out == 0) {
    cat("\nXls file converted to csv successfully!\n")
  } else {
    cat("\nError converting xls file to csv!\n")
    cat(paste0("Please check name of file and if command line tool ",
               "'Gnumeric' is installed"))
  }
  
  return(csv.file <- csv.file) # just to make sure it won't print in screen
}