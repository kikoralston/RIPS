year <- 1993

downloadStationbyUSAF(725200, year, save.to.temp = FALSE)
a <- gzfile(paste0("./725200-94823-", year, ".gz"))
b <- readLines(a)
close(a)
c <- b[grep("AL1", b)]
d <- sapply(c, FUN=function(x) { pos <- regexpr("AL1", x); as.numeric(substr(x,pos+3, pos+4)) })
names(d) <- NULL
e <- sapply(c, FUN=function(x) { substr(x,16, 27) })
names(d) <- e
