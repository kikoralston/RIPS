# *****************************************************************************
# Script to read regional population scenarios from NATIONAL CLIMATE ASSESMENT 
# (NCA)
# *****************************************************************************

library(shapefiles)

setwd("~/GoogleDrive/CMU/RIPS/socioeconomic/NCA")

# URL of population projections for US counties from NCA (only Southeast)
url.ncal <- "https://edg.epa.gov/data/Public/ORD/NCEA/Southeast.zip"

# downloads and reads county populations data from NCA
tt <- tempfile()
download.file(url.ncal, tt)
unzip(tt, files = "Southeast/SoutheastPop.dbf", junkpaths = TRUE)
unlink(tt)
data.nca <- read.dbf(dbf.name = "SoutheastPop.dbf")

# gets FIPS county codes from downloaded data 
fips <- as.character(data.nca$dbf$FIPS)
fips[nchar(fips) < 5] <- paste0("0", fips[nchar(fips) < 5])

# reads mapping of FIPS codes to county names
conn <- url(paste0("http://www2.census.gov/geo/docs/reference/codes/",
                   "files/national_county.txt"))
fips.original <- read.csv(conn, header = FALSE, 
                          colClasses = rep("character", 5))

fips.original$FIPS <- paste0(fips.original$V2,fips.original$V3)

# list of TVA states
tva.states <- c("AL", "GA", "KY", "MS", "NC", "TN", "VA")

# rows with counties in above states
idx.rows <- which(fips.original$V1 %in% tva.states)

# keeps only these counties
fips.tva <- fips.original[idx.rows, ]

# which rows in data.nca contain the counties from TVA
idx.rows.2 <- which(fips %in% fips.tva$FIPS)

# get idx of columns with annual data
cols.data <- grep("Y2", names(data.nca$dbf))

# filter rows and columns with population data from TVA
data.pop.tva <- data.nca$dbf[idx.rows.2, c(1, cols.data)]


