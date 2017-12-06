# ************************************************************************
# Read Data from Shared Socioeconomic Pathways
# ************************************************************************

library(tidyr)
library(ggplot2)

read.ssp <- function() {
  url.ssp <- paste0("https://tntcat.iiasa.ac.at/SspDb/download/",
                    "common/SspDb_country_data_2013-06-12.csv.zip")
  
  csv.file <- basename(url.ssp)
  csv.file <- substr(csv.file, 1, nchar(csv.file) - 4)
  
  tt <- tempfile()
  download.file(url.ssp, tt)
  ssp.data <- read.csv(unz(tt, csv.file))
  unlink(tt)
  
  pop <- 
    ssp.data[(ssp.data$VARIABLE == "Population" & ssp.data$REGION == "USA"), ]
  
  gdp <- 
    ssp.data[(ssp.data$VARIABLE == "GDP|PPP" & ssp.data$REGION == "USA"), ]
  
  pop.2 <- pop[pop$MODEL == "IIASA GDP", ]
  # The arguments to gather():
  # - data: Data object
  # - key: Name of new key column (made from names of data columns)
  # - value: Name of new value column
  # - ...: Names of source columns that contain values
  pop.2 <- gather(pop.2, YEAR, pop, X1950:X2150)
  pop.2$YEAR <- as.numeric(gsub("X", "", pop.2$YEAR))
  
  gdp.2 <- gdp[gdp$MODEL == "IIASA GDP", ]
  gdp.2 <- gather(gdp.2, YEAR, gdp, X1950:X2150)
  gdp.2$YEAR <- as.numeric(gsub("X", "", gdp.2$YEAR))
  
  if(FALSE) {
    g <- ggplot()
    g <- g + geom_line(data = pop.2, aes(x=YEAR, y=pop, colour=SCENARIO))
    g <- g + geom_point(data = pop.2, aes(x=YEAR, y=pop, colour=SCENARIO))
    g <- g + theme_bw() + ylab("US Population (Million)")
    png("~/GoogleDrive/CMU/RIPS/climatechange/slides_scenarios/pop_plot.png",
        width = 960, height = 600, res = 100)
    print(g)
    dev.off()
    
    g <- ggplot()
    g <- g + geom_line(data = gdp.2, aes(x=YEAR, y=gdp, colour=SCENARIO))
    g <- g + geom_point(data = gdp.2, aes(x=YEAR, y=gdp, colour=SCENARIO))
    g <- g + theme_bw() + ylab("US Real GDP (Billion 2005 US$)")
    png("~/GoogleDrive/CMU/RIPS/climatechange/slides_scenarios/gdp_plot.png",
        width = 960, height = 600, res = 100)
    print(g)
    dev.off()
  }
  
  return(list(pop.2, gdp.2))
  
}
