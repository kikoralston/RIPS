
table.markdown <- function(df, name.file=''){
  # writes data frame in markdown table format to a txt file
  
  df$mean.load <- round(df$mean.load, 0)
  df$name <- as.character(df$name)
  nc <- ncol(df)
  
  cat('| ', paste(names(df), collapse =" | "), ' |\n',
      '|', rep('----- |', nc), '\n', sep='', file = name.file)
  for (i in 1:nrow(df)) {
    cat('|', paste(df[i, ], collapse =" | "), '| \n', file = name.file)
  }
}

library('plyr')
library('dplyr')
library(lubridate)
library(tidyr)

library(ggplot2)
library(gridExtra)
library(grid)

companies <- read.csv(file = paste('/Users/kiko/Google Drive/CMU/RIPS/SERC companies/SERC.csv'),
                      stringsAsFactors = FALSE)
companies[companies$Name.ferc != "", c("Name")] <- companies[companies$Name.ferc != "", c("Name.ferc")]
companies <- unlist(companies %>% select(Name, FERC) %>% filter(FERC == 1) %>% select(Name))

source('readFERC.R')

ferc.zip <- downloadFERCdata(ferc.zip = 'form714-database.zip')

df.results <- data.frame(name=as.character(),
                         year.ini=as.numeric(),
                         year.end=as.numeric(),
                         n.points=as.numeric(),
                         mean.load=as.numeric())

pb <- txtProgressBar(min = 0, max = length(companies), style = 3)
setTxtProgressBar(pb, 0)

data.ferc <- NULL

for (i in 1:length(companies)) {
  
  aux <- readFercNew(get.online = FALSE,
                     ferc.zip = ferc.zip,
                     name.util = companies[i])
  aux$name <- companies[i]
  
  if(i==1){
    data.ferc <- aux
  } else{
    data.ferc <- rbind(data.ferc, aux)
  }
  
  setTxtProgressBar(pb, i)
}

df.summary <- data.ferc %>% group_by(name) %>% 
  summarize(year.ini=min(year(time)), year.end=max(year(time)),
            n.points=n(), mean.load=mean(load, na.rm = TRUE))

row.names(df.summary) <- NULL
table.markdown(df = df.summary)

xx <- data.ferc %>% group_by(time) %>% 
  summarise(hourly.load=sum(load, na.rm=TRUE)) %>%
  mutate(year=year(time))

annual.summary <- data.ferc %>% mutate(year=year(time)) %>% group_by(year) %>%
  summarize(annual.energy=sum(load, na.rm=TRUE)/1e3,
            n.companies=length(unique(name)))

zz <- left_join(xx, annual.summary, by='year')

# Read EIA data ----
source("excel.R")

cat('\n\nReading EIA data ...\n')
eia.url <- 'https://www.eia.gov/electricity/data/eia411/archive/net_energy_load_2015.xls'

eia.file <- tempfile()
download.file(url = eia.url, destfile = eia.file, method='curl')

f.csv <- tempfile(fileext = ".csv")
f.xls <- eia.file

convert.xls2csv(xls.file = f.xls, csv.file = f.csv)

list.csv <- list.files(path=tempdir(), pattern = ".csv", 
                       full.names = TRUE)
f.csv <- list.csv[1]
eia.data <- read.table(file=f.csv, sep = ",", 
                       skip=8, nrows=20, header = TRUE,
                       stringsAsFactors = FALSE)
eia.data <- eia.data[ ,c(-1)]

row.serc <- grep("SERC", eia.data[,1], ignore.case = TRUE)

eia.aux <- eia.data[row.serc, ]
eia.aux <- eia.aux[-1]
eia.aux <- data.frame(year=seq(1990, 2026), energy=unlist(eia.aux))

a <- as.POSIXct(paste0(eia.aux$year,"-01-01"), format="%Y-%m-%d")
b <- as.POSIXct(paste0(eia.aux$year,"-12-31"), format="%Y-%m-%d")

eia.data.serc <- rbind(data.frame(time=a, energy=eia.aux$energy),
                       data.frame(time=b, energy=eia.aux$energy))

eia.data.serc <- eia.data.serc[order(eia.data.serc$time), ]

eia.data.serc <- eia.data.serc %>% filter(year(time) %in% seq(2006, 2015))

eia.data.serc <- eia.data.serc %>% left_join(zz, by='time') %>% 
  select(time, energy, annual.energy)
names(eia.data.serc) <- c('time', 'EIA', 'FERC.Data')

eia.data.serc <- eia.data.serc %>% gather(key='type', value='energy', EIA:FERC.Data)

g1 <- ggplot(data=zz) + geom_line(aes(x=time, y=hourly.load/1e3)) + 
  ylab("SERC Hourly Load (GW)") +
  xlab("")

g2 <- ggplot() +
  geom_line(data=eia.data.serc, aes(x=time, y=energy/1e3, linetype=type)) +
  ylab(expression(paste("SERC Annual Consumption(", 10^3, " MWh)"))) +
  xlab("") + theme(legend.position=c(1,1), legend.justification=c(1,1))


g3 <- ggplot(data=zz) + geom_line(aes(x=time, y=n.companies)) + 
  ylab("# companies in data set") +
  xlab("time") +
  theme(axis.title.y=element_text(margin=margin(0,20,0,0)))

g1 <- ggplot_gtable(ggplot_build(g1))
g2 <- ggplot_gtable(ggplot_build(g2))
g3 <- ggplot_gtable(ggplot_build(g3))

maxWidth <- unit.pmax(g1$widths[2:3], g2$widths[2:3], g3$widths[2:3])

g1$widths[2:3] <- maxWidth
g2$widths[2:3] <- maxWidth
g3$widths[2:3] <- maxWidth

png(paste0(path.plot, 'serc_data.png'), width=11, height=8.5, units='in', 
    res=300 )
grid.arrange(g1, g2, g3)
dev.off()

