# This script computes the initial CO2 emission values from the results of
# a UCED run using 2015 fleet (or other current fleet). It reads two files 
# (genByPlantUCYYY.csv and genFleetUCYYY.csv) and prints out the total
# CO2 emissions for SERC that will be used by the CE model

library(ssh)
library(plyr)
library(dplyr)
library(jsonlite)
library(plotly)
library(tidyr)

# if TRUE will try to dowanload files from a server using SSH
download.files <- TRUE

# path with complete UC results (assumes it is in a server)
path.results <- '/sharedstorage/user/fralston/out/CE_RUNS/UC/UC/'

# path in local computer to save remote files
path.out <- '~/Documents/CE/out/'

if (download.files) {
  ssh.session <- ssh_connect(host='fralston@epp-htc-005.ece.local.cmu.edu')
  
  scp_download(ssh.session, 
               files=paste0(path.results, c('/genByPlantUC2015.csv',
                                            '/genFleetUC2015.csv'), collapse = ' '), 
               to=path.out)
  
  ssh_disconnect(ssh.session)  
}

# read UC files and copmute total CO2 emissions
gen.by.plant <- read.csv(paste0(path.out, c('/genByPlantUC2015.csv')), 
                         stringsAsFactors = FALSE) %>%
  gather(key='hour', value='gen', -genID)

gen.fleet <- read.csv(paste0(path.out, c('/genFleetUC2015.csv')), 
                      stringsAsFactors=FALSE) %>%
  mutate(genID = paste0(ORIS.Plant.Code, '+', Unit.ID)) %>%
  select(Plant.Name, ORIS.Plant.Code, genID, CO2EmRate.lb.MMBtu., 
         Heat.Rate..Btu.kWh., Capacity..MW.,  PlantType, Modeled.Fuels, 
         State.Name) %>%
  mutate(CO2EmRate.ton.MMBTU=CO2EmRate.lb.MMBtu./2000) %>%
  mutate(CO2EmRate.ton.MWh.=(CO2EmRate.ton.MMBTU)*(Heat.Rate..Btu.kWh./1000))

df.gen.total.uc <- gen.by.plant %>% left_join(gen.fleet, by= 'genID') %>%
  mutate(CO2Em.ton=gen*1000*CO2EmRate.ton.MWh.) %>%
  group_by(genID) %>% 
  summarise(gen=sum(gen), CO2Em.ton=sum(CO2Em.ton), Capacity=Capacity..MW.[1],
            PlantType=PlantType[1], Modeled.Fuels=Modeled.Fuels[1])

total.emissions <- signif(sum(df.gen.total.uc$CO2Em.ton), 6) + 1000
sum(df.gen.total.uc$gen)

cat('----------------------------------------------------------------\n',
    sprintf('TOTAL EMISSIONS IN 2015: %s tons CO2\n', 
            formatC(total.emissions, format= 'd', big.mark = ',')),
    '----------------------------------------------------------------\n',
    sep='')

write_json(data.frame(startYr=2015, startEms=total.emissions), 
           path=paste0(path.out, '/co2values.json'))

uc.generation <- df.gen.total.uc %>% 
  mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
  mutate(PlantType = ifelse(Modeled.Fuels %in% coal.types, 'Coal', PlantType)) %>%
  mutate(PlantType = ifelse(grepl('Natural Gas', Modeled.Fuels), 
                            'Natural Gas', PlantType)) %>%
  mutate(PlantType = ifelse(grepl('Fuel Oil', Modeled.Fuels), 
                            'Fuel Oil', PlantType)) %>%
  mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 
                            'Hydro', PlantType)) %>%
  mutate(PlantType = ifelse(PlantType == 'Solar PV', 'Solar', PlantType)) %>%
  select(-Modeled.Fuels) %>% group_by(PlantType) %>%
  summarise(gen=sum(gen, na.rm = TRUE)) %>% mutate(type='uc')

# compare uc generation in 2015 with EIA historical generation

# df.total comes from file '/RIPS/EIA_Data/get_power_plant_gen.R'
eia.gen <- df.total %>% group_by(year, type) %>% 
  summarise(gen=sum(gen, na.rm = TRUE)/1000) %>%
  filter(year==2015) %>% mutate(PlantType=type) %>%
  mutate(type='EIA') %>% ungroup() %>% select(PlantType, gen, type)

df.compare <- rbind(uc.generation, eia.gen)

ggplotly(ggplot(df.compare) + geom_col(aes(x=type, y=gen, fill=PlantType)) + 
           scale_fill_manual(values=col.pallete))

