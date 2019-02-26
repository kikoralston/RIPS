# This script computes the initial CO2 emission values from the results of
# a UCED run using 2015 fleet (or other current fleet). It reads two files 
# (genByPlantUCYYY.csv and genFleetUCYYY.csv) and prints out the total
# CO2 emissions for SERC that will be used by the CE model

library(ssh)
library(plyr)
library(dplyr)
library(jsonlite)

# if TRUE will try to dowanload files from a server using SSH
download.files <- FALSE

# path with complete UC results (assumes it is in a server)
path.results <- '/sharedstorage/user/fralston/out/CE_RUNS/AreacoreSERCCellsallCurtailEnvRegsCnoneSnor/UC/'

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
gen.by.plant <- read_csv(paste0(path.out, c('/genByPlantUC2015.csv'))) %>%
  gather(key='hour', value='gen', -genID)

gen.fleet <- read.csv(paste0(path.out, c('/genFleetUC2015.csv')), 
                      stringsAsFactors=FALSE) %>%
  mutate(genID = paste0(ORIS.Plant.Code, '+', Unit.ID)) %>%
  select(genID, CO2EmRate.lb.MMBtu., Heat.Rate..Btu.kWh.) %>%
  mutate(CO2EmRate.ton.MWh.=(CO2EmRate.lb.MMBtu./2000)*(Heat.Rate..Btu.kWh./1000))

df.gen.total <- gen.by.plant %>% left_join(gen.fleet, by= 'genID') %>%
  mutate(CO2Em.ton=gen*1000*CO2EmRate.ton.MWh.) %>%
  group_by(genID) %>% 
  summarise(gen=sum(gen), CO2Em.ton=sum(CO2Em.ton))

total.emissions <- signif(sum(df.gen.total$CO2Em.ton), 6) + 1000

cat('----------------------------------------------------------------\n',
    sprintf('TOTAL EMISSIONS IN 2015: %s tons CO2\n', 
            formatC(total.emissions, format= 'd', big.mark = ',')),
    '----------------------------------------------------------------\n',
    sep='')

write_json(data.frame(startYr=2015, startEms=total.emissions), 
           path=paste0(path.out, '/co2values.json'))


