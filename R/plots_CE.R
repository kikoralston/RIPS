# Jannuary 2019
# Set of functions that read CE model output and create plots for the analysis

library(rgeos)
library(plyr)
library(dplyr)
library(ggplot2)
library(ssh)
library(tidyr)
library(rgdal)
library(maps)
library(gdxrrw)
library(xtable)
library(lubridate)
library(maptools)
library(gridExtra)
library(grid)
library(mapdata)
library(sf)


# General Parameters ----

renewablePlantTypes <- c('Geothermal', 'Hydro', 'Pumped Storage', 'Wind', 
                        'Solar PV')
otherTypes <- c('Geothermal', 'Fuel Cell', 'Landfill Gas', 'Non-Fossil Waste',
                'Municipal Solid Waste', 'Biomass')
ngTypes <- c('Combined Cycle', 'Combustion Turbine', 'O/G Steam')

plant.types <- c('Hydro', 'Nuclear', 'Coal', 'Natural Gas', 'Solar', 'Wind', 
                 'Fuel Oil', 'Other')

coal.types <- c('Bituminous', 'Bituminous& Subbituminous', 'Lignite', 
                'Subbituminous')

col.pallete <- rev(c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1',
                     '#35978f','#01665e'))

shapes.pallete <- c(21, 22, 23, 24, 25, 15, 16, 17, 18, 19)

# Auxiliary Functions ----

set.year.retired.by.age <- function(year.online, lifetime, retirement.yr) {
  
  year.retired <- ifelse(retirement.yr == 9999,
                         year.online + lifetime,
                         retirement.yr)
  
  return(year.retired)
}

update.ID <- function(uniqueid, planttype, region, capacity) {
  
  ifelse(uniqueid == "",
         paste(uniqueid, gsub(" ", "",as.character(planttype)), 
               region, 
               as.character(capacity),  sep = "_"),
         uniqueid) 
  
}

parse.cooling.type <- function(x) {
  
  idx.cooling <- grepl('\\+OT|\\+DC|\\+RC', x, perl=TRUE)
  
  out <- x
  
  out <- gsub('\\.', ' ', out)
  
  out[idx.cooling] <- substr(out[idx.cooling], nchar(out[idx.cooling])-1, 
                             nchar(out[idx.cooling]))
  out[!idx.cooling] <- NA
  
  return(out)
}

parse.type <- function(x) {
  
  out <- x
  
  # remove digits with Wind/solar aggregation block
  out <- gsub('\\.\\d\\d', '', out)
  
  out <- gsub('\\+', ' ', out)
  out <- gsub('OT|DC|RC|CCS', '', out, perl = TRUE)
  out <- trimws(out)
  return(out)
}

parse.ccs <- function(x) {
  idx.ccs <- grepl('\\.CCS', x, perl=TRUE)
  return(idx.ccs)  
}

get.capacity.tech <- function(type, ccs, cool, df.new.techs) {
  
  out <- paste(type, ifelse(ccs, 'CCS', ''))
  
  out <- trimws(out)
  
  out <- ifelse(is.na(cool), out, paste(out, cool))
  
  new.techs <- df.new.techs %>% 
    mutate(Cooling = ifelse(Cooling == "dry cooling", 'DC',
                            ifelse(Cooling == "once through",
                                   'OT',
                                   ifelse(Cooling == "recirculating",
                                          'RC',
                                          NA)))) %>%
    mutate(Type = ifelse(is.na(Cooling), as.character(Type), 
                         paste(Type, Cooling)))
  
  out <- data.frame(Type=out, stringsAsFactors = FALSE) %>%
    left_join(new.techs, by='Type') %>% as.data.frame() %>%
    select(Capacity) %>% unlist(use.names = FALSE)
  
  return(out)
}

# Main Functions ----

plot.map.serc <- function(path.data.ce, path.plot){
  
  g <- get.map.serc(path.data.ce)
  
  pdf(paste0(path.plot, "/sercmap.pdf"), width=7*16/9, height=7)
  print(g)
  dev.off()
  
  system(sprintf('pdfcrop %s %s', paste0(path.plot, "/sercmap.pdf"), 
                 paste0(path.plot, "/sercmap.pdf")))
}

plot.map.serc.with.plants <- function(path.data.ce, path.plot, gen.fleet, 
                                      show.points=TRUE, width=7, 
                                      axis.labels=TRUE){
  
  g <- get.map.serc(path.data.ce)

  gen.fleet.2 <- gen.fleet %>% 
    select(Plant.Name, UniqueID_Final, PlantType, Region.Name, State.Name, County, 
           Capacity..MW., On.Line.Year, Retirement.Year,
           YearRetiredByAge, Modeled.Fuels, Longitude, Latitude) %>%
    mutate(PlantType = as.character(PlantType), Modeled.Fuels=as.character(Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
    # simplify fossil fuels in Modeled.Fuels
    mutate(Modeled.Fuels = ifelse(Modeled.Fuels %in% coal.types, 'Coal', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Natural Gas', Modeled.Fuels), 'Natural Gas', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Fuel Oil', Modeled.Fuels), 'Fuel Oil', Modeled.Fuels)) %>%
    # for fossil fuels, replace fuel in PlantType
    mutate(PlantType = ifelse(Modeled.Fuels %in% c("Coal", "Natural Gas", "Fuel Oil"), Modeled.Fuels, PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 'Hydro', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType == 'Solar PV', 'Solar', PlantType)) %>%
    mutate(Plant.Name = ifelse(Plant.Name == '', paste0(Longitude, '_', Latitude), 
                               Plant.Name)) %>%
    filter(Plant.Name != '' & Plant.Name != 'NA_NA') %>%
    select(-Modeled.Fuels) %>% group_by(Plant.Name) %>%
    dplyr::summarise(PlantType=PlantType[1], Capacity..MW.=sum(Capacity..MW., na.rm = TRUE),
                     Longitude=Longitude[1], Latitude=Latitude[1])

  if (show.points) {
    g_plants <- g + 
      geom_point(data=gen.fleet.2, aes(x=Longitude, y=Latitude, shape=PlantType), 
                 fill='yellow')
  } else {
    g_plants <- g + 
      geom_point(data=gen.fleet.2, aes(x=Longitude, y=Latitude, shape=PlantType), 
                 fill='yellow', alpha=0.0)
  } 
  
  g_plants <- g_plants +
    scale_shape_manual(values = shapes.pallete) + 
    theme(legend.position = "top", legend.title=element_blank(),
          legend.margin=margin(), legend.spacing=unit(0, 'points'),
          plot.margin=margin(),legend.box.spacing=unit(0.1, 'points'),
          legend.key.height = unit(0.8, 'line')) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
  if (!show.points) {
    g_plants <- g_plants + theme(legend.text = element_text(colour = 'white'))
  }

  if (!axis.labels) {
    g_plants <- g_plants + theme(axis.text=element_blank())
  }
  
  # assume ratio 16:9
  pdf(path.plot, width=width, height=width*9/16)
  print(g_plants)
  dev.off()
  system(sprintf('pdfcrop %s %s', path.plot, path.plot))
  
  return(g_plants)
}


get.map.serc <- function(path.data.ce) {
  
  bounds.longitude <- c(-91, -76)
  bounds.latitude <- c(29, 39)
  
  serc.zones <- c('S_SOU', 'S_C_TVA', 'S_VACA', 'S_C_KY')
  
#  s1 <- readOGR(paste0(path.data.ce, 'IPMRegionsShapeFile/'), 
#                'IPM_Regions_20121023_US')
  
  s1 <- st_read(paste0(path.data.ce, 'IPMRegionsShapeFile/'), 
                'IPM_Regions_20121023_US') 
  
  zones.map <- ifelse(s1$IPM_Region %in% serc.zones, 
                      s1$IPM_Region, NA)
  
  # s1.union.serc <- unionSpatialPolygons(s1, zones.map, avoidUnaryUnion=TRUE)
  
  s1.union.serc <- s1 %>% filter(IPM_Region %in% serc.zones) %>% filter(st_is_valid(.)) 
  
  l <- list()
  for (z in serc.zones) {
    l[[z]] <- s1.union.serc %>% filter(IPM_Region == z) %>% st_union() %>% st_simplify()
  }

  df.ipm.zones <- l %>% ldply(.fun = data.frame) %>% st_sf() %>% st_simplify()
  
  df.ipm.counties <- s1 %>% mutate(FIPS=as.numeric(as.character(FIPS)))
  
  map.us.counties <- map_data('county')
  
  map.us.counties <- map.us.counties %>% 
    mutate(polyname=paste(region, subregion, sep=',')) %>%
    left_join(county.fips, by='polyname') %>%
    dplyr::rename(FIPS=fips) %>% left_join(df.ipm.counties, by='FIPS')
  
  map.us.counties <- map.us.counties %>% 
    mutate(IPM_Region = as.character(IPM_Region)) %>%
    mutate(IPM_Region = ifelse(IPM_Region %in% serc.zones,
                               IPM_Region,
                               NA))
  
  # state boundaries & rivers
  state_map <- map_data("state")
  path.rivers <- "~/CMU/RIPS/paperCE/rivers/"
  rivers <- readOGR(path.rivers)
  rivers2 <- subset(rivers, COUNTRY == "USA")
  rivers3 <- st_transform(st_as_sf(rivers2), 
                          crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  rivers4 <- st_crop(x=rivers3, 
                     c(xmin=bounds.longitude[1], xmax=bounds.latitude[2], 
                       ymin=bounds.latitude[1], ymax=bounds.latitude[2]))
  
  labels_regions <- st_centroid(df.ipm.zones$geometry) %>% st_coordinates() %>% as.data.frame() %>%
    dplyr::rename(x=X, y=Y)
  labels_regions$label <- as.character(df.ipm.zones$.id)
  
#  labels_regions <- data.frame(x=as.numeric(),
#                               y=as.numeric(),
#                               label=as.character())
#  for (p in s1.union.serc@polygons) {
#    x <- p@labpt[1]
#    y <- p@labpt[2]
#    label <- levels(df.ipm.counties$IPM_Region)[as.integer(p@ID)]
#    
#    labels_regions <- rbind(labels_regions,
#                            data.frame(x=x,
#                                       y=y,
#                                       label=label))
#  }
  
  pal <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
           '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  
  # change position of labels and jsutification
  labels_regions <- labels_regions %>%
    mutate(hjust=0.5, vjust=0.5)
  
  labels_regions <- labels_regions %>%
    mutate(hjust=ifelse(label=='S_C_KY', 0, hjust),
           hjust=ifelse(label=='S_C_TVA', 0.75, hjust),
           hjust=ifelse(label=='S_VACA', 0, hjust),
           vjust=ifelse(label=='S_SOU', 0.8, vjust)) %>%
    mutate(x=ifelse(label=='S_VACA', -78, x))
  
  g <- ggplot() +
    geom_polygon(data=map.us.counties, aes(x=long, y=lat, group=group), 
                 colour='gray80') +
    geom_polygon(data=map.us.counties, aes(x=long, y=lat, group=group, 
                                           fill=IPM_Region)) +
    geom_sf(data=rivers4, colour='dodgerblue') +
    coord_sf(xlim=bounds.longitude, ylim = bounds.latitude) +
    geom_label(data=labels_regions, aes(x=x, y=y, label=label, hjust=hjust, 
                                        vjust=vjust), size=rel(2.3)) + 
    guides(fill=FALSE) +
#    coord_map("polyconic", xlim=bounds.longitude, ylim=bounds.latitude) + 
    scale_fill_manual(breaks=serc.zones, values=pal, na.value="gray95") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + 
    theme_bw()
  
  return(g)
}

plot.inital.fleet <- function(path.output.ce, path.data.ce, df.demand=NULL, 
                              path.out.plot=path.output.ce, width=7,
                              divide.by.zone=FALSE) {
  
  ini.gen.fleet.raw <- read.csv(file=paste0(path.output.ce, '/genFleetInitial.csv'))
  
  lifetimes <- read.csv(file=paste0(path.data.ce, 
                                    'NewPlantData/LifetimeValuesExistingPlants4Aug2016.csv'))
  
  ini.gen.fleet.2 <- ini.gen.fleet.raw %>% 
    left_join(lifetimes, by='PlantType') %>%
    select(Plant.Name, UniqueID_Final, PlantType, Region.Name, State.Name, County, 
           Capacity..MW., On.Line.Year, Retirement.Year, Lifetime.yrs., 
           YearRetiredByAge, Modeled.Fuels) %>%
    mutate(YearRetiredByAge=set.year.retired.by.age(On.Line.Year, 
                                                    Lifetime.yrs.,
                                                    Retirement.Year)) %>%
    mutate(YearRetiredByAge=ifelse(PlantType %in% renewablePlantTypes,
                                   9999,
                                   YearRetiredByAge)) %>%
    mutate(UniqueID_Final = as.character(UniqueID_Final)) %>%
    mutate(UniqueID_Final = update.ID(UniqueID_Final, PlantType, Region.Name, 
                                      Capacity..MW.)) %>%
    mutate(PlantType = as.character(PlantType), Modeled.Fuels=as.character(Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
    # simplify fossil fuels in Modeled.Fuels
    mutate(Modeled.Fuels = ifelse(Modeled.Fuels %in% coal.types, 'Coal', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Natural Gas', Modeled.Fuels), 'Natural Gas', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Fuel Oil', Modeled.Fuels), 'Fuel Oil', Modeled.Fuels)) %>%
    # for fossil fuels, replace fuel in PlantType
    mutate(PlantType = ifelse(Modeled.Fuels %in% c("Coal", "Natural Gas", "Fuel Oil"), Modeled.Fuels, PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 'Hydro', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType == 'Solar PV', 'Solar', PlantType))
  
  years <- expand.grid(UniqueID_Final=unique(ini.gen.fleet.2$UniqueID_Final), 
                       currYear=seq(2020, 2050, by=5))
  
  df.fleet.time <- ini.gen.fleet.2 %>% left_join(years, by='UniqueID_Final') %>% 
    arrange(UniqueID_Final)
  
  # create new column with online Cap
  
  df.fleet.time <- df.fleet.time %>% 
    mutate(online.cap = ifelse(YearRetiredByAge < currYear, 0, Capacity..MW.))
  
  # aggreagate by PlantType and region and currYear
  df.fleet.time.2 <- df.fleet.time %>% group_by(Region.Name, PlantType, currYear) %>%
    dplyr::summarise(online.cap = sum(online.cap)) %>% arrange(currYear) %>%
    ungroup() %>%
    mutate(PlantType=factor(PlantType, levels=rev(plant.types))) %>%
    mutate(Region.Name = as.character(Region.Name))
  
  list.cases <- NULL
  df.demand.aux <- NULL
  if(!is.null(df.demand)) {
    df.demand.aux <- df.demand %>% gather(key = "type", value = "value", 
                                          mean.load:peak.load)
    list.cases <- unique(df.demand$case)
  }
  ncases <- length(list.cases)
  
  g <- ggplot() + geom_col(data=df.fleet.time.2,
                           aes(x=currYear, y=online.cap/1e3, fill=PlantType),
                           width=2)
  
  g <- g + scale_fill_manual(values = col.pallete) + theme_bw() +
    xlab('Year') + ylab('Installed Capacity (GW)') +
    scale_x_continuous(breaks = seq(2020, 2050, 5)) +
    guides(fill=guide_legend(title=NULL))
    
  if (!is.null(df.demand.aux)) {
    g <- g + geom_line(data=df.demand.aux, 
                       aes(x=currYear, y=value/1e3, linetype=case, color=type), 
                       size=1) + 
      guides(linetype=guide_legend(title=NULL, reverse=TRUE))
  }
  
  if (divide.by.zone) {
    g <- g + facet_grid(Region.Name ~ ., scales = 'free_y') 
  }
  
  pdf(paste0(path.out.plot,"/supply_demand.pdf"), width=width, height=width*9/16)
  print(g)
  dev.off()
  
  return(df.fleet.time.2)
}


get.fleet <- function(fileGenfleet, currYear) {
  
  df.fleet <- read.csv(file=fileGenfleet, 
                       stringsAsFactors = FALSE) %>%
    # remove plants retired by age
    filter(YearRetiredByAge > currYear | is.na(YearRetiredByAge)) %>%
    # remove plants originally marked to retire in IPM
    filter(Retirement.Year > currYear) %>%
    # remove plants that were retired by the CE model due to low CF
    filter(is.na(YearRetiredByCE)) %>%
    mutate(PlantType = as.character(PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
    mutate(Modeled.Fuels = ifelse(PlantType == 'Other', 'Other', Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(Modeled.Fuels == 'Other', 'NA', PlantType)) %>%
    mutate(Modeled.Fuels = ifelse(Modeled.Fuels %in% coal.types, 'Coal', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(Modeled.Fuels == 'Nuclear Fuel', 'Nuclear', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Natural Gas', Modeled.Fuels), 'Natural Gas', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Fuel Oil', Modeled.Fuels), 'Fuel Oil', Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 
                              'Hydro', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType == 'Combustion Turbine', 'CT', PlantType)) %>%
    mutate(Modeled.Fuels = ifelse(PlantType == 'Hydro', 'Hydro', Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(grepl('Hydro|Wind|Solar', PlantType), 'NA', PlantType)) %>%
    mutate(cooling=ifelse(grepl('recirculating', Cooling.Tech),
                          'RC',
                          ifelse(grepl('once through', Cooling.Tech),
                                 'OT',
                                 ifelse(grepl('dry cooling', Cooling.Tech),
                                        'DC', 'NA')))) %>%
    mutate(cooling = ifelse(PlantType == 'Hydro', 'NA', cooling)) %>%
    mutate(cooling = ifelse(Modeled.Fuels == 'Other', 'NA', cooling)) %>%
    mutate(PlantType = ifelse(PlantType == 'NA', NA, PlantType)) %>%
    mutate(cooling = ifelse(cooling == 'NA', NA, cooling)) %>%
    mutate(id = paste0(ORIS.Plant.Code, "+", Unit.ID)) %>%
    select(Plant.Name, PlantType, Modeled.Fuels, Capacity..MW.,cooling, id)
  
  return(df.fleet)
  
}


get.new.additions <- function(y, path.gdx) {
  
  serc.zones <- read.csv(paste0(path.gdx,'/zoneNamesToNumbers.csv'))
  
  # map cell to zone
  map.cell2zone <- rgdx.param(gdxName = paste0(path.gdx, 'gdxOutYear', y, 
                                               '.gdx'), 'pCellzones')
  map.cell2zone$c <- as.character(map.cell2zone$c)
  
  tech.not.curtailed <- rgdx(gdxName = paste0(path.gdx, '/gdxOutYear', y, '.gdx'), 
                             requestList = list(name="vNnotcurtailed", 
                                                form='full', 
                                                field=c('l')))
  
  df.tech.not.curtailed <- reshape2::melt(tech.not.curtailed$val) %>%
    dplyr::rename(zone=z, type=technotcurtailed) %>% 
    mutate(ZoneNum=as.numeric(as.factor(zone))) %>% select(-zone) %>%
    left_join(serc.zones, by='ZoneNum') %>%
    arrange(ZoneNum)
  
  if(nrow(df.tech.not.curtailed) > 0) {
    df.tech.not.curtailed <- df.tech.not.curtailed %>% mutate(year=y)
  } else {
    df.tech.not.curtailed <- df.tech.not.curtailed %>% mutate(year=as.integer())
  }
  df.tech.not.curtailed <- df.tech.not.curtailed %>% 
    dplyr::select(year, Zone, type, value) %>% as.data.frame()
  
  tech.curtailed <- rgdx(gdxName = paste0(path.gdx, '/gdxOutYear', y, '.gdx'), 
                         requestList = list(name="vNcurtailed", 
                                            form='full',field=c('l')))
  df.tech.curtailed <- reshape2::melt(tech.curtailed$val) %>%
    mutate(c=as.character(c)) %>%
    dplyr::rename(type=techcurtailed) %>%
    left_join(map.cell2zone, by="c") %>%
    mutate(ZoneNum=pCellzones) %>% select(-pCellzones) %>%
    left_join(serc.zones, by='ZoneNum') %>%
    arrange(ZoneNum)

  if(nrow(df.tech.curtailed) > 0) {
    df.tech.curtailed <- df.tech.curtailed %>% mutate(year=y)
  } else {
    df.tech.curtailed <- df.tech.curtailed %>% mutate(year=as.integer())
  }
  df.tech.curtailed <- df.tech.curtailed %>% 
    dplyr::select(year, Zone, type, value) %>% as.data.frame()
  
  df.tech.curtailed.zone <- df.tech.curtailed %>% group_by(Zone, type) %>%
    dplyr::summarise(value=sum(value, na.rm=TRUE)) %>% mutate(year=y) %>%
    dplyr::select(year, Zone, type, value) %>% as.data.frame()
  
  tech.renew <- rgdx(gdxName = paste0(path.gdx, '/gdxOutYear', y, '.gdx'), 
                     requestList = list(name="vNrenew", 
                                        form='full',field=c('l')))
  
  df.tech.renew <- reshape2::melt(tech.renew$val) %>%
    dplyr::rename(zone=z, type=techrenew) %>% 
    mutate(ZoneNum=as.numeric(as.factor(zone))) %>% select(-zone) %>%
    mutate(type=ifelse(grepl('Wind', type), 'Wind', 'Solar PV')) %>%
    group_by(ZoneNum, type) %>% dplyr::summarize(value=sum(value)) %>% 
    ungroup() %>% left_join(serc.zones, by='ZoneNum') %>%
    arrange(ZoneNum)

  if(nrow(df.tech.renew) > 0) {
    df.tech.renew <- df.tech.renew %>% mutate(year=y)
  } else {
    df.tech.renew <- df.tech.renew %>% mutate(year=as.integer())
  }
  df.tech.renew <- df.tech.renew %>%  
    dplyr::select(year, Zone, type, value) %>% as.data.frame()
  
  df.complete.decision <- rbind(df.tech.curtailed.zone, df.tech.not.curtailed, 
                                df.tech.renew)
  
  return(df.complete.decision)
}

plot.new.additions <- function(path.output.ce,
                               gamsSysDir='/Applications/GAMS24.7/sysdir/',
                               path.rcps=c('RCP 4.5'=sprintf('%s/rcp45/',
                                                             path.output.ce),
                                           'RCP 8.5'=sprintf('%s/rcp85/',
                                                             path.output.ce)),
                               yearIni=2020, yearEnd=2050,
                               path.file=paste0(path.output.ce, "/newadditions.pdf"),
                               divide.zones=TRUE) {
  
  igdx(gamsSysDir)
  df.new.techs <- read.csv(paste0(path.rcps[1],
                                  sprintf('newTechsCE%4d.csv', yearIni)))
  #df.new.techs <- read.csv(file='./newtechs.csv', header = TRUE, stringsAsFactors = FALSE)
  #names(df.new.techs) <- df.new.techs[1,]
  #df.new.techs <- df.new.techs[-1,]
  
  df.new.techs <- df.new.techs %>% 
    transmute(Type=TechnologyType, Fuel=FuelType, Cooling=Cooling.Tech,
              Capacity=as.numeric(Capacity.MW.), 
              CAPEX=as.numeric(CAPEX.2012..MW.)/1e3,
              FixedOM=as.numeric(FOM.2012..MW.yr.)/1e3, 
              VarOM=as.numeric(VOM.2012..MWh.)) %>%
    arrange(Type, Cooling)
  
  # get investment decisions
  df.complete.decision <- data.frame(case=as.character(),
                                     year=as.numeric(), 
                                     Zone=as.character(),
                                     type=as.character(),
                                     value=as.numeric(),
                                     stringsAsFactors = FALSE)
  
  for (rcp in path.rcps) {
    #serc.zones <- read.csv(paste0(rcp,'zoneNamesToNumbers.csv'))
    case.name <- names(which(path.rcps == rcp))
    
    for (y in seq(yearIni, yearEnd, by=5)) {
      
      #path.gdx <- paste0(rcp, '/gdxOutYear', y, '.gdx')
      
      df <- get.new.additions(y, rcp)

      df.complete.decision <- rbind(df.complete.decision,
                                    data.frame(case=case.name, df, 
                                               stringsAsFactors = FALSE))
    }
  }
  
  # add columns for cooling type and CCS
  df.complete.decision.2 <- df.complete.decision %>%
    mutate(coolingType = parse.cooling.type(type),
           ccs = parse.ccs(type),
           Type=parse.type(type)) %>%
    mutate(capacity=get.capacity.tech(Type, ccs, coolingType, df.new.techs)) %>%
    select(case, year, Zone, Type, coolingType, ccs, value, capacity) %>%
    mutate(total.cap = value*capacity,
           Type = ifelse(Type == 'Combined Cycle', 'Natural Gas', Type)) %>%
    mutate(Type = ifelse(Type == 'Coal Steam', 'Coal', Type)) %>%
    mutate(Type = ifelse(Type == 'Solar PV', 'Solar', Type)) %>%
    mutate(Type=factor(Type, levels=rev(plant.types)))
  
  list.cases <- names(path.rcps) 
  ncases <- length(list.cases)
  
  if (ncases == 2) {
    width <- 1.5
    offset <- data.frame(case=list.cases,
                         off=c(-width/2, width/2))
  } else if (ncases == 3) {
    width <- 1
    offset <- data.frame(case=list.cases,
                         off=c(-1, 0, 1))
  }
  
  df.complete.decision.3 <- df.complete.decision.2 %>%
    left_join(offset, by='case') %>% mutate(year=year+off) %>%
    select(-off) %>%
    # remove CO2 from name of case (just for this plot --  to save space)
    mutate(case=trimws(gsub('CO2', '', case))) %>%
    mutate(case=ifelse(case=='Reference', 'Ref.', case))
    
  
  g <- ggplot() + geom_col(data=df.complete.decision.3, 
                           aes(x=year, y=total.cap/1e3, fill=Type), 
                           width=width-0.1)
  
  list.cases <- names(path.rcps) %>%
    gsub(pattern='CO2', replacement='', x=.) %>%
    trimws() %>%
    gsub(pattern='Reference', replacement='Ref.', x=.)
    
  breaks <- sort(unique(df.complete.decision.3$year))
  labels.up <- rep(list.cases, 7)
  
  g <- g + scale_fill_manual(limits = rev(plant.types), 
                             values = col.pallete) +
    theme_bw() + ylab('New Capacity (GW)') +
    guides(fill=guide_legend(title=NULL)) +
    theme_bw() +
    scale_x_continuous(sec.axis = dup_axis(breaks = breaks,
                                           labels = labels.up)) +
    theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"),
          axis.title.x = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0))
  
  if (divide.zones) {
    g <- g + facet_grid(Zone ~ ., scales = 'free_y')
  }
  
#   g_param <- ggplot_build(g)
#   x.range <- g_param$layout$panel_params[[1]]$x.range
#   year.values <- unique(df.complete.decision.3$year)
#   
#   df.top.axis <- df.complete.decision.3 %>%
#     mutate(xpos=year) %>%
#     select(case, xpos) %>%
#     distinct() %>%
#     # remove CO2 from name of case (just for this plot --  to save space)
#     mutate(case=trimws(gsub('CO2', '', case))) %>%
#     mutate(case=ifelse(case=='Reference', 'Ref.', case))
#     
#   g.top.axis <- ggplot(df.top.axis) + 
#     geom_text(aes(x=xpos, y=0, label=case), 
#               angle=90, hjust=0, vjust=0.5, size=3.3, 
#               nudge_y = 0.2) +
#     coord_cartesian(xlim = x.range, ylim = c(0, 1), expand = FALSE) +
#     scale_x_continuous(breaks = df.top.axis$xpos) +
#     theme(plot.margin = unit(c(0, 0, 0.15, 0), "cm"),
#           axis.title=element_blank(),
#           axis.text = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.ticks.length = unit(-0.25, "cm"),
#           axis.line.x = element_line(), 
#           panel.grid = element_blank(),
#           panel.background = element_rect(fill='white')
# #          panel.border = element_rect(fill = NA)
#     )
#   
#   lay <- rbind(c(rep(1,19), rep(2,329), rep(3, 52)),
#                c(rep(4,400)))
#   # 
#   gg <- arrangeGrob(rectGrob(gp=gpar(fill="white", lty='blank')),
#                     g.top.axis,
#                     #rectGrob(gp=gpar(fill="white", lty='blank')),
#                     rectGrob(gp=gpar(fill="white", lty='blank')),
#                     g, layout_matrix=lay,
#                     heights=c(2,20),
#                     padding=unit(0, units='mm'))
  
  # plot new additions by plant type
  pdf(path.file, width=7*16/9, height=7)
  print(g)
  #grid.draw(gg)
  dev.off()
  
  return(df.complete.decision.2)
}

plot.updated.fleet <- function(df.fleet, df.new.additions, path.output.ce,
                               df.demand,
                               gamsSysDir='/Applications/GAMS24.7/sysdir/',
                               path.rcps=c('RCP 4.5'=sprintf('%s/rcp45/',
                                                             path.output.ce),
                                           'RCP 8.5'=sprintf('%s/rcp85/',
                                                             path.output.ce)),
                               yearIni=2020, yearEnd=2050,
                               path.out.plot=path.output.ce) {
  
  list.cases <- names(path.rcps)
  ncases <- length(list.cases)
  
  for (i in 1:ncases) {
    if (i == 1) {
      df.fleet.updated <- data.frame(case=list.cases[i], df.fleet, 
                                     stringsAsFactors = FALSE)
    } else {
      df.fleet.updated <- rbind(df.fleet.updated,
                                data.frame(case=list.cases[i], df.fleet, 
                                           stringsAsFactors = FALSE))
    }
  }
  
  df.fleet.updated <- df.fleet.updated %>% 
    select(case, Region.Name, PlantType, currYear, online.cap)
  
  df.complete.decision.3 <- df.new.additions %>% 
    dplyr::rename(online.cap=total.cap, PlantType=Type, currYear=year, Region.Name=Zone) %>%
    select(case, Region.Name, PlantType, currYear, online.cap) %>%
    group_by(case, Region.Name, PlantType, currYear) %>%
    dplyr::summarise(online.cap=sum(online.cap))
  
  df.complete.decision.3 <- df.complete.decision.3 %>%
    arrange(case, Region.Name, PlantType, currYear) %>%
    group_by(case, Region.Name, PlantType) %>%
    mutate(online.cap = cumsum(online.cap)) 
  
  df.complete.decision.3 <- df.complete.decision.3 %>% as.data.frame()
  
  df.fleet.updated <- rbind(df.fleet.updated, df.complete.decision.3)  

  if (ncases == 2) {
    width <- 1.5
    offset <- data_frame(case=list.cases,
                         off=c(-width/2, width/2))
  } else if (ncases == 3) {
    width <- 1
    offset <- data_frame(case=list.cases,
                         off=c(-1, 0, 1))
  }
  
  df.aux1 <- df.fleet.updated %>%
    left_join(offset, by='case') %>% mutate(currYear=currYear+off) %>%
    select(-off)
  
  df.aux2 <- df.demand %>%
    left_join(offset, by='case') %>% mutate(currYear=currYear+off) %>%
    select(-off, -mean.load)
  
  g <- ggplot() + geom_col(data=df.aux1, 
                           aes(x=currYear, y=online.cap/1e3, fill=PlantType), 
                           width=width-0.1)
  
  g <- g + geom_point(data=df.aux2, 
                     aes(x=currYear, y=peak.load/1e3), shape=24, colour='black', 
                     fill='white', size=2)
  
  g <- g + scale_fill_manual(limits = rev(plant.types), values = col.pallete) +
    theme_bw() + ylab('Capacity (GW)') +
    guides(fill=guide_legend(title=NULL),
           shape=guide_legend(title=NULL)) + theme_bw() +
    theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "lines"),
          axis.title.x = element_blank()) +
    facet_grid(Region.Name ~ ., scales = 'free_y')

  g_param <- ggplot_build(g)
  x.range <- g_param$layout$panel_params[[1]]$x.range
  year.values <- unique(df.complete.decision.3$year)
  
  df.top.axis <- df.aux1 %>%
    mutate(xpos=currYear) %>%
    select(case, xpos) %>%
    distinct()
  
  g.top.axis <- ggplot(df.top.axis) + 
    geom_text(aes(x=xpos, y=0, label=case), 
              angle=90, hjust=0, vjust=0.5, size=3.3, 
              nudge_y = 0.2) +
    coord_cartesian(xlim = x.range, ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = df.top.axis$xpos) +
    theme(plot.margin = unit(c(0, 0, 0.15, 0), "cm"),
          axis.title=element_blank(),
          axis.text = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(-0.25, "cm"),
          axis.line.x = element_line(), 
          panel.grid = element_blank(),
          panel.background = element_rect(fill='white')
          #panel.border = element_rect(fill = NA)
    )
  
  lay <- rbind(c(rep(1,9), rep(2,165), rep(3, 26)),
               c(rep(4,200)))
  # 
  gg <- arrangeGrob(rectGrob(gp=gpar(fill="white", lty='blank')),
                    g.top.axis,
                    #rectGrob(gp=gpar(fill="white", lty='blank')),
                    rectGrob(gp=gpar(fill="white", lty='blank')),
                    g, layout_matrix=lay,
                    heights=c(2,20),
                    padding=unit(0, units='mm'))
  
    
  # plot updated fleet
  pdf(paste0(path.out.plot,"/updatedfleet.pdf"), width=7*16/9, height=7)
  #print(g)
  grid.draw(gg)
  dev.off()
}

get.demand <- function(path.output.ce, path.data.ce,
                       path.rcps=c('RCP 4.5'=sprintf('%s/rcp45/',
                                                     path.output.ce),
                                   'RCP 8.5'=sprintf('%s/rcp85/',
                                                     path.output.ce))) {
  
  df.demand <- data.frame(case=as.character(),
                          Region.Name=as.character(),
                          currYear=as.numeric(),
                          mean.load=as.numeric(),
                          peak.load=as.numeric(), stringsAsFactors = FALSE)
  
  for (rcp in path.rcps) {
    case.name <- names(which(path.rcps == rcp))
    for (y in seq(2020, 2050, by=5)){
      ok <- FALSE
      fname <- paste0(paste0(rcp, 'demandFullYrZonalCE', y,'.csv'))
      print(fname)

      if (file.exists(fname)) {
        df <- read.csv(file=fname, stringsAsFactors = FALSE, header = FALSE)
        df$V1 <- NULL
        df <- df %>% gather('hour', 'demand', -V2) %>% 
          mutate(hour=as.numeric(gsub('V', '', hour))-2) %>% arrange(V2, hour) %>%
          dplyr::rename(Region.Name=V2)
        
        df <- df %>% group_by(Region.Name) %>% 
          dplyr::summarise(mean.load=mean(demand), peak.load=max(demand)) %>%
          mutate(currYear=y, case=case.name) %>% 
          select(case, Region.Name, currYear, mean.load, peak.load)
        
        df.demand <- rbind(df.demand, df)
      }
    }
    
  }
  
  return(df.demand)
  
  #ssh_disconnect(ssh.con)
}

plot_renewable_cfs <- function(path.output.ce,
                          gamsSysDir='/Applications/GAMS24.7/sysdir/',
                          path.rcps=c('RCP 4.5'=sprintf('%s/rcp45/',
                                                        path.output.ce),
                                      'RCP 8.5'=sprintf('%s/rcp85/',
                                                        path.output.ce)),
                          yearIni=2020, yearEnd=2050,
                          path.out.plot=path.output.ce) {
  
  igdx(gamsSysDir)
  
  rcp <- path.rcps[2]
  
  df.out <- data.frame(year=as.numeric(),
                       z=as.character(),
                       techrenew=as.character(),
                       hour.of.day=as.numeric(),
                       type.CF=as.character(),
                       value=as.numeric())
  serc.zones <- read.csv(paste0(rcp,'zoneNamesToNumbers.csv'))
  
  for (y in seq(2020, 2050, by=5)) {
    x <- rgdx.param(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
                    symName = 'pCf') %>%
      mutate(ZoneNum=as.numeric(as.factor(z))) %>%
      left_join(serc.zones, by='ZoneNum') %>%
      select(g, techrenew, h, pCf, Zone) %>%
      mutate(Zone=as.character(Zone), techrenew=as.character(techrenew))
    
    x <- x %>% mutate(h = as.numeric(gsub('h', '', as.character(h)))-1) %>%
      mutate(time=ymd_h(sprintf('%4d-01-01 %2d', y, 0), tz='EST') + hours(h)) %>%
      mutate(hour.of.day=hour(time)) %>% group_by(Zone, techrenew, hour.of.day) %>%
      dplyr::summarize(CF=mean(pCf)) %>% mutate(year=y) %>%
      select(year, Zone, techrenew, hour.of.day, CF) %>% as.data.frame() 
    
    w <- x %>% group_by(year, Zone, techrenew) %>% dplyr::summarise(mean.CF=sum(CF)/24)
    
    x <- x %>% 
      left_join(w, by=c('year'='year' , 'Zone'='Zone' , 'techrenew'='techrenew')) %>%
      gather(key = "type.CF", value = "value", CF:mean.CF)
      
    
    df.out <- rbind(df.out, x)
  }
  
  
  ggplot(df.out) + geom_line(aes(x=hour.of.day, y=value, color=as.factor(year),
                                 linetype=type.CF)) + theme_bw() + 
    facet_grid(rows = vars(Zone), cols = vars(techrenew))
}


plot.donut.chart <- function(fileGenfleet, currYear=2050, GW=TRUE, 
                             labels.size=3, center.label=5, 
                             adjust.labels.1 = c('Natural Gas'='Natural\nGas',
                                                 'Combined Cycle'='Combined\nCycle'),
                             adjust.labels.2 = NULL) {
  
  df.fleet <- read.csv(file=fileGenfleet, 
                       stringsAsFactors = FALSE) %>%
    # remove plants retired by age
    filter(YearRetiredByAge > currYear | is.na(YearRetiredByAge)) %>%
    # remove plants originally marked to retire in IPM
    filter(Retirement.Year > currYear) %>%
    # remove plants that were retired by the CE model due to low CF
    filter(is.na(YearRetiredByCE)) %>%
    mutate(PlantType = as.character(PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
    mutate(Modeled.Fuels = ifelse(PlantType == 'Other', 'Other', Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(Modeled.Fuels == 'Other', 'NA', PlantType)) %>%
    mutate(Modeled.Fuels = ifelse(Modeled.Fuels %in% coal.types, 'Coal', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(Modeled.Fuels == 'Nuclear Fuel', 'Nuclear', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Natural Gas', Modeled.Fuels), 'Natural Gas', Modeled.Fuels)) %>%
    mutate(Modeled.Fuels = ifelse(grepl('Fuel Oil', Modeled.Fuels), 'Fuel Oil', Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 
                              'Hydro', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType == 'Combustion Turbine', 'CT', PlantType)) %>%
    mutate(Modeled.Fuels = ifelse(PlantType == 'Hydro', 'Hydro', Modeled.Fuels)) %>%
    mutate(PlantType = ifelse(grepl('Hydro|Wind|Solar', PlantType), 'NA', PlantType)) %>%
    mutate(cooling=ifelse(grepl('recirculating', Cooling.Tech),
                          'RC',
                          ifelse(grepl('once through', Cooling.Tech),
                                 'OT',
                                 ifelse(grepl('dry cooling', Cooling.Tech),
                                        'DC', 'NA')))) %>%
    mutate(cooling = ifelse(PlantType == 'Hydro', 'NA', cooling)) %>%
    mutate(cooling = ifelse(Modeled.Fuels == 'Other', 'NA', cooling)) %>%
    mutate(PlantType = ifelse(PlantType == 'NA', NA, PlantType)) %>%
    mutate(cooling = ifelse(cooling == 'NA', NA, cooling)) %>%
    select(Plant.Name, PlantType, Modeled.Fuels, Capacity..MW.,cooling)
  
  df.fleet <- df.fleet %>% group_by(Modeled.Fuels, PlantType, cooling) %>% dplyr::summarize(sum=sum(Capacity..MW.)) %>% as.data.frame()
  
  if (GW) {
    total <- sum(df.fleet$sum)/1000
    label.center <- paste0("SERC\n", formatC(total, big.mark=',', digits=0,  format = 'f'), ' GW')
    size.center <- 3
  } else {
    total <- sum(df.fleet$sum)
    label.center <- paste0("SERC\n", formatC(total, big.mark=',', digits=0,  format = 'f'), ' MW')
    size.center <- 2.5
  }
  
  size.labels <- 2.9
  
  lvl0 <- data.frame(name = label.center, value = total, level = 0, fill = NA) %>%
    mutate(xmin=0, xmax=1, ymin=0, ymax=value) %>%
    mutate(x.avg=0, y.avg=0, colour=FALSE, label=name, size=size.center)
  
  lvl1 <- df.fleet %>% group_by(Modeled.Fuels) %>% dplyr::summarize(sum=sum(sum)) %>%
    ungroup() %>%
    mutate(name = Modeled.Fuels, value = sum, level = 1, fill = name) %>%
    select(name, value, level, fill) %>%
    arrange(fill, name) %>%
    mutate(xmin=1.05, xmax=2, ymin=0, ymax=cumsum(value)) %>%
    mutate(ymin=lag(ymax, default=0),
           x.avg=(xmin+xmax)/2,
           y.avg=(ymin+ymax)/2,
           colour=TRUE) %>%
    mutate(label=ifelse(value<700 | is.na(name), '', name), size=size.labels)
  
  lvl2 <- df.fleet %>% group_by(Modeled.Fuels, PlantType) %>%
    dplyr::summarize(value=sum(sum)) %>%
    ungroup() %>%
    mutate(name=PlantType) %>%
    mutate(fill=Modeled.Fuels) %>%
    select(name, value, fill) %>% mutate(level = 2) %>%
    select(name, value, level, fill) %>%
    arrange(fill, name) %>%
    mutate(fill=ifelse(is.na(name),paste0(fill, '_NA'), fill)) %>%
    mutate(xmin=2.05, xmax=3, ymin=0, ymax=cumsum(value)) %>%
    mutate(ymin=lag(ymax, default=0),
           x.avg=(xmin+xmax)/2,
           y.avg=(ymin+ymax)/2,
           colour=ifelse(grepl('_NA', fill), FALSE, TRUE)) %>%
    mutate(label=ifelse(value<700 | is.na(name), '', name), size=size.labels)
  
  
  lvl3 <- df.fleet %>% group_by(Modeled.Fuels, PlantType, cooling) %>%
    dplyr::summarize(value=sum(sum)) %>%
    arrange(Modeled.Fuels, PlantType, cooling) %>%
    ungroup() %>%
    mutate(name=cooling) %>%
    mutate(fill=Modeled.Fuels) %>%
    select(name, value, fill) %>% mutate(level = 3) %>%
    select(name, value, level, fill) %>%
    mutate(name=ifelse(grepl("Hydro", fill) & !is.na(name), NA, name)) %>%
    mutate(fill=ifelse(is.na(name),paste0(fill, '_NA'), fill)) %>%
    mutate(xmin=3.05, xmax=4, ymin=0, ymax=cumsum(value)) %>%
    mutate(ymin=lag(ymax, default=0),
           x.avg=(xmin+xmax)/2,
           y.avg=(ymin+ymax)/2,
           colour=ifelse(grepl('_NA', fill), FALSE, TRUE)) %>%
    mutate(label=ifelse(value<600 | is.na(name), '', name), size=size.labels)
    
  new.pallete <- rev(col.pallete)
  names(new.pallete) <- plant.types
  missing.cases <- rep("#FFFFFF00", length(new.pallete))
  names(missing.cases) <- paste0(plant.types, '_NA')
  
  new.pallete.2 <- c(new.pallete, missing.cases)
  
  df.final <- bind_rows(lvl0, lvl1, lvl2, lvl3) %>%
    mutate(name = as.factor(name)) %>%
    arrange(fill, name) 
  
  # remove 2nd and 3rd layers of techs that have small capacity (< 600).
  # (this was created because of small coal in the end)
  df.final <- df.final %>% mutate(aux=ifelse(value < 600 & level == 1, TRUE, FALSE)) %>%
    group_by(fill) %>% arrange(fill, level) %>% mutate(aux=aux[1]) %>% ungroup() %>% #print(n=100)
    mutate(fill=ifelse(level > 1 & aux, paste0(fill, "_NA"), fill),
           colour=ifelse(grepl('_NA', fill), FALSE, colour)) %>%
    select(-aux)
  
  df.final <- df.final %>% mutate(level = as.factor(level))
  
  # create column for adjusting hjust and vjust of labels
  df.final <- df.final %>%
    mutate(hjust=0.5, vjust=0.5) %>%
    mutate(hjust=ifelse(name=='DC' & fill == 'Natural Gas', 1, hjust),
           vjust=ifelse(name=='DC' & fill == 'Natural Gas', 0, vjust)) %>%
    mutate(hjust=ifelse(name=='OT' & fill == 'Natural Gas', 0, hjust),
           vjust=ifelse(name=='OT' & fill == 'Natural Gas', 1, vjust)) %>%
    mutate(hjust=ifelse(name=='Other', 1, hjust),
           vjust=ifelse(name=='Other', 0, vjust)) %>%
    mutate(hjust=ifelse(name=='Solar', 0.2, hjust),
           vjust=ifelse(name=='Solar', 0.9, vjust))
  
  # adjust some labels
  if (!is.null(adjust.labels.1)) {
    for (lab in names(adjust.labels.1)) {
      df.final <- df.final %>% 
        mutate(label=ifelse(label == lab, adjust.labels.1[lab], label))
    }
  }
  
  if (!is.null(adjust.labels.2)) {
    for (lab in names(adjust.labels.2)) {
      df.final <- df.final %>% 
        mutate(label=ifelse(label == lab, adjust.labels.2[lab], label))
    }
  }
  
  df.final <- df.final %>% mutate(level= as.numeric(as.character(level)))
  
  g <- ggplot(data=df.final, aes(fill = fill)) +  
    geom_rect(aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, color=colour), 
              size = 0.1) +
    scale_fill_manual(values = new.pallete.2, na.translate = FALSE) +
    scale_color_manual(values = c('TRUE'='gray20', 'FALSE'='#FFFFFF00'), 
                       guide = F, na.translate = FALSE) +
    # label for center
    geom_text(data=df.final %>% filter(level == 0), 
              aes(x = x.avg, y = y.avg, label = label, hjust=hjust, 
                  vjust=vjust), size = rel(center.label)) +
    # label for layers
    geom_text(data=df.final %>% filter(level > 0),
              aes(x = x.avg, y = y.avg, label = label, hjust=hjust, 
                  vjust=vjust), size = rel(labels.size)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,4)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() + 
    guides(fill=FALSE, size=FALSE) +
    theme(panel.border = element_blank(),
          plot.margin = unit(c(-2, -2, -2, -2), "lines"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.key = element_blank()) +
    coord_polar(theta = "y", clip='off')
  
  #pdf("~/CMU/RIPS/paperCE/longdraft/images/sunburst.pdf", width=5, height=5)
  #print(g)
  #dev.off()
  #system('pdfcrop ~/CMU/RIPS/paperCE/longdraft/images/sunburst.pdf ~/CMU/RIPS/paperCE/longdraft/images/sunburst.pdf')
  
  return(g)
}


get.generation <- function(y, path.data) {

  # y <- 2050
  # path.data <- '/Volumes/RIPS/CE/results/base_case/'
  
  # get fleet
  df.fleet <- get.fleet(paste0(path.data, "/genFleetAfterCE", y,".csv"), y) %>%
    select(id, Modeled.Fuels) %>% mutate(Modeled.Fuels=tolower(Modeled.Fuels))
    
  path.gdx <- paste0(path.data, '/gdxOutYear', y,'.gdx')
  # path.gdx <- "/Volumes/RIPS/CE/results//rcp85//gdxOutYear2045.gdx"
  # path.gdx <- "/Users/kiko/Documents/CE/out/CE/gdxOutYear2025.gdx"

  # get ids from existing solar, wind, hydro, pumped hydro
  ids.solar <- rgdx.set(gdxName = paste0(path.gdx), symName = "solaregu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='solar')
  ids.wind <- rgdx.set(gdxName = paste0(path.gdx), symName = "windegu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='wind')
  ids.hydro <- rgdx.set(gdxName = paste0(path.gdx), symName = "hydroegu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='hydro')
  ids.pumped <- rgdx.set(gdxName = paste0(path.gdx), symName = "pumphydroegu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='pumped')
  
  #types.exist <- rbind(ids.solar, ids.hydro, ids.pumped, ids.wind)
  
  # **existing** renewables are aggregated within each zone and given a **dummy code**
  # I get the **dummy codes** here in order to map generation data of existing 
  # wind/solar to their types
  # OBS: after the investment decisions are made, these dummy codes are given to new plants added
  # to the fleet.
  types.exist.renew <- rbind(ids.solar, ids.wind)
  
  # get hours of each season
  hours.summer <- rgdx.set(gdxName = path.gdx, symName = "summerh") %>% mutate(season='summer')
  hours.winter <- rgdx.set(gdxName = path.gdx, symName = "winterh") %>% mutate(season='winter')
  hours.spring <- rgdx.set(gdxName = path.gdx, symName = "springh") %>% mutate(season='spring')
  hours.fall <- rgdx.set(gdxName = path.gdx, symName = "fallh") %>% mutate(season='fall')
  hours.special <- rgdx.set(gdxName = path.gdx, symName = "specialh") %>% mutate(season='special')
  
  df.hours.season <- rbind(hours.summer, hours.winter, hours.fall, hours.spring) %>%
    mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
    mutate(date=as.POSIXct(ymd(paste0(y,"0101")) + hours(h2))) %>%
    arrange(g, season, date) %>%
    group_by(g, season) %>% mutate(hour.in.season= row_number()) %>% ungroup()
  
  # for special hours put first peak demand than peak curtailment
  df.out <- get.hours.max.curtailment(path.data, y)
  
  hours.special <- hours.special %>% mutate(g=as.character(g), h=as.character(h))
  
  if (!is.null(df.out)) {
    df.out$order <- 2
    hours.special <- hours.special %>% left_join(df.out, by=c("g", "h")) %>%
      mutate(order=ifelse(is.na(order), 1, order))
  } else {
    hours.special$order <- 1
  }
  hours.special <- hours.special %>% group_by(g) %>%
    arrange(g, order, h) %>% mutate(hour.in.season= row_number()) %>% ungroup() %>%
    mutate(h2=as.numeric(gsub("h", "", h)),
           date=as.POSIXct(ymd(paste0(y,"0101")) + hours(h2))) %>%
    select(g, h, season, h2, date, hour.in.season)
    
  df.hours.season <- rbind(df.hours.season, hours.special)
      
  pDemand <- rgdx.param(gdxName = path.gdx, symName = "pDemand") %>% 
    group_by(g, h) %>% summarise(value = sum(pDemand)) %>%
    right_join(df.hours.season, by=c('g','h'))
    
  charge.ph <- rgdx(gdxName = path.gdx, requestList = list(name='vCharge', field='l'))
  if (!is.null(charge.ph)) {
    charge.ph <- convert.gdx.result.df(charge.ph) %>% 
      group_by(g, h) %>% summarise(value = sum(value))  %>%
      ungroup() %>% right_join(df.hours.season, by=c('g','h')) %>%
      mutate(value = ifelse(is.na(value), 0, value))
  }
  
  pDemand <- pDemand %>% left_join(charge.ph %>% select(g, h, value), 
                                      by=c('g','h')) %>%
    mutate(value=value.x+value.y) %>% 
    mutate(type="demand") %>%
    select(g, type, h, value, season, h2, date, hour.in.season) %>% ungroup()

  # get generation
  gexist <- rgdx(gdxName = path.gdx, 
                 requestList = list(name='vPegu', field='l', compress=FALSE,
                                    form = 'full'))
  gexist <- convert.gdx.result.df(gexist) 
  
  # first get wind/solar
  gexist.renew <-  gexist %>% 
    mutate(egu=as.character(egu)) %>%
    inner_join(types.exist.renew, by= 'egu') %>%
    group_by(g, type, h) %>% summarise(value = sum(value)) %>%
    ungroup() %>% right_join(df.hours.season, by=c('g','h'))

  # others (not renewables)
  gexist.others <- gexist %>% 
    mutate(egu=as.character(egu)) %>%
    anti_join(types.exist.renew, by= 'egu') %>%
    left_join(df.fleet, by= c('egu' = "id")) %>%
    dplyr::rename(type=Modeled.Fuels) %>%
    group_by(g, type, h) %>% summarise(value = sum(value)) %>%
    ungroup() %>% right_join(df.hours.season, by=c('g','h'))
  
  gexist <- rbind(gexist.renew, gexist.others)
  
  gtechcurt <- rgdx(gdxName = path.gdx, 
                    requestList = list(name='vPtechcurtailed', field='l', compress=FALSE,
                                       form = 'full'))
  gtechcurt  <- convert.gdx.result.df(gtechcurt) 

  if (!is.null(gtechcurt) && nrow(gtechcurt) > 0){
    gtechcurt <- gtechcurt %>% 
      group_by(g, techcurtailed, h) %>% summarise(value = sum(value)) %>%
      right_join(df.hours.season, by=c('g','h')) %>%
      ungroup() %>%
      mutate(techcurtailed=as.character(techcurtailed)) %>%
      mutate(techcurtailed=sapply(techcurtailed, FUN = {function(x){strsplit(x, "+", fixed=TRUE)[[1]][1]}})) %>%
      mutate(techcurtailed=simplify.thermal.names(techcurtailed)) %>%
      dplyr::rename(type=techcurtailed) %>%
      select(g, type, h, value, season, h2, date, hour.in.season)
  }
  
  gtechnotcurt <- rgdx(gdxName = path.gdx, 
                       requestList = list(name='vPtechnotcurtailed', field='l', 
                                          compress=FALSE, form = 'full'))
  gtechnotcurt  <- convert.gdx.result.df(gtechnotcurt)
  if (!is.null(gtechnotcurt) && nrow(gtechnotcurt) > 0){
    gtechnotcurt <- gtechnotcurt %>% 
      group_by(g, h) %>% summarise(value = sum(value)) %>%
      right_join(df.hours.season, by=c('g','h')) %>% 
      ungroup() %>%
      mutate(technotcurtailed=as.character(technotcurtailed)) %>%
      mutate(technotcurtailed=sapply(technotcurtailed, FUN = {function(x){strsplit(x, "+", fixed=TRUE)[[1]][1]}})) %>%
      mutate(technotcurtailed=simplify.thermal.names(technotcurtailed)) %>%
      dplyr::rename(type=technotcurtailed) %>%
      select(g, type, h, value, season, h2, date, hour.in.season)
  }
  
  gtechrenew <- rgdx(gdxName = path.gdx, 
                     requestList = list(name='vPtechrenew',field='l', 
                                        compress=FALSE, form = 'full'))
  gtechrenew  <- convert.gdx.result.df(gtechrenew)
  if (!is.null(gtechrenew) && nrow(gtechrenew) > 0){
    gtechrenew <- gtechrenew %>% 
      mutate(techrenew = as.character(techrenew)) %>%
      mutate(techrenew = ifelse(grepl('Solar PV\\+\\d\\d',techrenew, perl = TRUE), 
                                'solar', 'wind')) %>%
      group_by(g, techrenew, h) %>% summarise(value = sum(value)) %>%
      right_join(df.hours.season, by=c('g','h')) %>% ungroup() %>%
      dplyr::select(g, type=techrenew, h, value, season, h2, date, hour.in.season) %>%
      ungroup()
  }
  
  df.total <- gexist
  if (!is.null(gtechcurt) && nrow(gtechcurt) > 0) {
    df.total <- rbind(df.total, gtechcurt, stringsAsFactors = FALSE)
  }
  if (!is.null(gtechrenew) && nrow(gtechrenew) > 0) {
    df.total <- rbind(df.total, gtechrenew, stringsAsFactors = FALSE)
  }
  if (!is.null(gtechnotcurt) && nrow(gtechnotcurt) > 0) {
    df.total <- rbind(df.total, gtechnotcurt, stringsAsFactors = FALSE)
  }
  df.total <- rbind(df.total, pDemand, stringsAsFactors = FALSE)
  
  df.total <- df.total %>% group_by(g, h, type) %>% summarise(value = sum(value)) %>% ungroup() %>%
    right_join(df.hours.season, by=c('g','h')) %>% ungroup() 

  df.total$seasonweights <- 1
  
  seasons <- df.total %>% dplyr::filter(season != 'special') %>% select(season) %>% unlist() %>%
    unique()
  for (s in seasons) {
    w.season <- rgdx.scalar(gdxName = path.gdx, symName = sprintf("pWeight%s", s))
    df.total$seasonweights <- ifelse(df.total$season == s, w.season, df.total$seasonweights)
  }
  
  return(df.total) 
  
}

plot.generation <- function(path.rcps, name.file, year=2050,
                           width=7, height=7*9/16){
  
  list.gen <- list()
  for (rcp in names(path.rcps)) {
    list.gen[[rcp]] <- get.generation(year, path.rcps[rcp])
  }
  df.gen2 <- ldply(list.gen, data.frame) %>% group_by(.id, season, hour.in.season, type) %>%
    summarise(value = mean(value)) %>% mutate(type=tolower(type))
  
  pDemand.2 <- df.gen2 %>% filter(type == "demand")
  df.gen2 <- df.gen2 %>% filter(type != "demand")
  
  df.gen2 <- df.gen2 %>% ungroup() %>%
    mutate(type = factor(type, levels = rev(c("other", "hydro","nuclear","coal", "natural gas", "solar", "wind")))) %>%
    mutate(season = factor(season, levels = c("winter", "spring", "summer", "fall", "special")))
  
  pDemand.2 <- pDemand.2 %>% ungroup() %>% mutate(season = factor(season, levels = c("winter", "spring", "summer", "fall", "special")))
  
  g <- ggplot(df.gen2) +
    geom_area(aes(x=hour.in.season, y=value, fill=type)) +
    geom_line(data = pDemand.2, aes(x=hour.in.season, y=value, linetype='demand'),
              colour='black') +
    guides(colour=FALSE) + theme_minimal() +
    ylab('GWh') + xlab('hour in simulation period') + 
    scale_fill_brewer(type='div', palette = 'Spectral') +
    scale_linetype_manual(values=c('dashed')) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0,'lines'),
          panel.border = element_rect(fill = NA),
          legend.position = 'top') +
    guides(fill=guide_legend(title=NULL, nrow = 1, reverse = TRUE),
           linetype=guide_legend(title=NULL, nrow = 1)) +
    facet_grid(cols = vars(season), rows=vars(.id), scales = 'free_x')
  
  pdf(name.file, width=width, height=height)
  print(g)
  dev.off()
  
  return(df.gen2)
}

simplify.thermal.names <- function(x) {
  
  y <- ifelse(grepl("Combined Cycle", x), "natural gas", x)
  y <- ifelse(grepl("Coal", y), "coal", y)
  y <- ifelse(grepl("Nuclear", y), "nuclear", y)
  
  return(y)
}


convert.gdx.result.df <- function(list.gdx) {
  
  require(reshape2)
  
  if (list.gdx$form == 'full') {
    df.gdx <- reshape2::melt(list.gdx$val)
    
  } else { 
    df.gdx <- as.data.frame(list.gdx$val)
    
    if (nrow(df.gdx) > 0) {
      names(df.gdx) <- c(list.gdx$domains, 'value')
      
      for (i in 1:length(list.gdx$domains)) {
        dom.values <- list.gdx$uels[[i]]
        df.gdx[[i]] <- factor(df.gdx[[i]], levels = 1:length(dom.values), labels=c(dom.values))
      }
    } else {
      df.gdx <- NULL
    }
  }
  
  
  return (df.gdx)
  
}

define.season <- function(dt) {
  
  s <- ifelse(month(dt) %in% c(3,4,5), 'spring',
              ifelse(month(dt) %in% c(6,7,8), 'summer',
                     ifelse(month(dt) %in% c(9,10,11), 'fall', 'winter')))
  
  return(s)
  
} 


check.pumpedhydro <- function(y, path.gdx) {
  
  #path.gdx <- "/Volumes/Seagate_RIPS/CE/results//rcp85//gdxOutYear2045.gdx"
  path.gdx <- "/Users/kiko/Documents/CE/out/CE/gdxOutYear2025.gdx"

  x <- gdxInfo(gdxName = paste0(path.gdx), dump=FALSE, returnList=TRUE, returnDF=FALSE)
  
  ids.pumped <- rgdx.set(gdxName = path.gdx, symName = "pumphydroegu") %>% 
    mutate(egu=as.character(egu)) %>% unlist()  
  
  # get hours of each season
  hours.summer <- rgdx.set(gdxName = path.gdx, symName = "summerh") %>% mutate(season='summer')
  hours.winter <- rgdx.set(gdxName = path.gdx, symName = "winterh") %>% mutate(season='winter')
  hours.spring <- rgdx.set(gdxName = path.gdx, symName = "springh") %>% mutate(season='spring')
  hours.fall <- rgdx.set(gdxName = path.gdx, symName = "fallh") %>% mutate(season='fall')
  df.hours.season <- rbind(hours.summer, hours.winter, hours.fall, hours.spring) %>%
    mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
    mutate(date=as.POSIXct(ymd("20200101") + hours(h2)))
  
  gexist <- rgdx(gdxName = path.gdx, requestList = list(name='vPegu', field='l'))
  gexist <- convert.gdx.result.df(gexist) %>% filter(egu %in% ids.pumped) %>%
    mutate(egu=as.character(egu)) %>% 
    group_by(g, h) %>% summarise(value = sum(value)) %>%
    rename(value.exist = value)
  
  charge.ph <- rgdx(gdxName = path.gdx, requestList = list(name='vCharge', field='l'))
  if (!is.null(charge.ph)) {
    charge.ph <- convert.gdx.result.df(charge.ph) %>% 
      group_by(g, h) %>% summarise(value = sum(value)) %>%
      rename(value.charge = value)
  }
  
  soc.ph <- rgdx(gdxName = path.gdx, requestList = list(name='vSoc', field='l'))
  if (!is.null(soc.ph)) {
    soc.ph <- convert.gdx.result.df(soc.ph) %>% 
      group_by(g, h) %>% summarise(value = sum(value)) %>%
      rename(value.soc = value)
  }
  
  
  df.ph.total <- df.hours.season %>% left_join(soc.ph, by=c('g', 'h')) %>%
    left_join(gexist, by=c('g', 'h')) %>%
    left_join(charge.ph, by=c('g', 'h')) %>%
    mutate_at(.vars=vars(value.soc, value.exist, value.charge), 
              .funs={function(x){ifelse(is.na(x), 0, x)}}) %>%
    gather(key='type', value='value', value.soc, value.exist, value.charge)
  
  ggplot(df.ph.total) + geom_line(aes(x=date, y=value, linetype=type), colour='red') + 
    guides(colour=FALSE) + theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0,'lines'),
          panel.border = element_rect(fill = NA)) + 
    facet_grid(cols = vars(season), scales = 'free_x')
  
}


compute.annual.emmission <- function(y, path.gdx) {
  
  y <- 2020
  
  path.gdx <- "/Users/kiko/Documents/CE/out/CE/gdxOutYear2020.gdx"

  ids.solar <- rgdx.set(gdxName = paste0(path.gdx), symName = "solaregu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='solar')
  ids.wind <- rgdx.set(gdxName = paste0(path.gdx), symName = "windegu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='wind')
  ids.hydro <- rgdx.set(gdxName = paste0(path.gdx), symName = "hydroegu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='hydro')
  ids.pumped <- rgdx.set(gdxName = paste0(path.gdx), symName = "pumphydroegu") %>% 
    mutate(egu=as.character(egu)) %>% mutate(type='pumped')
  
  nrow(ids.hydro) + nrow(ids.pumped) + nrow(ids.wind)+ nrow(ids.solar)
  
  x <- rgdx(gdxName = path.gdx, requestList = list(name="vPegu", 
                                                   form='full',field=c('l')))
  
  pegu <- convert.gdx.result.df(x) %>% mutate(egu=as.character(egu))
  
  emrate.egu <- rgdx.param(gdxName = path.gdx, 'pCO2emrate', squeeze=FALSE) %>%
    transmute(egu=as.character(egu), pCO2emrate=pCO2emrate)
  
  pegu <- pegu %>% left_join(emrate.egu, by= 'egu') %>% 
    mutate(pCO2emrate=ifelse(is.na(pCO2emrate), 0, pCO2emrate))
  
  #pegu %>% filter(is.na(pCO2emrate)) %>% select(egu) %>% unique() %>% nrow()
  
  emrate.tech <- rgdx.param(gdxName = path.gdx, 'pCO2emratetechcurt', squeeze=FALSE) %>%
    transmute(tech=as.character(techcurtailed), pCO2emrate=pCO2emratetechcurt) %>% group_by(tech) %>%
    summarise(pCO2emrate=mean(pCO2emrate))
  
  x <- rgdx(gdxName = path.gdx, requestList = list(name="vPtechcurtailed", 
                                                   form='sparse',
                                                   field=c('l')))
  ptechcurt <- convert.gdx.result.df(x)
  if(!is.null(ptechcurt)) {
    ptechcurt <-  ptechcurt %>%
      mutate(genid=paste0(c,'_',techcurtailed), tech=techcurtailed) %>%
      left_join(emrate.tech, by='tech') %>%
      select(genid, h, gen, pCO2emrate) %>%
      mutate(pCO2emrate=ifelse(is.na(pCO2emrate), 0, pCO2emrate))
  }
  
  x <- rgdx(gdxName = path.gdx, requestList = list(name="vPtechrenew", 
                                                   form='full',
                                                   field=c('l')))
  ptechrenew <- convert.gdx.result.df(x)
  if(!is.null(ptechrenew)) {
    ptechrenew <- ptechrenew %>%
      mutate(genid=paste0(z,'_',techrenew), tech=techrenew) %>%
      left_join(emrate.tech, by='tech') %>%
      select(genid, h, gen, pCO2emrate) %>%
      mutate(pCO2emrate=ifelse(is.na(pCO2emrate), 0, pCO2emrate))
  }
  
  x <- rgdx(gdxName = path.gdx, requestList = list(name="vPtechnotcurtailed", 
                                                   form='full',field=c('l')))
  ptechnotcurt <- convert.gdx.result.df(x)
  if(!is.null(ptechnotcurt)) {
    ptechnotcurt <- ptechnotcurt %>%
      mutate(genid=paste0(z,'_','Nuclear'), tech='Nuclear') %>%
      left_join(emrate.tech, by='tech') %>%
      select(genid, h, gen, pCO2emrate) %>%
      mutate(pCO2emrate=ifelse(is.na(pCO2emrate), 0, pCO2emrate))
  }
  
  # row bind all dfs into one
  
  df.gen.ce <- pegu
  if (!is.null(ptechcurt)) {
    df.gen.ce <- rbind(df.gen.ce, ptechcurt)
  }
  if (!is.null(ptechnotcurt)) {
    df.gen.ce <- rbind(df.gen.ce, ptechnotcurt)
  }
  if (!is.null(ptechrenew)) {
    df.gen.ce <- rbind(df.gen.ce, ptechrenew)
  }

  df.gen.ce <- df.gen.ce %>% mutate(CO2emission=value*pCO2emrate)
  
  hours <- data.frame(h=as.character(), season=as.character())
  season.weigths <- data.frame(season=as.character(), weigth=as.numeric())
  for(s in c('summer', 'winter', 'fall', 'spring', 'special')) {
    zz <- rgdx.set(gdxName = path.gdx, sprintf('%sh', s))
    
    hours <- rbind(hours, data.frame(h=as.character(zz$h), 
                                     season=s))
    if (s != 'special') {
      aux <- rgdx.scalar(gdxName = path.gdx,sprintf('pWeight%s', s))
      season.weigths <- rbind(season.weigths,
                              data.frame(season=s,
                                         weights=aux))
    }
  }
  hours <- hours %>% mutate(h=as.character(h),
                            season=as.character(season))
  season.weigths <- season.weigths %>% mutate(season=as.character(season))
  
  df.gen.ce <- df.gen.ce %>% left_join(hours, by='h')
  
  ce.generation <- df.gen.ce %>% 
    # group by seasons and genid, sum season gen by unit, add columns with weights to normallize to total seasonal value
    group_by(g, season) %>% summarize(gen=sum(value, na.rm=TRUE),
                                      emission=sum(CO2emission, na.rm=TRUE)) %>%
    left_join(season.weigths, by='season') %>% filter(!is.na(season)) %>%
    mutate(weights=ifelse(is.na(weights), 1, weights)) %>%
    mutate(gen=gen*weights, emission=emission*weights) %>% ungroup()
  
  ce.generation %>% group_by(g) %>% summarize(annual.emission=sum(emission)) %>% select(annual.emission) %>% unlist() %>% mean()
  
}


get.hours.max.curtailment <- function(path.data, y) {
  
  df.out <- NULL
  if (file.exists(paste0(path.data, '/curtailmentsHourlyCombined', y, '.csv'))) {
    df <- read.csv(paste0(path.data, '/curtailmentsHourlyCombined', y, '.csv'), header = FALSE, stringsAsFactors = FALSE) %>%
      pivot_longer(-c(V1), names_to = "hour", values_to = "value") %>%
      mutate(hour=as.numeric(gsub("V", "", hour))-1) %>%
      group_by(V1) %>% dplyr::arrange(value) %>%
      summarise(hour=hour[1], value=value[1]) %>% ungroup() %>%
      mutate(hour.day=ifelse(hour %% 24 == 0, 24, hour %% 24))
    
    df.out <- NULL  
    for(i in unique(df$V1)) {
      x <- df %>% dplyr::filter(V1 == i)
      initial.hour.day <- x$hour - x$hour.day + 1
      
      hours.days <- paste0("h", seq(from=initial.hour.day, length.out = 24))
      
      if (is.null(df.out)) {
        df.out <- data.frame(g=i, h=hours.days)
      } else {
        df.out <- rbind(df.out, data.frame(g=i, h=hours.days))
      }
    }
    
    df.out <- df.out %>% mutate(g=as.character(g), h=as.character(h))
  }
  
  return(df.out)
}

# # check supply demand equations
# # meetdemand <- rgdx(gdxName = path.gdx, requestList = list(name='meetdemand', field='l'))
# rgdx.param(gdxName = path.gdx, symName = 'pLinesources')
# 
# 
# x <- gdxInfo(gdxName = paste0(path.gdx), dump=FALSE, returnList=TRUE, returnDF=FALSE)
# 
# # get capacity from all existing generators
# pCapacExist <- rgdx.param(gdxName = paste0(path.gdx), symName = "pCapac")
# 
# # get max energy from existing solar generators
# pMaxgensolar <- rgdx.param(gdxName = path.gdx, symName = "pMaxgensolar") %>%
#   group_by(g, h) %>% summarise(value = sum(pMaxgensolar))
# pCapacSolar <- pCapacExist %>% filter(egu %in% ids.solar) %>% group_by(h) %>% summarise(value=sum(pCapac))
# 
# pMaxgensolar <- pMaxgensolar %>% mutate(capac = mean(pCapacSolar$value),
#                                         cf = value/capac) %>%
#   mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
#   mutate(hour.of.day=hour(ymd("20200101") + hours(h2)))
# 
# # get max energy from existing wind generators
# pMaxgenwind <- rgdx.param(gdxName = paste0(path.gdx), symName = "pMaxgenwind") %>%
#   group_by(g, h) %>% summarise(value = sum(pMaxgenwind))
# 
# ids.non.thermal <- c(ids.wind, ids.solar, ids.hydro, ids.pumped)
# 
# 
# # remove existing solar, wind, hydro, pumped hydro from data set
# capac.thermal <- rgdx.param(gdxName = paste0(path.gdx), symName = "pCapac") %>%
#   mutate(egu=as.character(egu)) %>% filter(!(egu %in% ids.non.thermal)) %>% 
#   group_by(g, h) %>% summarise(value = sum(pCapac))
# 
# ggplot(capac.thermal) + geom_line(aes(x=h, y=value, colour=g, group=g))
# 
# 
# x <- pMaxgensolar %>% group_by(hour.of.day) %>% 
#   summarise(cf.mean=mean(cf), lb=quantile(cf, 0.025), ub=quantile(cf, 0.975),
#             sd=sd(cf))
# 
# ggplot(x) + geom_line(aes(x=hour.of.day, y=cf.mean)) + geom_errorbar(aes(x=hour.of.day, ymin=lb, ymax=ub))
# ggplot(x) + geom_line(aes(x=hour.of.day, y=cf.mean)) + geom_errorbar(aes(x=hour.of.day, ymin=cf.mean - sd*1.96, ymax=cf.mean + sd*1.96))

#
# lifetimes <- read.csv(paste0(path.data.ce, '/NewPlantData/',
#                              'LifetimeValuesExistingPlants4Aug2016.csv'))
# 
# df.new.techs <- read.csv(paste0(path.rcps[1],sprintf('newTechsCE%4d.csv', y)))
# df.new.techs <- df.new.techs %>% 
#   transmute(Type=TechnologyType, Fuel=FuelType, 
#             Cooling=as.character(Cooling.Tech),
#             Capacity=as.numeric(Capacity.MW.), 
#             CAPEX=as.numeric(CAPEX.2012..MW.)/1e3,
#             FixedOM=as.numeric(FOM.2012..MW.yr.)/1e3, 
#             VarOM=as.numeric(VOM.2012..MWh.)) %>%
#   arrange(Type, Cooling) %>%
#   mutate(Cooling=ifelse(Type == 'Nuclear', NA, Cooling)) %>%
#   mutate(cool2=ifelse(is.na(Cooling), '',
#                       ifelse(Cooling=='dry cooling', 'DC',
#                              ifelse(Cooling=='recirculating', 'RC',
#                                     ifelse(Cooling == 'once-through', 
#                                            'OT', ''))))) %>%
#   mutate(type.2=ifelse(is.na(Cooling),
#                        gsub(' ', '.', Type),
#                        gsub(' ', '.', paste0(Type, '.', cool2)))) 
# 
# 
# df <- get.new.additions(2050, paste0(path.output.ce, '/rcp45/')) %>%
#   mutate(type2=trimws(gsub('\\.',' ', gsub('DC|RC', '', type, perl = TRUE),
#                            perl = TRUE))) %>% 
#   left_join(lifetimes, by=c('type2' = 'PlantType'))
# 
# df.2 <- df %>% left_join(df.new.techs, by=c('type' = 'type.2')) %>%
#   select(year, Zone, type, value, Fuel, Capacity, CAPEX, FixedOM, VarOM, Lifetime.yrs.) %>%
#   group_by(year, type, Fuel) %>%
#   summarise(value=sum(value),
#             Capacity=mean(Capacity),
#             CAPEX=mean(CAPEX),
#             FixedOM=mean(FixedOM),
#             VarOM=mean(VarOM),
#             n=mean(Lifetime.yrs.)) %>%
#   mutate(Total.MW = value*Capacity) %>%
#   mutate(annual.fix.cost=Total.MW*(CAPEX*(rate*(1+rate)^n/((1+rate)^n-1))+FixedOM))
# 
# 
# 


# map.cell2zone <- rgdx.param(gdxName = paste0(rcp, 'gdxOutYear', y, 
#                                              '.gdx'), 'pCellzones')
# map.cell2zone$c <- as.character(map.cell2zone$c)
# 
# tech.not.curtailed <- rgdx(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
#                            requestList = list(name="vNnotcurtailed", 
#                                               form='full', 
#                                               field=c('l')))
# 
# df.tech.not.curtailed <- data.frame(zone=tech.not.curtailed$uels$z,
#                                     tech.not.curtailed$val, 
#                                     row.names = NULL) %>%
#   mutate(ZoneNum=as.numeric(as.factor(zone))) %>%
#   select(-zone) %>%
#   gather(key='type', value='value', Nuclear) %>%
#   left_join(serc.zones, by='ZoneNum') %>%
#   arrange(ZoneNum) %>% mutate(year=y) %>%
#   select(year, Zone, type, value) %>% as.data.frame()
# 
# tech.curtailed <- rgdx(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
#                        requestList = list(name="vNcurtailed", 
#                                           form='full',field=c('l')))
# 
# df.tech.curtailed <- data.frame(cells=rownames(tech.curtailed$val),
#                                 tech.curtailed$val, row.names = NULL,
#                                 stringsAsFactors = FALSE) %>%
#   left_join(map.cell2zone, by=c("cells" = "c")) %>%
#   mutate(ZoneNum=pCellzones) %>% select(-pCellzones) %>%
#   left_join(serc.zones, by='ZoneNum') %>%
#   gather(key='type', value='value', -c(cells, ZoneNum, Zone)) %>% 
#   arrange(ZoneNum) %>% mutate(year=y)
# 
# df.tech.curtailed.zone <- df.tech.curtailed %>% group_by(Zone, type) %>%
#   dplyr::summarise(value=sum(value, na.rm=TRUE)) %>% mutate(year=y) %>%
#   select(year, Zone, type, value) %>% as.data.frame()
# 
# tech.renew <- rgdx(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
#                    requestList = list(name="vNrenew", 
#                                       form='full',field=c('l')))
# 
# df.tech.renew <- data.frame(zone=tech.renew$uels$z,
#                             tech.renew$val, row.names = NULL) %>%
#   mutate(ZoneNum=as.numeric(as.factor(zone))) %>%
#   select(-zone) %>%
#   gather(key='type', value='value', -ZoneNum) %>%
#   left_join(serc.zones, by='ZoneNum') %>%
#   arrange(ZoneNum) %>% mutate(type=gsub('\\.\\d\\d', '', type)) %>% 
#   group_by(Zone, type) %>%
#   dplyr::summarise(value=sum(value, na.rm=TRUE)) %>% mutate(year=y) %>%
#   select(year, Zone, type, value) %>% ungroup() %>% as.data.frame()
