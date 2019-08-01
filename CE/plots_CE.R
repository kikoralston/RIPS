# Jannuary 2019
# Set of functions that read CE model output and create plots for the analysis


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
  
  idx.cooling <- grepl('\\.OT|\\.DC|\\.RC', x, perl=TRUE)
  
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
  
  out <- gsub('\\.', ' ', out)
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
    mutate(Type = ifelse(Type == 'Nuclear' | is.na(Cooling), as.character(Type), 
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

plot.map.serc.with.plants <- function(path.data.ce, path.plot, 
                                      gen.fleet=NULL, show.points=TRUE, 
                                      width=7){
  
  g <- get.map.serc(path.data.ce)

  if (is.null(gen.fleet)) {
    gen.fleet <- read.csv(file=paste0(path.output.ce, '/base_case/genFleetInitial.csv'))  
  }
  
  gen.fleet.2 <- gen.fleet %>% 
    select(Plant.Name, UniqueID_Final, PlantType, Region.Name, State.Name, County, 
           Capacity..MW., On.Line.Year, Retirement.Year,
           YearRetiredByAge, Modeled.Fuels, Longitude, Latitude) %>%
    mutate(PlantType = as.character(PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
    mutate(PlantType = ifelse(Modeled.Fuels %in% coal.types, 'Coal', PlantType)) %>%
    mutate(PlantType = ifelse(grepl('Natural Gas', Modeled.Fuels), 
                              'Natural Gas', PlantType)) %>%
    mutate(PlantType = ifelse(grepl('Fuel Oil', Modeled.Fuels), 
                              'Fuel Oil', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 
                              'Hydro', PlantType)) %>%
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
  
  # assume ratio 16:9
  pdf(paste0(path.plot, "/sercmap_fleet.pdf"), width=width, height=width*9/16)
  print(g_plants)
  dev.off()
  system(sprintf('pdfcrop %s %s', paste0(path.plot, "/sercmap_fleet.pdf"), 
                 paste0(path.plot, "/sercmap_fleet.pdf")))
  
  return(g_plants)
}


get.map.serc <- function(path.data.ce) {
  
  bounds.longitude <- c(-91, -76)
  bounds.latitude <- c(29, 39)
  
  serc.zones <- c('S_SOU', 'S_C_TVA', 'S_VACA', 'S_C_KY')
  
  s1 <- readOGR(paste0(path.data.ce, 'IPMRegionsShapeFile/'), 
                'IPM_Regions_20121023_US')
  
  zones.map <- ifelse(s1@data$IPM_Region %in% serc.zones, 
                      s1@data$IPM_Region, NA)
  
  s1.union.serc <- unionSpatialPolygons(s1, zones.map)
  
  df.ipm.counties <- s1@data %>% 
    mutate(FIPS=as.numeric(as.character(FIPS)))
  
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
  
  labels_regions <- data.frame(x=as.numeric(),
                               y=as.numeric(),
                               label=as.character())
  
  for (p in s1.union.serc@polygons) {
    x <- p@labpt[1]
    y <- p@labpt[2]
    label <- levels(df.ipm.counties$IPM_Region)[as.integer(p@ID)]
    
    labels_regions <- rbind(labels_regions,
                            data.frame(x=x,
                                       y=y,
                                       label=label))
  }
  
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
    mutate(PlantType = as.character(PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% otherTypes, 'Other', PlantType)) %>%
    mutate(PlantType = ifelse(Modeled.Fuels %in% coal.types, 'Coal', PlantType)) %>%
    mutate(PlantType = ifelse(grepl('Natural Gas', Modeled.Fuels), 
                              'Natural Gas', PlantType)) %>%
    mutate(PlantType = ifelse(grepl('Fuel Oil', Modeled.Fuels), 
                              'Fuel Oil', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType %in% c('Hydro', 'Pumped Storage'), 
                              'Hydro', PlantType)) %>%
    mutate(PlantType = ifelse(PlantType == 'Solar PV', 'Solar', PlantType)) %>%
    select(-Modeled.Fuels)
  
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
    mutate(Region.Name == as.character(Region.Name))
  
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

get.new.additions <- function(y, path.gdx) {
  
  serc.zones <- read.csv(paste0(path.gdx,'zoneNamesToNumbers.csv'))
  
  # map cell to zone
  map.cell2zone <- rgdx.param(gdxName = paste0(path.gdx, 'gdxOutYear', y, 
                                               '.gdx'), 'pCellzones')
  map.cell2zone$c <- as.character(map.cell2zone$c)
  
  tech.not.curtailed <- rgdx(gdxName = paste0(path.gdx, '/gdxOutYear', y, '.gdx'), 
                             requestList = list(name="vNnotcurtailed", 
                                                form='full', 
                                                field=c('l')))
  
  df.tech.not.curtailed <- data.frame(zone=tech.not.curtailed$uels$z,
                                      tech.not.curtailed$val, 
                                      row.names = NULL) %>%
    mutate(ZoneNum=as.numeric(as.factor(zone))) %>%
    select(-zone) %>%
    gather(key='type', value='value', Nuclear) %>%
    left_join(serc.zones, by='ZoneNum') %>%
    arrange(ZoneNum) %>% mutate(year=y) %>%
    select(year, Zone, type, value) %>% as.data.frame()
  
  tech.curtailed <- rgdx(gdxName = paste0(path.gdx, '/gdxOutYear', y, '.gdx'), 
                         requestList = list(name="vNcurtailed", 
                                            form='full',field=c('l')))
  
  df.tech.curtailed <- data.frame(cells=rownames(tech.curtailed$val),
                                  tech.curtailed$val, row.names = NULL,
                                  stringsAsFactors = FALSE) %>%
    left_join(map.cell2zone, by=c("cells" = "c")) %>%
    mutate(ZoneNum=pCellzones) %>% select(-pCellzones) %>%
    left_join(serc.zones, by='ZoneNum') %>%
    gather(key='type', value='value', -c(cells, ZoneNum, Zone)) %>% 
    arrange(ZoneNum) %>% mutate(year=y)
  
  df.tech.curtailed.zone <- df.tech.curtailed %>% group_by(Zone, type) %>%
    dplyr::summarise(value=sum(value, na.rm=TRUE)) %>% mutate(year=y) %>%
    select(year, Zone, type, value) %>% as.data.frame()
  
  tech.renew <- rgdx(gdxName = paste0(path.gdx, '/gdxOutYear', y, '.gdx'), 
                     requestList = list(name="vNrenew", 
                                        form='full',field=c('l')))
  
  df.tech.renew <- data.frame(zone=tech.renew$uels$z,
                              tech.renew$val, row.names = NULL) %>%
    mutate(ZoneNum=as.numeric(as.factor(zone))) %>%
    select(-zone) %>%
    gather(key='type', value='value', -ZoneNum) %>%
    left_join(serc.zones, by='ZoneNum') %>%
    arrange(ZoneNum) %>% mutate(type=gsub('\\.\\d\\d', '', type)) %>% 
    group_by(Zone, type) %>%
    dplyr::summarise(value=sum(value, na.rm=TRUE)) %>% mutate(year=y) %>%
    select(year, Zone, type, value) %>% ungroup() %>% as.data.frame()
  
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
                               path.out.plot=path.output.ce) {
  
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
    serc.zones <- read.csv(paste0(rcp,'zoneNamesToNumbers.csv'))
    case.name <- names(which(path.rcps == rcp))
    
    for (y in seq(yearIni, yearEnd, by=5)) {
      
      path.gdx <- paste0(rcp, '/gdxOutYear', y, '.gdx')
      
      df <- get.new.additions(y, path.gdx)

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
    select(-off)
  
  g <- ggplot() + geom_col(data=df.complete.decision.3, 
                           aes(x=year, y=total.cap/1e3, fill=Type), 
                           width=width-0.1)
  
  g <- g + scale_fill_manual(limits = rev(plant.types), 
                             values = col.pallete) +
    theme_bw() + ylab('New Capacity (GW)') +
    guides(fill=guide_legend(title=NULL)) +
    theme_bw() +
    theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "lines"),
          axis.title.x = element_blank(),
          panel.spacing = unit(0.1, "lines")) +
    facet_grid(Zone ~ ., scales = 'free_y')

  g_param <- ggplot_build(g)
  x.range <- g_param$layout$panel_params[[1]]$x.range
  year.values <- unique(df.complete.decision.3$year)
  
  df.top.axis <- df.complete.decision.3 %>%
    mutate(xpos=year) %>%
    select(case, xpos) %>%
    distinct() %>%
    mutate(case=ifelse(case=='Reference', 'Ref.', case))
    
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
#          panel.border = element_rect(fill = NA)
    )
  
  lay <- rbind(c(rep(1,11), rep(2,163), rep(3, 26)),
               c(rep(4,200)))
  # 
  gg <- arrangeGrob(rectGrob(gp=gpar(fill="white", lty='blank')),
                    g.top.axis,
                    #rectGrob(gp=gpar(fill="white", lty='blank')),
                    rectGrob(gp=gpar(fill="white", lty='blank')),
                    g, layout_matrix=lay,
                    heights=c(2,20),
                    padding=unit(0, units='mm'))
  
  # plot new additions by plant type
  pdf(paste0(path.out.plot,"/newadditions.pdf"), width=7*16/9, height=7)
  #print(g)
  grid.draw(gg)
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


plot.donut.chart <- function(fileGenfleet) {
  
  df.fleet <- read.csv(file=fileGenfleet, 
                       stringsAsFactors = FALSE) %>%
    # remove plants retired
    filter(YearRetiredByAge >= 2050 | is.na(YearRetiredByAge)) %>%
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
  
  total <- sum(df.fleet$sum)
  
  lvl0 <- data.frame(name = paste0("SERC\n", formatC(total, big.mark=',', 
                                                     digits=0, 
                                                     format = 'f'),
                                   ' MW'), 
                     value = total, level = 0, fill = NA) %>%
    mutate(xmin=0, xmax=1, ymin=0, ymax=value) %>%
    mutate(x.avg=0, y.avg=0, colour=FALSE, label=name)
  
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
    mutate(label=ifelse(value<700 | is.na(name), '', name))
  
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
    mutate(label=ifelse(value<700 | is.na(name), '', name))
  
  
  lvl3 <- df.fleet %>% group_by(Modeled.Fuels, PlantType, cooling) %>%
    dplyr::summarize(value=sum(sum)) %>%
    arrange(Modeled.Fuels, PlantType, cooling) %>%
    ungroup() %>%
    mutate(name=cooling) %>%
    mutate(fill=Modeled.Fuels) %>%
    select(name, value, fill) %>% mutate(level = 3) %>%
    select(name, value, level, fill) %>%
    mutate(fill=ifelse(is.na(name),paste0(fill, '_NA'), fill)) %>%
    mutate(xmin=3.05, xmax=4, ymin=0, ymax=cumsum(value)) %>%
    mutate(ymin=lag(ymax, default=0),
           x.avg=(xmin+xmax)/2,
           y.avg=(ymin+ymax)/2,
           colour=ifelse(grepl('_NA', fill), FALSE, TRUE)) %>%
    mutate(label=ifelse(value<600 | is.na(name), '', name))
  
  #col.pallete <- brewer.pal(9, 'YlOrBr')
  
  new.pallete <- rev(col.pallete)
  names(new.pallete) <- plant.types
  missing.cases <- rep("#FFFFFF00", length(new.pallete))
  names(missing.cases) <- paste0(plant.types, '_NA')
  
  new.pallete.2 <- c(new.pallete, missing.cases)
  
  df.final <- bind_rows(lvl0, lvl1, lvl2, lvl3) %>%
    mutate(name = as.factor(name)) %>%
    arrange(fill, name) %>%
    mutate(level = as.factor(level))
  
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
  df.final <- df.final %>% 
    mutate(label=ifelse(label == 'Natural Gas', 'Natural\nGas', label)) %>%
    mutate(label=ifelse(label == 'Combined Cycle', 'Combined\nCycle', label))
  
  g <- ggplot(data=df.final, aes(fill = fill)) +  
    geom_rect(aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, color=colour), 
              size = 0.1) +
    scale_fill_manual(values = new.pallete.2, na.translate = FALSE) +
    scale_color_manual(values = c('TRUE'='gray20', 'FALSE'='#FFFFFF00'), 
                       guide = F, na.translate = FALSE) +
    geom_text(aes(x = x.avg, y = y.avg, label = label, hjust=hjust, 
                  vjust=vjust), size = rel(3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,4)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() + 
    guides(fill=FALSE) +
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


levelized.cost.energy <- function(y, path.output.ce,
                                  gamsSysDir='/Applications/GAMS24.7/sysdir/',
                                  path.out.plot=path.output.ce) {
  
  rate <- 0.07
  igdx(gamsSysDir, silent=TRUE)
  
  lifetimes <- rgdx.param(gdxName = paste0(path.output.ce, '/gdxOutYear', y, '.gdx'), 
                          'pLife')
  
#  lifetimes <- read.csv(paste0(path.data.ce, '/NewPlantData/',
#                               'LifetimeValuesExistingPlants4Aug2016.csv'))  
  
  df.new.techs <- read.csv(paste0(path.output.ce,sprintf('newTechsCE%4d.csv', y)))
  df.new.techs <- df.new.techs %>%
    transmute(Type=TechnologyType, Fuel=FuelType,
              Cooling=as.character(Cooling.Tech),
              CAPEX=as.numeric(CAPEX.2012..MW.)/1e3,
              FixedOM=as.numeric(FOM.2012..MW.yr.)/1e3) %>%
    arrange(Type, Cooling) %>%
    mutate(Cooling=ifelse(Type %in% c('Nuclear', 'Solar PV', 'Wind'), 'other', Cooling)) %>%
    mutate(cool2=ifelse(Cooling == 'other', '',
                        ifelse(Cooling=='dry cooling', 'DC',
                               ifelse(Cooling=='recirculating', 'RC',
                                      ifelse(Cooling == 'once through',
                                             'OT', ''))))) %>%
    mutate(type.2=ifelse(Cooling == 'other',
                         gsub(' ', '.', Type),
                         gsub(' ', '.', paste0(Type, '.', cool2)))) %>%
    mutate(type.3=ifelse(cool2 != '', paste0(Type,'+',cool2), as.character(Type))) %>%
    left_join(lifetimes, by=c('type.3' = 'tech')) %>%
    select(-type.3) %>% dplyr::rename(n=pLife)

  # Fixed annualized costs of new plants
  vFixNew <- rgdx(gdxName = paste0(path.output.ce, 'gdxOutYear', y, '.gdx'), 
                  requestList = list(name="vIc", form='full', field=c('l')))$val
  # Variable annualized costs of ALL plants (existing and new)
  vVarTotal <- rgdx(gdxName = paste0(path.output.ce, 'gdxOutYear', y, '.gdx'), 
                    requestList = list(name="vVc", form='full', field=c('l')))$val
  
  fname <- paste0(path.output.ce, 'demandFullYrZonalCE', y,'.csv')
  df.dem <- read.csv(file=fname, stringsAsFactors = FALSE, header = FALSE) %>%
    select(-V1) %>% gather('hour', 'demand', -V2) %>% 
    mutate(hour=as.numeric(gsub('V', '', hour))-2) %>% arrange(V2, hour) %>%
    dplyr::rename(Region.Name=V2)
  
  annual.demand <- sum(df.dem$demand)
  
  # get existing generators
  fname <- paste0(path.output.ce, 'genFleetAfterCE', y,'.csv')
  df.egu <- read.csv(file=fname, stringsAsFactors = FALSE, header = TRUE) %>%
    select(Plant.Name, UniqueID_Final, ORIS.Plant.Code, Unit.ID, PlantType, Region.Name, 
           State.Name, Capacity..MW., Heat.Rate..Btu.kWh., Cooling.Tech, On.Line.Year, 
           Retirement.Year, Modeled.Fuels, FOM...MW.yr., VOM...MWh., YearAddedCE, 
           YearRetiredByCE, YearRetiredByAge) %>%
    mutate(Cooling.Tech=ifelse(grepl('once through', Cooling.Tech), 'once through',
                               ifelse(grepl('recirculating', Cooling.Tech), 'recirculating',
                                      ifelse(grepl('dry cooling', Cooling.Tech), 'dry cooling', 'other')))) %>%
    mutate(YearAddedCE=ifelse(is.na(YearAddedCE), 9999, YearAddedCE),
           YearRetiredByCE=ifelse(is.na(YearRetiredByCE), 9999, YearRetiredByCE),
           YearRetiredByAge=ifelse(is.na(YearRetiredByAge), 9999, YearRetiredByAge)) %>%
    # remove those retired by the model in years <= y
    filter(YearRetiredByCE > y & YearRetiredByAge > y)
  
  # existing generators not added by the model
  df.egu.original <- df.egu %>% filter(YearAddedCE == 9999) %>%
    mutate(annual.fix.cost=Capacity..MW.*(FOM...MW.yr.)/1e3)
    
  # existing generators added by the model (add capital costs)
  df.egu.added.ce <- df.egu %>% filter(YearAddedCE != 9999) %>% 
    left_join(df.new.techs, by=c('PlantType' = 'Type', 'Cooling.Tech' = 'Cooling')) %>% 
    mutate(annual.fix.cost=Capacity..MW.*(CAPEX*(rate*(1+rate)^n/((1+rate)^n-1))+FixedOM))
  
  # "annualized" levelized cost of energy of total fleet in year y
  
  df.out <-  data.frame(year=y,
                        type=c('Fixed Costs', 'Variable Costs'),
                        absolute.value=c(vFixNew+sum(df.egu.original$annual.fix.cost)+sum(df.egu.added.ce$annual.fix.cost), vVarTotal),
                        annual.demand=annual.demand) %>%
    mutate(value.per.MWh=absolute.value*1e3/annual.demand)
  
  return(df.out)
}


check.available.energy <- function(y, path.gdx) {

  #path.gdx <- "/Volumes/Seagate_RIPS/CE/results//rcp85//gdxOutYear2045.gdx"
  path.gdx <- "/Users/kiko/Documents/CE/out/CE/gdxOutYear2025.gdx"
  
  
  pDemand <- rgdx.param(gdxName = path.gdx, symName = "pDemand") %>% 
    group_by(g, h) %>% summarise(value = sum(pDemand)) %>%
    mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
    mutate(date=as.POSIXct(ymd("20200101") + hours(h2))) %>%
    mutate(season=define.season(date))
    
  ggplot(pDemand) + geom_line(aes(x=date, y=value), colour='red') + 
    guides(colour=FALSE) + theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0,'lines'),
          panel.border = element_rect(fill = NA)) + 
    facet_grid(cols = vars(season), scales = 'free_x')
  
  
  # check generation
  gexist <- rgdx(gdxName = path.gdx, requestList = list(name='vPegu', field='l'))
  gexist <- convert.gdx.result.df(gexist) %>% 
    group_by(g, h) %>% summarise(value = sum(value)) %>%
    mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
    mutate(date=as.POSIXct(ymd("20200101") + hours(h2))) %>%
    mutate(season=define.season(date)) %>%
    rename(value.exist = value)

  ggplot(gexist) + geom_line(aes(x=date, y=value.exist), colour='red') + 
    guides(colour=FALSE) + theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0,'lines'),
          panel.border = element_rect(fill = NA)) + 
    facet_grid(cols = vars(season), scales = 'free_x')
  
  gtechcurt <- rgdx(gdxName = path.gdx, requestList = list(name='vPtechcurtailed', field='l'))
  gtechcurt  <- convert.gdx.result.df(gtechcurt)
  
  if (!is.null(gtechcurt)){
    gtechcurt <- gtechcurt %>% 
      group_by(g, h) %>% summarise(value = sum(value)) %>%
      mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
      mutate(date=as.POSIXct(ymd("20200101") + hours(h2))) %>%
      mutate(season=define.season(date)) %>%
      rename(value.tech.curt = value)
  }
  
  gtechnotcurt <- rgdx(gdxName = path.gdx, requestList = list(name='vPtechnotcurtailed', field='l'))
  gtechnotcurt  <- convert.gdx.result.df(gtechnotcurt)

  if (!is.null(gtechnotcurt)){
    gtechnotcurt <- gtechnotcurt %>% 
      group_by(g, h) %>% summarise(value = sum(value)) %>%
      mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
      mutate(date=as.POSIXct(ymd("20200101") + hours(h2))) %>%
      mutate(season=define.season(date)) %>%
      rename(value.tech.not.curt = value)
  }
  
  gtechrenew <- rgdx(gdxName = path.gdx, requestList = list(name='vPtechrenew', field='l'))
  gtechrenew  <- convert.gdx.result.df(gtechrenew)
  
  if (!is.null(gtechrenew)){
    gtechrenew <- gtechrenew %>% 
      group_by(g, h) %>% summarise(value = sum(value)) %>%
      mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
      mutate(date=as.POSIXct(ymd("20200101") + hours(h2))) %>%
      mutate(season=define.season(date)) %>%
      rename(value.tech.renew = value)
  }

  ggplot(gtechrenew) + geom_line(aes(x=date, y=value.tech.renew), colour='red') + 
    guides(colour=FALSE) + theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0,'lines'),
          panel.border = element_rect(fill = NA)) + 
    facet_grid(cols = vars(season), scales = 'free_x')
  
  gtotal <- pDemand %>% select(g, date, value) %>% 
    left_join(gexist %>% select(g, date, value.exist), by=c('g', 'date'))

  if (!is.null(gtechcurt)) {
    gtotal <- gtotal %>% 
      left_join(gtechcurt %>% select(g, date, value.tech.curt), by=c('g', 'date'))
  }

  if (!is.null(gtechnotcurt)) {
    gtotal <- gtotal %>% 
      left_join(gtechnotcurt %>% select(g, date, value.tech.not.curt), by=c('g', 'date'))
  }

  if (!is.null(gtechrenew)) {
    gtotal <- gtotal %>% 
      left_join(gtechrenew %>% select(g, date, value.tech.renew), by=c('g', 'date'))
  }
  
  gtotal <- gtotal %>% mutate_at(.vars=vars(value.exist, value.tech.curt, value.tech.renew),
                                 .fun=function(x){ifelse(is.na(x), 0, x)}) %>%
    mutate(net=value.exist + value.tech.curt + value.tech.renew - value)
  
  range(gtotal$value.dem - (gtotal$value.genexist+gtotal$value))
  
  x <- gdxInfo(gdxName = paste0(path.gdx), dump=FALSE, returnList=TRUE, returnDF=FALSE)
  
  # get capacity from all existing generators
  pCapacExist <- rgdx.param(gdxName = paste0(path.gdx), symName = "pCapac")
  
  # get max energy from existing solar generators
  ids.solar <- rgdx.set(gdxName = paste0(path.gdx), symName = "solaregu") %>% 
    mutate(egu=as.character(egu)) %>% unlist()
  pMaxgensolar <- rgdx.param(gdxName = path.gdx, symName = "pMaxgensolar") %>%
    group_by(g, h) %>% summarise(value = sum(pMaxgensolar))
  pCapacSolar <- pCapacExist %>% filter(egu %in% ids.solar) %>% group_by(h) %>% summarise(value=sum(pCapac))

  pMaxgensolar <- pMaxgensolar %>% mutate(capac = mean(pCapacSolar$value),
                                          cf = value/capac) %>%
    mutate(h2 = as.numeric(gsub('h', '', as.character(h)))) %>%
    mutate(hour.of.day=hour(ymd("20200101") + hours(h2)))
  
  # get max energy from existing wind generators
  pMaxgenwind <- rgdx.param(gdxName = paste0(path.gdx), symName = "pMaxgenwind") %>%
    group_by(g, h) %>% summarise(value = sum(pMaxgenwind))
    
  # get ids from existing solar, wind, hydro, pumped hydro
  ids.wind <- rgdx.set(gdxName = paste0(path.gdx), symName = "windegu") %>% 
    mutate(egu=as.character(egu)) %>% unlist()
  ids.hydro <- rgdx.set(gdxName = paste0(path.gdx), symName = "hydroegu") %>% 
    mutate(egu=as.character(egu)) %>% unlist()
  ids.pumped <- rgdx.set(gdxName = paste0(path.gdx), symName = "pumphydroegu") %>% 
    mutate(egu=as.character(egu)) %>% unlist()
  
  ids.non.thermal <- c(ids.wind, ids.solar, ids.hydro, ids.pumped)
  
  
  # remove existing solar, wind, hydro, pumped hydro from data set
  capac.thermal <- rgdx.param(gdxName = paste0(path.gdx), symName = "pCapac") %>%
    mutate(egu=as.character(egu)) %>% filter(!(egu %in% ids.non.thermal)) %>% 
    group_by(g, h) %>% summarise(value = sum(pCapac))
  
  ggplot(capac.thermal) + geom_line(aes(x=h, y=value, colour=g, group=g))
  
  
  x <- pMaxgensolar %>% group_by(hour.of.day) %>% 
    summarise(cf.mean=mean(cf), lb=quantile(cf, 0.025), ub=quantile(cf, 0.975),
              sd=sd(cf))
  
  ggplot(x) + geom_line(aes(x=hour.of.day, y=cf.mean)) + geom_errorbar(aes(x=hour.of.day, ymin=lb, ymax=ub))
  ggplot(x) + geom_line(aes(x=hour.of.day, y=cf.mean)) + geom_errorbar(aes(x=hour.of.day, ymin=cf.mean - sd*1.96, ymax=cf.mean + sd*1.96))
  
}

convert.gdx.result.df <- function(list.gdx) {
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
  
  return (df.gdx)
  
}

define.season <- function(dt) {
  
  s <- ifelse(month(dt) %in% c(3,4,5), 'spring',
              ifelse(month(dt) %in% c(6,7,8), 'summer',
                     ifelse(month(dt) %in% c(9,10,11), 'fall', 'winter')))
  
  return(s)
  
} 


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
