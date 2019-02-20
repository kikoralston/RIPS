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

plot.map.serc <- function(path.data.ce, path.output.ce) {
  
  s1 <- readOGR(paste0(path.data.ce, 'IPMRegionsShapeFile/'), 
                'IPM_Regions_20121023_US')
  
  serc.zones <- read.csv(paste0(path.output.ce, 'zoneNamesToNumbers.csv'))
  
  zones.map <- ifelse(s1@data$IPM_Region %in% serc.zones$Zone, 
                      s1@data$IPM_Region, NA)
  
  s1.union.serc <- unionSpatialPolygons(s1 ,zones.map)
  
  df.ipm.counties <- s1@data %>% 
    mutate(FIPS=as.numeric(as.character(FIPS)))
  
  map.us.counties <- map_data('county')
  
  map.us.counties <- map.us.counties %>% 
    mutate(polyname=paste(region, subregion, sep=',')) %>%
    left_join(county.fips, by='polyname') %>%
    rename(FIPS=fips) %>% left_join(df.ipm.counties, by='FIPS')
  
  map.us.counties <- map.us.counties %>% 
    mutate(IPM_Region = as.character(IPM_Region)) %>%
    mutate(IPM_Region = ifelse(IPM_Region %in% serc.zones$Zone,
                               IPM_Region,
                               NA))
  
  state_map <- map_data('state')
  
  # state boundaries
  state_map <- map_data("state")
  
  pal <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
           '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  
  g <- ggplot() +
    geom_polygon(data=map.us.counties, aes(x=long, y=lat, group=group), 
                 colour='gray80') +
    geom_polygon(data=map.us.counties, aes(x=long, y=lat, group=group, 
                                           fill=IPM_Region)) +
    geom_path(data=state_map, aes(x=long, y=lat, group=group)) +
    coord_map("polyconic", xlim=c(-95, -75), ylim=c(27.5, 40)) + 
    scale_fill_manual(breaks=serc.zones$Zone, values=pal, na.value="gray90") +
    theme_bw() + theme(axis.title.x=element_blank(), 
                       axis.title.y=element_blank()) + 
    guides(fill=guide_legend(title=NULL))
  
  pdf("./slides/sercmap.pdf", width=7*16/9, height=7)
  print(g)
  dev.off()
}

plot.inital.fleet <- function(path.output.ce, path.data.ce, df.demand, 
                              path.out.plot=path.output.ce) {
  
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
    summarise(online.cap = sum(online.cap)) %>% arrange(currYear) %>%
    ungroup() %>%
    mutate(PlantType=factor(PlantType, levels=rev(plant.types))) %>%
    mutate(Region.Name == as.character(Region.Name))
  
  df.demand.aux <- df.demand %>% gather(key = "type", value = "value", 
                                    mean.load:peak.load)

  list.cases <- unique(df.demand$case)
  ncases <- length(list.cases)
  
  g <- ggplot() + geom_col(data=df.fleet.time.2,
                           aes(x=currYear, y=online.cap/1e3, fill=PlantType),
                           width=2) +
    geom_line(data=df.demand.aux, aes(x=currYear, y=value/1e3, linetype=case,
                                      color=type), size=1)
  
  g <- g + scale_fill_manual(values = col.pallete) + theme_bw() +
    xlab('Year') + ylab('Installed Capacity (GW)') +
    scale_x_continuous(breaks = seq(2020, 2050, 5)) +
    guides(fill=guide_legend(title=NULL),
           linetype=guide_legend(title=NULL, reverse=TRUE)) +
    facet_grid(Region.Name ~ ., scales = 'free_y') 
    
  pdf(paste0(path.out.plot,"/supply_demand.pdf"), width=7*16/9, height=7)
  print(g)
  dev.off()
  
  return(df.fleet.time.2)
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
  df.new.techs <- read.csv(paste0(path.rcps[1],sprintf('newTechsCE%4d.csv',
                                                       yearIni)))
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
      # map cell to zone
      map.cell2zone <- rgdx.param(gdxName = paste0(rcp, 'gdxOutYear', y, 
                                                   '.gdx'), 'pCellzones')
      map.cell2zone$c <- as.character(map.cell2zone$c)
      
      tech.not.curtailed <- rgdx(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
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
      
      tech.curtailed <- rgdx(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
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
        summarise(value=sum(value, na.rm=TRUE)) %>% mutate(year=y) %>%
        select(year, Zone, type, value) %>% as.data.frame()
      
      tech.renew <- rgdx(gdxName = paste0(rcp, '/gdxOutYear', y, '.gdx'), 
                         requestList = list(name="vNrenew", 
                                            form='full',field=c('l')))
      
      df.tech.renew <- data.frame(zone=tech.renew$uels$z,
                                  tech.renew$val, row.names = NULL) %>%
        mutate(ZoneNum=as.numeric(as.factor(zone))) %>%
        select(-zone) %>%
        gather(key='type', value='value', -ZoneNum) %>%
        left_join(serc.zones, by='ZoneNum') %>%
        arrange(ZoneNum) %>% mutate(year=y) %>%
        select(year, Zone, type, value) %>% as.data.frame()
      
      df.complete.decision <- rbind(df.complete.decision,
                                    data.frame(case=case.name,
                                               df.tech.curtailed.zone, 
                                               stringsAsFactors = FALSE),
                                    data.frame(case=case.name,
                                               df.tech.not.curtailed, 
                                               stringsAsFactors = FALSE),
                                    data.frame(case=case.name,
                                               df.tech.renew,
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
    offset <- data_frame(case=list.cases,
                         off=c(-width/2, width/2))
  } else if (ncases == 3) {
    width <- 1
    offset <- data_frame(case=list.cases,
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
          axis.title.x = element_blank()) +
    facet_grid(Zone ~ ., scales = 'free_y')

  g_param <- ggplot_build(g)
  x.range <- g_param$layout$panel_params[[1]]$x.range
  year.values <- unique(df.complete.decision.3$year)
  
  df.top.axis <- df.complete.decision.3 %>%
    mutate(xpos=year) %>%
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
    rename(online.cap=total.cap, PlantType=Type, currYear=year, Region.Name=Zone) %>%
    select(case, Region.Name, PlantType, currYear, online.cap) %>%
    group_by(case, Region.Name, PlantType, currYear) %>%
    summarise(online.cap=sum(online.cap))
  
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
          rename(Region.Name=V2)
        
        df <- df %>% group_by(Region.Name) %>% 
          summarise(mean.load=mean(demand), peak.load=max(demand)) %>%
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
      summarize(CF=mean(pCf)) %>% mutate(year=y) %>%
      select(year, Zone, techrenew, hour.of.day, CF) %>% as.data.frame() 
    
    w <- x %>% group_by(year, Zone, techrenew) %>% summarise(mean.CF=sum(CF)/24)
    
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
  
  df.fleet <- df.fleet %>% group_by(Modeled.Fuels, PlantType, cooling) %>% summarize(sum=sum(Capacity..MW.)) %>% as.data.frame()
  
  total <- sum(df.fleet$sum)
  
  lvl0 <- data.frame(name = paste0("SERC\n", formatC(total, big.mark=',', 
                                                     digits=0, 
                                                     format = 'f'),
                                   ' MW'), 
                     value = total, level = 0, fill = NA) %>%
    mutate(xmin=0, xmax=1, ymin=0, ymax=value) %>%
    mutate(x.avg=0, y.avg=0, colour=FALSE, label=name)
  
  lvl1 <- df.fleet %>% group_by(Modeled.Fuels) %>% summarize(sum=sum(sum)) %>%
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
    summarize(value=sum(sum)) %>%
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
    summarize(value=sum(sum)) %>%
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
    mutate(label=ifelse(value<400 | is.na(name), '', name))
  
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
    mutate(hjust=ifelse(name=='DC' & fill == 'Natural Gas', 0, hjust),
           vjust=ifelse(name=='DC' & fill == 'Natural Gas', 1, vjust)) %>%
    mutate(hjust=ifelse(name=='Other', 1, hjust),
           vjust=ifelse(name=='Other', 0, vjust)) %>%
    mutate(hjust=ifelse(name=='Solar', 0.2, hjust),
           vjust=ifelse(name=='Solar', 0.9, vjust))
  
  g <- ggplot(data=df.final, aes(fill = fill)) +  
    geom_rect(aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, color=colour), 
              size = 0.1) +
    scale_fill_manual(values = new.pallete.2, na.translate = FALSE) +
    scale_color_manual(values = c('TRUE'='gray20', 'FALSE'='#FFFFFF00'), 
                       guide = F, na.translate = FALSE) +
    geom_text(aes(x = x.avg, y = y.avg, label = label, hjust=hjust, 
                  vjust=vjust), size = rel(2.5)) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() + guides(fill=FALSE) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "lines"),
          panel.grid = element_blank(),
          legend.position = 'none') +
    coord_polar(theta = "y")
  
#  pdf("~/CMU/RIPS/paperCE/longdraft/images/sunburst.pdf", width=5, height=5)
#  print(g)
#  dev.off()
#  system('pdfcrop ~/CMU/RIPS/paperCE/longdraft/images/sunburst.pdf ~/CMU/RIPS/paperCE/longdraft/images/sunburst.pdf')
  
  return(g)
}
