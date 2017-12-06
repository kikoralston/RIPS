library(maps)
library(mapdata)

library(ggplot2)
library(plyr)
library(dplyr)

usa <- map_data("county")

usa2 <- usa %>% filter(region %in% c('mississippi', 'alabama', 'georgia',
                                     'north carolina', 'tennessee', 
                                     'virginia', 'kentucky',
                                     'missouri', 'arkansas')) 

# coord_map("mercator",xlim = c(-90.5, -80),ylim = c(32, 38))

g <- ggplot() + geom_polygon(data = usa2, 
                             aes(x=long, y = lat, group = group), 
                             fill='white', color = "gray") + 
  coord_fixed(1.3425, xlim = c(-90.7, -80), ylim = c(32, 38), expand = FALSE) + 
  theme_bw() +
  scale_y_continuous(breaks=seq(30, 39, by=0.5)) +
  scale_x_continuous(breaks=seq(-95, -75, by=1))

img <- readJPEG(source = '~/Downloads/TVA-Service-Area-Map.jpg')
gg <- rasterGrob(image = img, gp = gpar(alpha=0.9))

g <- g + annotation_custom(gg)

g <- g + geom_polygon(data = usa2, 
                      aes(x=long, y = lat, group = group), 
                      fill='white', color = "gray", alpha=0.4)
