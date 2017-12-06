# May 2017
# simplified analysis of impact on the generation cost due to the change
# in the demand

# get data frame with data from TVA power plants
source("supplycurve.R")

# function to create demand projection
source("rips_demand.R")

source("dispatch_LP.R")

utility.data <- get.utility.data()

df.data <- read.data.frame()
lm.model <- fit.lm.model(df=df.data)

# summer demand 2015
load.2015 <- (dplyr::filter(df.data, year==2015) %>% 
                group_by(year, season, hour.of.day) %>% 
                summarize(mean.value = mean(load, na.rm=TRUE)) %>% 
                filter(season == "Summer"))$mean.value

# load 2015 using simulated data from UW

load.2015.proj <- (project.demand(lm.hourly=lm.model, 
                                     year.proj=2015) %>% 
                        filter(season == "Summer"))$mean.value

sol1 <- run.LP(utility.data=utility.data, hourly.demand=load.2015)
obj.val1 <- compute.solution(sol1$solution, plants.data = utility.data,
                             H=24)
df.solution1 <- parse.solution(sol1$solution, plants.data = utility.data,
                               H=24)
g1 <- plot.dispatch.by.type(df.solution = df.solution1)
g1 <- g1 + ylim(0, 30000)

# summer demand 2040
load.2040 <- project.demand(lm.hourly=lm.model)
load.2040.summer <-  (load.2040 %>% filter(season == "Summer"))$mean.value

# comparing summer loads
compare.loads <- data.frame(type=c(rep("2015 hist", 24), 
                                   rep("2015 sim", 24),
                                   rep("2040 sim", 24)),
                            hour.of.day=rep(1:24, 3),
                            value=c(load.2015, load.2015.proj, 
                                    load.2040.summer))

shape.scale <- c(21, 22, 23, 24, 25)

g <- ggplot(data=compare.loads)
g <- g + geom_line(aes(x=hour.of.day, y=value, group=type))
g <- g + geom_point(aes(x=hour.of.day, y=value, shape=type), size=2, fill='white')
g <- g + scale_shape_manual(values = shape.scale)

path.plot <- '/Users/kiko/GoogleDrive/CMU/RIPS/reports/20170601/'

png(filename=paste(path.plot, 'compare_loads.png', sep = '/'))
print(g)
dev.off()

sol2 <- run.LP(utility.data=utility.data, hourly.demand=load.2040.summer)
obj.val2 <- compute.solution(sol2$solution, plants.data = utility.data,
                             H=24)
df.solution2 <- parse.solution(sol2$solution, plants.data = utility.data,
                              H = 24)
g2 <- plot.dispatch.by.type(df.solution = df.solution2)
g2 <- g2 + ylim(0, 30000)

g.compare <- plot.compare.total.gen.by.type(df.solution1, df.solution2, 
                                            case.names = c('2015', '2040'))

png(filename = paste(path.plot, 'curve2015.png', sep = '/'))
print(g1)
dev.off()

png(filename = paste(path.plot, 'curve2040.png', sep = '/'))
print(g2)
dev.off()

png(filename = paste(path.plot, 'compare.png', sep = '/'))
print(g.compare)
dev.off()

cost1 <- compute.cost.by.type(sol1$solution, utility.data)
cost1$case <- '2015'
cost2 <- compute.cost.by.type(sol2$solution, utility.data)
cost2$case <- '2040'

df <- rbind(cost1, cost2)

df$type <- factor(df$type, levels = rev(c("HYDRO", "NUCLEAR", 
                                          "COAL", "GAS", "OIL")))
df$case <- as.factor(df$case)

g <- ggplot(data = df) + 
  geom_bar(aes(x=case, y=value, fill=type), stat='identity', width = 0.6) +
  theme_bw()
g <- g + xlab("Case") + ylab("Cost (MM $)")
g

png(filename = paste(path.plot, 'compare_cost.png', sep = '/'))
print(g)
dev.off()

xx <- df %>% group_by(case) %>% summarize(total=sum(value))
sprintf('Change in generation cost: %5.2f%%',(xx$total[2]/xx$total[1]-1)*100)

file.copy(from='./supplycurve.png', 
          to=paste(path.plot, 'supplycurve.png', sep='/'))
