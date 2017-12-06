source('readUW.R')

year.ini <- 1990
year.end <- 2099

path.in <- '~/Downloads/'
path.out <- '~/Documents/GCM/'

name.zip <- 'Jackson.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='MS', path.out=path.out)

name.zip <- 'Nashville.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='TN', path.out=path.out)

name.zip <- 'Montgomery.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='AL', path.out=path.out)

name.zip <- 'Atlanta.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='GA', path.out=path.out)

name.zip <- 'Columbia.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='SC', path.out=path.out)

name.zip <- 'Charlotte.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='NC', path.out=path.out)

name.zip <- 'Louisville.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='KY', path.out=path.out)

name.zip <- 'Little Rock.zip' 
convert.txt2rds(name.zip = file.path(path.in, name.zip), yinf = year.ini, ysup=year.end, prefix='AR', path.out=path.out)
