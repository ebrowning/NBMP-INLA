## INLA models - using presence absence at spots/walks 

## 
setwd("C:/Users/ellab/Dropbox/PhD/NBMP_data")
library(data.table); library(INLA); library(rgdal);library(dplyr);library(sf);library(raster);library(fields); library(viridis)
source('funcs_plotting_inla.R')
# read in data used in GAM models 
field <- fread('field_chpt3.csv')
field <- field[, -c(1)] # remove stupid V1 columns
colnames(field)

# set factor and numeric columns
str(field)
field[, c(1,6,7,8,10,11,16,31,48,41,49, 50,54,55,67)] <- lapply(field[, c(1,6,7,8,10, 11,16,31,48,49,49, 50, 54,55,67)], as.factor)
#field[, c(18:21, 23:26)] <- lapply(field[, c(18:21, 23:26)], as.numeric)

# Make new year column where the first year is 1 rather than 1997
field$CountYear2 <- field$CountYear - 1997

# add intercept 
field$intercept <- 1

# remove surveys with fewer than 12 spots
field <- field[field$max_spot >11, ]
table(field$max_spot)

# Make a modified site column, which combines detector and ssites. At some repeatedly surveyed sites the detector used changes, and could be a reason for the increases observed in some species
## unique site and detector combinations = modified site
field$mod_ssites <- as.factor(as.numeric(as.factor(paste(field$ssites, field$Detector, sep = "_"))))

x = count(field, "mod_ssites")
nrow(x[x$freq == 1,])


field$log10_prop_forest <- log10(field$prop_forest + 1)
field$log10_prop_broad <- log10(field$prop_broad + 1)
field$log10_prop_needle <- log10(field$prop_needle + 1)
field$log10_prop_grass <- log10(field$prop_grass + 1)
field$log10_prop_agri <- log10(field$prop_agri + 1)
field$log10_prop_urban <- log10(field$prop_urban + 1)
field$log10_prop_water <- log10(field$prop_water + 1)
field$log10_prop_other <- log10(field$prop_other + 1)


# set hyperparamters for fixed effects (gaussian assumed)
hyper.fix = list(theta1 = list(initial = log(100), fixed = TRUE))
# NB 'pc.prec' = penalized complexity prior'
#
my.init = NULL
inla.setOption("enable.inla.argument.weights", TRUE)
# family 
family1 = 'binomial'
control.family1 = list(control.link=list(model="logit"))
## iid hyper-pars
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1))) 

######################
## Binomial models - number of spots or walks where a bat was counted

##########
# Pipistrellus pipistrellus

f_pip <- pip_spot_count ~ -1 + intercept + scale(Temperature) + WindStrength + 
  VolSkillSelfAssessment  + mins_after_set + Detector_GroupID +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))

M1_pip <- inla(f_pip, 
               data = field,
               family = family1,
               control.family = control.family1,
               Ntrials = 12,
               weights = field$weight, 
               control.fixed=list(prec=1),
               control.inla = list(int.strategy='eb' ),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M1_pip)
plot(M1_pip)

M1_pip$cpo$failure[M1_pip$cpo$failure > 0]

plot(field2$pip_spot_count, inla.link.invlogit(M1_pip$summary.fitted.values$`0.5quant`))
# plot the distribution of random effects
plot(density(M1_pip$summary.random$CountYear$mean))
lines(0+c(-2, 2)*M1_pip$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")

plot(density(M1_pip$summary.random$ssites$mean))
lines(0+c(-2, 2)*M1_pip$summary.hyperpar$`0.5quant`[2]^(-0.5) , c(0,0), col="blue")


###############

f_pyg <- pyg_spot_count ~ -1 + intercept + scale(Temperature) + WindStrength  + Detector_GroupID +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))

M1_pyg <- inla(f_pyg, 
               data = field,
               family = 'binomial',
               control.family = control.family1,
               Ntrials = 12,
               weights = field$weight,
               control.fixed=list(expand.factor.strategy = 'inla',prec=1),
               control.inla = list(int.strategy='eb'),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M1_pyg)
plot(M1_pyg)

plot(field$pyg_spot_count, inla.link.invlogit(M1_pyg$summary.fitted.values$`0.5quant`))
# plot the distribution of random effects
plot(density(M1_pyg$summary.random$CountYear$mean))
lines(0+c(-2, 2)*M1_pyg$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")

plot(density(M1_pyg$summary.random$ssites$mean))
lines(0+c(-2, 2)*M1_pyg$summary.hyperpar$`0.5quant`[2]^(-0.5) , c(0,0), col="blue")

##############

f_noc <- noc_walk_count ~ -1 + intercept + scale(Temperature) + scale(Duration) + 
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M1_noc <- inla(f_noc, 
               data = field,
               family = 'binomial',
               control.family = control.family1,
               Ntrials = 12,
               weights = field$weight,
               control.fixed=list(expand.factor.strategy = 'inla',prec=1),
               control.inla = list(int.strategy='eb'),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M1_noc)
plot(M1_noc)

plot(field$noc_walk_count, inla.link.invlogit(M1_noc$summary.fitted.values$`0.5quant`))
# plot the distribution of random effects
plot(density(M1_noc$summary.random$CountYear$mean))
lines(0+c(-2, 2)*M1_noc$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")

plot(density(M1_noc$summary.random$ssites$mean))
lines(0+c(-2, 2)*M1_noc$summary.hyperpar$`0.5quant`[2]^(-0.5) , c(0,0), col="blue")

##############
### subset field data to Eng/wal
field_ew <- field[!field$country == 'scot', ]

f_ser <- ser_walk_count ~ -1 + intercept + scale(Temperature) + scale(Duration) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))
## iid hyper-pars
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5))) 
M1_ser <- inla(f_ser, 
               data = field_ew,
               family = 'binomial',
               control.family = control.family1,
               Ntrials = 12,
               control.fixed=list(expand.factor.strategy = 'inla',prec=1),
               control.inla = list(int.strategy='eb', strategy = "laplace", npoints = 21),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M1_ser)
plot(M1_ser)
# plot the distribution of random effects
plot(density(M1_ser$summary.random$CountYear2$mean))
lines(0+c(-2, 2)*M1_ser$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")

plot(density(M1_ser$summary.random$ssites$mean))
lines(0+c(-2, 2)*M1_ser$summary.hyperpar$`0.5quant`[2]^(-0.5) , c(0,0), col="blue")


#########################################################
#########################################################
# --- Accounting for space using the SPDE
###
## Make meshes, SPDEs and stacks 

library(rgdal)
shape6 <- readOGR(dsn = ".",
                  layer = "GB_coarse_osgb")
#convert to km from metres using the sf package 
projection(shape6)
# transfrom into an af object
shape6 <- st_as_sf(shape6)
# change the units from metres to km 
shape6 <- st_transform(shape6, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=km +no_defs +ellps=airy +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894")
# convert back to a spatial polygon
shape6 <- as(shape6, 'Spatial')
plot(shape6)

max.edge = 5
mesh_poly = inla.mesh.2d(boundary = shape6, 
                          max.edge = c(1, 5)*max.edge,
                          cutoff = 0.5, 
                         offset = c(10, 20))
plot(mesh_poly, main = '')


# Set range prioirs for the spde model
prior.median.sd = 1; prior.median.range = 50
spde = inla.spde2.pcmatern(mesh_poly, prior.range = c(prior.median.range, 0.5),
                            prior.sigma = c(prior.median.sd, 0.1), constr = T)

indexs <- inla.spde.make.index("s", spde$n.spde)

locs = data.matrix(field[ , c('easting', 'northing')]/1000)
points(locs, col = 'red')

d= pointDistance(cbind(field$easting, field$northing), lonlat = FALSE)

 A = inla.spde.make.A(mesh = mesh_poly, loc = locs)

## E. serotinus requires a separate mesh and SPDE to the other three species as only surveyed for in England/Wales
### subset field data to Eng/wal
field_ew <- field[!field$country == 'scot', ]
### get shapefile for Eng/Wal
shape7 <- readOGR(dsn = ".", 
                  layer = "EngWal_coarse")
# transform to osgb
shape7 <- spTransform(shape7, CRS("+init=epsg:27700"))
 
plot(shape7)

max.edge = 5
mesh_poly2 = inla.mesh.2d(boundary = shape7, 
                          max.edge = c(1, 5)*max.edge,
                          cutoff = 0.5, 
                          offset = c(10, 20))
plot(mesh_poly2, main = '')

d= pointDistance(cbind(field_ew$easting, field_ew$northing), lonlat = FALSE)
hist(d/1000)
# Set range prioirs for the spde model
prior.median.sd = 1; prior.median.range = 50
spde2 = inla.spde2.pcmatern(mesh_poly2, prior.range = c(prior.median.range, 0.5),
                            prior.sigma = c(prior.median.sd, 0.1), constr = T)
#index for stack
indexs2 <- inla.spde.make.index("s", spde2$n.spde)
# A matrix and stack
locs2 = data.matrix(field_ew[, c("easting", "northing")]/1000)


A2 = inla.spde.make.A(mesh = mesh_poly2, loc = locs2)

#################
# -- Make stacks for each species

colnames(field) # check which column is needed for each species
mod_covs <- c("intercept", "Temperature", "Duration", "mins_after_set", "WindStrength", "VolSkillSA_Poor", "Detector_GroupID",
              "log10_prop_broad", "log10_prop_needle", "log10_prop_grass", 
              "log10_prop_agri", "log10_prop_water", "log10_prop_urban", 
              "tas_2_1km", "tas_4_1km", "precip_1_1km", 
              "CountYear2", "mod_ssites", "ssites")

## P. pipistrellus
stack_pip <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = field$pip_spot_count), 
                        effects=list(s=indexs, # spatial
                                     data.frame(field[ ,..mod_covs])),
                        A =list(A,1))


## P. pygmaeus
stack_pyg <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = field$pyg_spot_count), 
                        effects=list(s=indexs, # spatial
                                     data.frame(field[ ,..mod_covs])),
                        A =list(A,1))
## N. noctula
stack_noc <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = field$noc_walk_count),
                        effects=list(s=indexs, # spatial
                                     data.frame(field[ ,..mod_covs])),
                        A =list(A,1))
## E. serotinus
stack_ser <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = field_ew$ser_walk_count), 
                        effects=list(s=indexs2, # spatial
                                     data.frame(field_ew[ ,..mod_covs])),
                        A =list(A2,1))


#################################
# -- Run models

my.init = NULL

## P. pipistrellus 

f_pip6 <- count ~ -1 + intercept + scale(Temperature) + mins_after_set + 
  WindStrength + VolSkillSA_Poor + Detector_GroupID +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde)

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))

M6_pip <- inla(f_pip6, 
               data = inla.stack.data(stack_pip),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(prec=1),
               control.inla = list(int.strategy='eb', npoints = 21),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_pip), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M6_pip)
plot(M6_pip)

cov_plot(M6_pip$summary.fixed, ylims = c(-2, 1), plot_title = 'P. pipistrellus', y_text_pos = -2)
yr_bintrend_plot(M6_pip$summary.random$CountYear2, ylims=c(50,220), plot_title = '')
local.plot.field(M6_pip$summary.random$s$mean, mesh_poly, zlim = c(-3,3))

plot(density(M6_pip$summary.random$CountYear2$mean))
lines(0+c(-2, 2)*M6_pip$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")

## P. pygmaeus


hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
f_pyg6 <- count ~ -1 + intercept + scale(Temperature) + 
  VolSkillSA_Poor + Detector_GroupID +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde)

M6_pyg <- inla(f_pyg6, 
               data = inla.stack.data(stack_pyg),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(prec=1),
               control.inla = list(int.strategy='eb', npoints = 21),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_pyg), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M6_pyg)
plot(M6_pyg)

local.plot.field(M6_pyg$summary.random$s$mean, mesh_poly, zlim= c(-3,3))
cov_plot(M6_pyg$summary.fixed, ylims = c(-1, 1), plot_title = 'P. pygmaeus', y_text_pos = -1)
yr_bintrend_plot(M6_pyg$summary.random$CountYear2, ylims=c(50,220), plot_title = '')

plot(density(M6_pyg$summary.random$CountYear2$mean))
lines(0+c(-2, 2)*M6_pyg$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")


### N. noctula 
f_noc6 <- count ~ -1 + intercept + scale(Temperature)+ scale(Duration) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)+
  f(s, model = spde)

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M6_noc <- inla(f_noc6, 
               data = inla.stack.data(stack_noc),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(prec=1),
               control.inla = list(int.strategy='eb', npoints = 21),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_noc), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE), verbose = TRUE)
summary(M6_noc)
plot(M6_noc)

local.plot.field(M6_noc$summary.random$s$mean, mesh_poly, zlim = c(-3,3))
cov_plot(M6_noc$summary.fixed, ylims = c(-1, 1), plot_title = 'N. noctula', y_text_pos = -1)
yr_bintrend_plot(M6_noc$summary.random$CountYear2, ylims=c(50,220), plot_title = '')


### E serotinus 

f_ser6 <- count ~ -1 + intercept + scale(Temperature)+ scale(Duration) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)+
  f(s, model = spde2)

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M6_ser <- inla(f_ser6, 
               data = inla.stack.data(stack_ser),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(prec=1),
               control.inla = list(int.strategy='eb'),
               control.mode=list(restart=TRUE, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_ser), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M6_ser)
plot(M6_ser)

local.plot.field(M6_ser$summary.random$s$mean, mesh_poly2, zlim = c(-3,3))
cov_plot(M6_ser$summary.fixed, ylims = c(-3, 3), plot_title = 'E. serotinus', y_text_pos = -3.5)
yr_bintrend_plot(M6_ser$summary.random$CountYear2, ylims=c(50,220), plot_title = '')

############################################################################################################
############################################################################################################
# -- Acounting for envrionment

## P. pip
f_pip5 <- count ~ -1 + intercept + scale(Temperature) + mins_after_set + 
  WindStrength + VolSkillSA_Poor + Detector_GroupID +
  scale(log10_prop_broad) + scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + scale(log10_prop_urban) +
  scale(tas_2_1km) + 
  scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde) 

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))

M5_pip <- inla(f_pip5, 
               data = inla.stack.data(stack_pip),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(expand.factor.strategy = 'inla', prec=1),
               num.threads=1,
               control.inla = list(int.strategy='eb',npoints = 21),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_pip), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_pip)
plot(M5_pip)

cov_plot(M5_pip$summary.fixed, ylims = c(-1, 2), plot_title = 'P. pipistrellus', y_text_pos = -1)
yr_bintrend_plot(M5_pip$summary.random$CountYear, ylims=c(50,220), plot_title = '')
local.plot.field(M5_pip$summary.random$s$mean, mesh_poly)#, zlim = c(0,5))

plot(density(M5_pip$summary.random$s$mean))
lines(0+c(-2, 2)*M5_pip$summary.hyperpar$`0.5quant`[3]^(-0.5) , c(0,0), col="blue")

## P. pygmaeus


f_pyg5 <- count ~ -1 + intercept + scale(Temperature) + 
  #WindStrength + 
  VolSkillSA_Poor + Detector_GroupID +
  scale(log10_prop_broad) + scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + 
  scale(log10_prop_urban) +
  scale(tas_2_1km) + scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde) 

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))

M5_pyg <- inla(f_pyg5, 
               data = inla.stack.data(stack_pyg),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(expand.factor.strategy = 'inla', prec=1),
               num.threads=1,
               control.inla = list(int.strategy='eb'),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_pyg), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_pyg)
plot(M5_pyg)

local.plot.field(M5_pyg$summary.random$s$mean, mesh_poly, zlim= c(-3,3))
cov_plot(M5_pyg$summary.fixed, ylims = c(-1, 1), plot_title = 'P. pygmaeus', y_text_pos = -1)
yr_bintrend_plot(M5_pyg$summary.random$CountYear, ylims=c(50,150), plot_title = '')

### N. noctula 

f_noc5 <- count ~ -1 + intercept + scale(Temperature)+ scale(Duration) +
  scale(log10_prop_broad) + #scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + #scale(log10_prop_urban) +
  #scale(tas_2_1km) + scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M5_noc <- inla(f_noc5, 
               data = inla.stack.data(stack_noc),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(prec=1),
               num.threads=1,
               control.inla = list(int.strategy='eb'),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_noc), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_noc)
plot(M5_noc)

local.plot.field(M5_noc$summary.random$s$mean, mesh_poly, zlim = c(-3,3))
cov_plot(M5_noc$summary.fixed, ylims = c(-1, 1), plot_title = 'N. noctula', y_text_pos = -1)
yr_bintrend_plot(M5_noc$summary.random$CountYear2, ylims=c(50,220), plot_title = '')

plot(density(M5_noc$summary.random$CountYear$mean))
lines(0+c(-2, 2)*M5_noc$summary.hyperpar$`0.5quant`[1]^(-0.5) , c(0,0), col="blue")

tmp = inla.tmarginal(function(x) x, M5_noc$marginals.hyperpar[[3]]) 
plot(tmp, type = "l", xlab = expression(sigma[u]), ylab = "Density")

### E serotinus 

f_ser4 <- count ~ -1 + intercept + scale(Temperature) + scale(Duration) +# WindStrength +
  # scale(log10_prop_broad) + #scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + #scale(log10_prop_urban) +
  scale(tas_2_1km) + scale(tas_4_1km) + #scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)+
  f(s, model = spde2) 

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M5_ser <- inla(f_ser4, 
               data = inla.stack.data(stack_ser),
               family = 'binomial',
               Ntrials = 12,
               control.fixed=list(prec=1),
               control.inla = list(int.strategy='eb'),#, npoints = 21),
               control.mode=list(restart=TRUE, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_ser), link = 1, compute=T),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_ser)
plot(M5_ser)

local.plot.field(M5_ser$summary.random$s$mean, mesh_poly2, zlim = c(-3,3))
cov_plot(M5_ser$summary.fixed, ylims = c(-1, 1), plot_title = 'E. serotinus', y_text_pos = -1)
yr_bintrend_plot(M5_ser$summary.random$CountYear, ylims=c(50,220), plot_title = '')

save(M1_pip, M1_pyg, M1_noc, M1_ser, M5_pip, M5_pyg, M5_noc, M5_ser, M6_pip, M6_pyg, M6_noc, M6_ser, file = "chpt3_inla_binomial_mods_2.RData")

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

## -- Disaggregating trends to country level 

#### -- England only models

f_eng <- field[field$country == 'eng', ]

## make new mesh for eng


library(rgdal)
shape8 <- readOGR(dsn = ".",
                  layer = "eng_coarse")
#convert to km from metres using the sf package 
projection(shape8)
# transfrom into an af object
shape8 <- st_as_sf(shape8)
# change the units from metres to km 
shape8 <- st_transform(shape8, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=km +no_defs +ellps=airy +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894")
# convert back to a spatial polygon
shape8 <- as(shape8, 'Spatial')
plot(shape8)

max.edge = 10
mesh_poly = inla.mesh.2d(boundary = shape8, 
                         max.edge = c(1, 5)*max.edge,
                         cutoff = 0.5, 
                         offset = c(10, 20))
plot(mesh_poly, main = '')


# Set range prioirs for the spde model
prior.median.sd = 1; prior.median.range = 50
spde = inla.spde2.pcmatern(mesh_poly, prior.range = c(prior.median.range, 0.5),
                           prior.sigma = c(prior.median.sd, 0.5), constr = T)
indexs <- inla.spde.make.index("s", spde$n.spde)

locs = data.matrix(f_eng[ , c('easting', 'northing')]/1000)
points(locs, col = 'red')

A = inla.spde.make.A(mesh = mesh_poly, loc = locs)
colnames(f_eng)
mod_covs <- c("intercept", "Temperature", "Duration", "mins_after_set", "WindStrength", "VolSkillSA_Poor", "Detector_GroupID",
              "log10_prop_broad", "log10_prop_needle", "log10_prop_grass", 
              "log10_prop_agri", "log10_prop_water", "log10_prop_urban", 
              "tas_2_1km", "tas_4_1km", "precip_1_1km", 
              "CountYear2", "mod_ssites", "ssites")

## P. pipistrellus
stack_pip <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_eng$pip_spot_count ), 
                        effects=list(s=indexs, # spatial
                                     data.frame(f_eng[ ,..mod_covs])),
                        A =list(A,1))
## P. pygmaeus
stack_pyg <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_eng$pyg_spot_count),
                        effects=list(s=indexs, # spatial
                                     data.frame(f_eng[ ,..mod_covs])),
                        A =list(A,1))
## N. noctula
stack_noc <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_eng$noc_walk_count ),
                        effects=list(s=indexs, # spatial
                                     data.frame(f_eng[ ,..mod_covs])),
                        A =list(A,1))
## E. serotinus 
stack_ser <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_eng$ser_walk_count), 
                        effects=list(s=indexs, # spatial
                                     data.frame(f_eng[ ,..mod_covs])),
                        A =list(A,1))

## models
f_pip5 <- count ~ -1 + intercept + scale(Temperature) + mins_after_set + 
  WindStrength + VolSkillSA_Poor + Detector_GroupID +
  scale(log10_prop_broad) + #scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + scale(log10_prop_urban) +
  #scale(tas_2_1km) + 
  scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde) 

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))

M5_eng_pip <- inla(f_pip5, 
                   data = inla.stack.data(stack_pip),
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_pip), link = 1, compute=T),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_eng_pip)
plot(M5_eng_pip)

cov_plot(M5_eng_pip$summary.fixed, ylims = c(-1, 2), plot_title = 'P. pipistrellus', y_text_pos = -1)
yr_bintrend_plot(M5_eng_pip$summary.random$CountYear, ylims=c(50,220), plot_title = '')
local.plot.field(M5_eng_pip$summary.random$s$sd, mesh_poly, zlim = c(0,3))


## P. pygmaeus

f_pyg5 <- count ~ -1 + intercept + scale(Temperature) + 
  #WindStrength + 
  VolSkillSA_Poor + Detector_GroupID +
  scale(log10_prop_broad) + scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + 
  scale(log10_prop_urban) +
  scale(tas_2_1km) + scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde) 

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))

M5_eng_pyg <- inla(f_pyg5, 
                   data = inla.stack.data(stack_pyg),
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb',npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_pyg), link = 1, compute=T),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_eng_pyg)
plot(M5_eng_pyg)

local.plot.field(M5_eng_pyg$summary.random$s$sd, mesh_poly, zlim= c(0,3))
cov_plot(M5_eng_pyg$summary.fixed, ylims = c(-1, 1), plot_title = 'P. pygmaeus', y_text_pos = -1)
yr_bintrend_plot(M5_eng_pyg$summary.random$CountYear, ylims=c(50,220), plot_title = '')

### N. noctula 
f_noc5 <- count ~ -1 + intercept + scale(Temperature)+ scale(Duration) +
  scale(log10_prop_broad) + #scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + #scale(log10_prop_urban) +
  #scale(tas_2_1km) + scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M5_eng_noc <- inla(f_noc5, 
                   data = inla.stack.data(stack_noc),
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb',npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_noc), link = 1, compute=T),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_eng_noc)
plot(M5_eng_noc)

local.plot.field(M5_eng_noc$summary.random$s$sd, mesh_poly, zlim = c(0,3))
local.plot.field(M5_eng_noc$summary.random$s$mean, mesh_poly, zlim = c(-3,3))
cov_plot(M5_eng_noc$summary.fixed, ylims = c(-1, 1), plot_title = 'N. noctula', y_text_pos = -1)
yr_bintrend_plot(M5_eng_noc$summary.random$CountYear2, ylims=c(50,220), plot_title = '')

### E serotinus 

f_ser <- ser_walk_count ~ -1 + intercept + scale(Temperature) + scale(Duration) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M1_eng_ser <- inla(f_ser, 
                   data = f_eng,
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb',  npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(compute = TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M1_eng_ser)
plot(M1_eng_ser)

cov_plot(M1_eng_ser$summary.fixed, ylims = c(-2, 1), plot_title = 'E. serotinus', y_text_pos = -1)
yr_bintrend_plot(M1_eng_ser$summary.random$CountYear, ylims=c(50,220), plot_title = '')

#############################################
#############################################

### -- Scotland only 
f_scot <- field[field$country == 'scot', ]


## make new mesh for scot


shape9 <- readOGR(dsn = ".",
                  layer = "scot_simple_main")
#convert to km from metres using the sf package 
projection(shape9)
# transfrom into an af object
shape9 <- st_as_sf(shape9)
# change the units from metres to km 
shape9 <- st_transform(shape9, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=km +no_defs +ellps=airy +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894")
# convert back to a spatial polygon
shape9 <- as(shape9, 'Spatial')
plot(shape9)

max.edge = 10
mesh_poly = inla.mesh.2d(boundary = shape9, 
                         max.edge = c(1, 5)*max.edge,
                         cutoff = 0.5, 
                         offset = c(10, 20))
plot(mesh_poly, main = '')

# Set range prioirs for the spde model
prior.median.sd = 1; prior.median.range = 50
spde = inla.spde2.pcmatern(mesh_poly, prior.range = c(prior.median.range, 0.5),
                           prior.sigma = c(prior.median.sd, 0.5), constr = T)
indexs <- inla.spde.make.index("s", spde$n.spde)

locs = data.matrix(f_scot[ , c('easting', 'northing')]/1000)
points(locs, col = 'red')

A = inla.spde.make.A(mesh = mesh_poly, loc = locs)

colnames(f_scot)

## P. pipistrellus
stack_pip <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_scot$pip_spot_count), 
                        effects=list(s=indexs, # spatial
                                     data.frame(f_scot[ ,..mod_covs])),
                        A =list(A,1))
## P. pygmaeus
stack_pyg <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_scot$pyg_spot_count),
                        effects=list(s=indexs, # spatial
                                     data.frame(f_scot[ ,..mod_covs])),
                        A =list(A,1))
## N. noctula
stack_noc <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                        data = list(count = f_scot$noc_walk_count),
                        effects=list(s=indexs, # spatial
                                     data.frame(f_scot[ ,..mod_covs])),
                        A =list(A,1))



## models

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M5_scot_pip <- inla(f_pip5, 
                   data = inla.stack.data(stack_pip),
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_pip), link = 1, compute=T),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_scot_pip)
plot(M5_scot_pip)

cov_plot(M5_scot_pip$summary.fixed, ylims = c(-1, 2), plot_title = 'P. pipistrellus', y_text_pos = -1)
yr_bintrend_plot(M5_scot_pip$summary.random$CountYear, ylims=c(50,220), plot_title = '')
local.plot.field(M5_scot_pip$summary.random$s$sd, mesh_poly, zlim = c(0,3))
local.plot.field(M5_scot_pip$summary.random$s$mean, mesh_poly, zlim = c(-3,3))


## P. pygmaeus

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M5_scot_pyg <- inla(f_pyg5, 
                   data = inla.stack.data(stack_pyg),
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb',npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_pyg), link = 1, compute=T),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_scot_pyg)
plot(M5_scot_pyg)

local.plot.field(M5_scot_pyg$summary.random$s$sd, mesh_poly, zlim= c(0,3))
local.plot.field(M5_scot_pyg$summary.random$s$mean, mesh_poly, zlim= c(-3,3))
cov_plot(M5_scot_pyg$summary.fixed, ylims = c(-1, 1), plot_title = 'P. pygmaeus', y_text_pos = -1)
yr_bintrend_plot(M5_scot_pyg$summary.random$CountYear, ylims=c(50,220), plot_title = '')

### N. noctula 

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))

M5_scot_noc <- inla(f_noc5, 
                   data = inla.stack.data(stack_noc),
                   family = 'binomial',
                   control.family = control.family1,
                   Ntrials = 12,
                   control.fixed=list(expand.factor.strategy = 'inla',prec=1),
                   control.inla = list(int.strategy='eb',npoints = 21),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_noc), link = 1, compute=T),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,  config = TRUE))
summary(M5_scot_noc)
plot(M5_scot_noc)

local.plot.field(M5_scot_noc$summary.random$s$sd, mesh_poly, zlim = c(0,3))
local.plot.field(M5_scot_noc$summary.random$s$mean, mesh_poly, zlim = c(-3,3))
cov_plot(M5_scot_noc$summary.fixed, ylims = c(-1, 1), plot_title = 'N. noctula', y_text_pos = -1)
yr_bintrend_plot(M5_scot_noc$summary.random$CountYear2, ylims=c(50,220), plot_title = '')

save(M5_eng_noc, M5_eng_pip, M5_eng_pyg, M1_eng_ser, M5_scot_noc, M5_scot_pip, M5_scot_pyg, file = "chpt3_binmods_disag_2.RData")



