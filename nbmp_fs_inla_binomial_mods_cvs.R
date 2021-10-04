#######################################################################

### Leave-one-out cross-validation

#######################################################################

# Load libraries
pacman::p_load("dplyr", "rgdal", "tidyr", "data.table",
               "rgdal", "INLA", "sf", "sp",
               "MASS", "Metrics", "INLAutils")

setwd("C:/Users/ellab/Dropbox/PhD/NBMP_data")


#Read in data - only relevant columns 

field <- fread('field_chpt3.csv', select = c('ssites', 'mod_ssites', 'CountYear',  'CountYear2', 'country', 'reg', 'easting', 'northing', 'intercept', 'Temperature', 'WindStrength', 
                                             'VolSkillSA_Poor', 'mins_after_set', 'Duration', "Detector_GroupID",
                                             'log10_prop_broad', 'log10_prop_needle', 'log10_prop_grass', 'log10_prop_agri', 'log10_prop_urban', 
                                             'tas_2_1km', 'tas_4_1km', 'precip_1_1km', 'pip_spot_count', 'pyg_spot_count', 'noc_walk_count', 'ser_walk_count', 'max_spot', 'max_walk'))

# remove surveys with fewer than 12 spots
field <- field[field$max_spot >11, ]
table(field$max_spot)

str(field)
field[, c(1,2,5,6,11,12,15)] <- lapply(field[, c(1,2,5,6,11,12,15)], as.factor)

#Define regions 
reg = unique(field$reg)

# Define meshes from shapefiles

## GB mesh
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



### - England Wales mesh 
shape7 <- readOGR(dsn = ".", 
                  layer = "EngWal_lesscoarse_outline")
# transform to osgb
shape7 <- spTransform(shape7, CRS("+init=epsg:27700"))
#convert to km from metres using the sf package 
# transfrom into an sf object
shape7 <- st_as_sf(shape7)
# change the units from metres to km 
shape7 <- st_transform(shape7, "+init=epsg:27700 +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=km +no_defs +ellps=airy +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894")
# convert back to a spatial polygon
shape7 <- as (shape7, 'Spatial')
plot(shape7)

max.edge = 5
mesh_poly2 = inla.mesh.2d(boundary = shape7, 
                          max.edge = c(1, 5)*max.edge,
                          cutoff = 0.5, 
                          offset = c(10, 20))
plot(mesh_poly2, main = '')


#Things for all models 
my.init = NULL
control.family2 = list(control.link=list(model="log")) 
## iid hyper-pars
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1))) 
pip_hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))
pyg_hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))
noc_hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.5, 0.5)))
ser_hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.1)))



##################################################
##################################################
##################################################
## M1

# Formulas 
f_pip <- y ~ -1 + intercept + scale(Temperature) + WindStrength + 
  VolSkillSA_Poor  + mins_after_set + Detector_GroupID + 
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(mod_ssites, model = "iid", hyper=hyper.iid, constr=T)

f_pyg <- y ~ -1 + intercept + scale(Temperature) + VolSkillSA_Poor  + Detector_GroupID + 
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(mod_ssites, model = "iid", hyper=hyper.iid, constr=T)

f_noc <- y ~ -1 + intercept + scale(Temperature) + scale(Duration) + 
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)

f_ser <- y ~ -1 + intercept + scale(Temperature) + scale(Duration) + WindStrength +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)

#Loop through regions
for (i in 1:length(reg)){
  print(reg[i])
  # assign data
  data <- field
  
  #### Loop through species models
  species <- c("P. pipistrellus", "P. pygmaeus", "N. noctula", "E. serotinus")
  
  for (j in 1:length(species)){
    print(species[j])
    species_x <- species[j]
    # condition y and formula on species and set filename for saving outputs
    if(species_x=="P. pipistrellus"){
      data$y = data$pip_spot_count
      mod_filename = "pip"
      form_x = f_pip}
    if(species_x =="P. pygmaeus"){
      data$y = data$pyg_spot_count
      mod_filename = "pyg" 
      form_x = f_pyg}
    if(species_x=="N. noctula"){
      data$y = data$noc_walk_count
      mod_filename = "noc" 
      form_x = f_noc }
    if(species_x == "E. serotinus"){
      data$y = data$ser_walk_count
      mod_filename = "ser"
      form_x = f_ser }
    
    # condition subset of data and mesh on species - England/Wales only for E.ser
    if(j=="E. serotinus"){
      df = data[!data$county == 'scot']
      n = nrow(df)
    } else {
      df = data
      n = nrow(df)
    }
    
    
    # assign new df for INLA modelling and set y cases to NA to predict
    df_x <- df
    df_x$y[df_x$reg == reg[i]] <- NA
    
    # Run model
    mod <- inla(form_x, 
                data = df_x,
                family = 'binomial',
                Ntrials = 12,
                control.fixed=list(expand.factor.strategy = 'inla', prec=1),
                control.inla = list(int.strategy='eb', npoints = 21),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(compute=T),
                control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE))
    print(mod$waic$waic)
    
    # Save model
    save(mod, file = paste("./cv_models/", mod_filename, "M1", reg[i], ".RData", sep = "_"))
    
    ## Model posterior predictive distribution
    s <- 1000
    samples = inla.posterior.sample(1000, mod)
    samples.m <- sapply(samples, function(x) x$latent[1:n])
    samples.m <- exp(samples.m)
    kappa <- as.vector(sapply(samples, function(x) x$hyperpar[1]))
    pred.samples=matrix(NA,nrow=dim(samples.m)[1],ncol=s)
    
    if(j == "P. pipistrellus") {
      
      for (l in 1:dim(samples.m)[1])
      {
        # sample from neg bin
        pred.samples[l,]=rpois(s,samples.m[l,],kappa)
      }
    } else {
      
      for (l in 1:dim(samples.m)[1])
      {
        # sample from neg bin
        pred.samples[l,]=rnegbin(s,samples.m[l,],kappa)
      }
    }
    
    model_pred <- data.frame(Species    = rep(j, nrow(df)),
                             Model      = rep(i, nrow(df)), 
                             CountYear  = df$CountYear,
                             reg        = df$reg,
                             ssites     = df$ssites,
                             easting    = df$easting,
                             northin    = df$northing, 
                             Observed   = df$y,
                             Predicted  = apply(pred.samples,1,mean),
                             lci        = apply(pred.samples,1,quantile,probs=c(0.025),na.rm=T),
                             uci        = apply(pred.samples,1,quantile,probs=c(0.975),na.rm=T)) %>%
      # Subset to predicted region
      subset(reg == reg[i])
    
    save(model_pred, file = paste("./cv_models/predictions/", mod_filename, "M1", reg[i], "pred", ".RData", sep = "_"))
    
  }
}

##################################################
##################################################
##################################################

## M6

# Formulas 
f_pip6 <- count ~ -1 + intercept + scale(Temperature) + mins_after_set + 
  WindStrength + VolSkillSA_Poor + Detector_GroupID +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde_x)

f_pyg6 <- count ~ -1 + intercept + scale(Temperature) + 
  VolSkillSA_Poor + Detector_GroupID +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde_x)

f_noc6 <- count ~ -1 + intercept + scale(Temperature)+ scale(Duration) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)+
  f(s, model = spde_x)

f_ser6 <- count ~ -1 + intercept + scale(Temperature) + scale(Duration) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T) +
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)+
  f(s, model = spde_x)
#Loop through regions
for (i in 1:length(reg)){
  print(reg[i])
  # assign data
  data <- field
  
  #### Loop through species models
  species <- c("P. pipistrellus", "P. pygmaeus", "N. noctula", "E. serotinus")
  
  for (j in 1:length(species)){
    print(species[j])
    species_x <- species[j]
    # condition y and formula on species and set filename for saving outputs
    if(species_x=="P. pipistrellus"){
      data$y = data$pip_spot_count
      mod_filename = "pip"
      form_x = f_pip6
      hyper.rw = pip_hyper.rw}
    if(species_x =="P. pygmaeus"){
      data$y = data$pyg_spot_count
      mod_filename = "pyg" 
      form_x = f_pyg6
      hyper.rw = pyg_hyper.rw}
    if(species_x=="N. noctula"){
      data$y = data$noc_walk_count
      mod_filename = "noc" 
      form_x = f_noc6 
      hyper.rw = noc_hyper.rw}
    if(species_x == "E. serotinus"){
      data$y = data$ser_walk_count
      mod_filename = "ser"
      form_x = f_ser6 
      hyper.rw = ser_hyper.rw}
    
    # condition subset of data and mesh on species - England/Wales only for E.ser
    if(j=="E. serotinus"){
      df = data[!data$county == 'scot']
      n = nrow(df)
      mesh_x = mesh_poly2
    } else {
      df = data
      n = nrow(df)
      mesh_x = mesh_poly
    }
    
    
    # condition model family on species - P pip = poisson, else nbinomial
    if(j == "P.pipistrellus"){
      mod_fam = "Poisson"
    } else {
      mod_fam = "nbinomial"
    }
    
    # assign new df for INLA modelling and set y cases to NA to predict
    df_x <- df
    df_x$y[df_x$reg == reg[i]] <- NA
    # create A matrix
    locs = cbind(df_x$easting/1000, df$northing/1000) # using km so not such large numbers
    A_x = inla.spde.make.A(mesh = mesh_x, loc = locs)
    
    # Set range prioirs for the spde model
    prior.median.sd = 1; prior.median.range = 50
    spde_x = inla.spde2.pcmatern(mesh_x, prior.range = c(prior.median.range, 0.5),
                                prior.sigma = c(prior.median.sd, 0.5), constr = T)
    
    indexs <- inla.spde.make.index("s", spde_x$n.spde)
    
    # create stack
    stack_x <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                          data = list(count = df_x$y), 
                          effects=list(s=indexs, # spatial
                                       data.frame(df_x[,-"y"])),
                          A =list(A_x,1))
    
    # Run model
    mod <- inla(form_x, 
                data = inla.stack.data(stack_x),
                family = mod_fam,
                Ntrials = 12,
                control.fixed=list(expand.factor.strategy = 'inla', prec=1),
                control.inla = list(int.strategy='eb', npoints = 21),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_x), link = 1, compute=T),
                control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE))
    print(mod$waic$waic)
    
    # Save model
    save(mod, file = paste("./cv_models/", mod_filename, "M6", reg[i], ".RData", sep = "_"))
    
    ## Model posterior predictive distribution
    s <- 1000
    samples = inla.posterior.sample(1000, mod)
    samples.m <- sapply(samples, function(x) x$latent[1:n])
    samples.m <- exp(samples.m)
    kappa <- as.vector(sapply(samples, function(x) x$hyperpar[1]))
    pred.samples=matrix(NA,nrow=dim(samples.m)[1],ncol=s)
    
    if(j == "P. pipistrellus") {
      
      for (l in 1:dim(samples.m)[1])
      {
        # sample from neg bin
        pred.samples[l,]=rpois(s,samples.m[l,],kappa)
      }
    } else {
      
      for (l in 1:dim(samples.m)[1])
      {
        # sample from neg bin
        pred.samples[l,]=rnegbin(s,samples.m[l,],kappa)
      }
    }
    
    model_pred <- data.frame(Species    = rep(j, nrow(df)),
                             Model      = rep(i, nrow(df)), 
                             CountYear  = df$CountYear,
                             reg        = df$reg,
                             ssites     = df$ssites,
                             easting    = df$easting,
                             northin    = df$northing, 
                             Observed   = df$y,
                             Predicted  = apply(pred.samples,1,mean),
                             lci        = apply(pred.samples,1,quantile,probs=c(0.025),na.rm=T),
                             uci        = apply(pred.samples,1,quantile,probs=c(0.975),na.rm=T)) %>%
      # Subset to predicted region
      subset(reg == reg[i])
    
    save(model_pred, file = paste("./cv_models/predictions/", mod_filename, "M6", reg[i], "pred", ".RData", sep = "_"))
    
  }
}

##################################################
##################################################
##################################################



# M5

f_pip5 <- count ~ -1 + intercept + scale(Temperature) + mins_after_set + 
  WindStrength + VolSkillSA_Poor + Detector_GroupID +
  scale(log10_prop_broad) + #scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + scale(log10_prop_urban) +
  #scale(tas_2_1km) + 
  scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde_x) 

f_pyg5 <- count ~ -1 + intercept + scale(Temperature) + 
  #WindStrength + 
  VolSkillSA_Poor + Detector_GroupID +
  scale(log10_prop_broad) + scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + scale(log10_prop_urban) +
  scale(tas_2_1km) + scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde_x)

f_noc5 <- count ~ -1 + intercept + scale(Temperature)+ scale(Duration) +
  scale(log10_prop_broad) + #scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + #scale(log10_prop_urban) +
  #scale(tas_2_1km) + scale(tas_4_1km) + scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T) +
  f(s, model = spde_x)

f_ser5 <- count ~ -1 + intercept + scale(Temperature) + scale(Duration) +# WindStrength +
  #scale(log10_prop_broad) + scale(log10_prop_needle)+ #scale(log10_prop_forest) +
  scale(log10_prop_grass) + scale(log10_prop_agri) + #scale(log10_prop_urban) +
  scale(tas_2_1km) + scale(tas_4_1km) + #scale(precip_1_1km) +
  f(CountYear2, model = "rw2", hyper = hyper.rw, constr = T, scale.model = T)+ 
  f(ssites, model = "iid", hyper=hyper.iid, constr=T)+
  f(s, model = spde_x)

#Loop through regions
for (i in 1:length(reg)){
  print(reg[i])
  # assign data
  data <- field
  
  #### Loop through species models
  species <- c("P. pipistrellus", "P. pygmaeus", "N. noctula", "E. serotinus")
  
  for (j in 1:length(species)){
    print(species[j])
    species_x <- species[j]
    # condition y and formula on species and set filename for saving outputs
    if(species_x=="P. pipistrellus"){
      data$y = data$pip_spot_count
      mod_filename = "pip"
      form_x = f_pip5
      hyper.rw = pip_hyper.rw}
    if(species_x =="P. pygmaeus"){
      data$y = data$pyg_spot_count
      mod_filename = "pyg" 
      form_x = f_pyg5
      hyper.rw = pyg_hyper.rw}
    if(species_x=="N. noctula"){
      data$y = data$noc_walk_count
      mod_filename = "noc" 
      form_x = f_noc5 
      hyper.rw = noc_hyper.rw}
    if(species_x == "E. serotinus"){
      data$y = data$ser_walk_count
      mod_filename = "ser"
      form_x = f_ser5 
      hyper.rw = ser_hyper.rw}
    
    # condition subset of data and mesh on species - England/Wales only for E.ser
    if(j=="E. serotinus"){
      df = data[!data$county == 'scot']
      n = nrow(df)
      mesh_x = mesh_poly2
    } else {
      df = data
      n = nrow(df)
      mesh_x = mesh_poly
    }
    
    
    # condition model family on species - P pip = poisson, else nbinomial
    if(j == "P.pipistrellus"){
      mod_fam = "Poisson"
    } else {
      mod_fam = "nbinomial"
    }
    
    # assign new df for INLA modelling and set y cases to NA to predict
    df_x <- df
    df_x$y[df_x$reg == reg[i]] <- NA
    # create A matrix
    locs = cbind(df_x$easting/1000, df$northing/1000) # using km so not such large numbers
    A_x = inla.spde.make.A(mesh = mesh_x, loc = locs)
    
    # Set range prioirs for the spde model
    prior.median.sd = 1; prior.median.range = 50
    spde_x = inla.spde2.pcmatern(mesh_x, prior.range = c(prior.median.range, 0.5),
                                 prior.sigma = c(prior.median.sd, 0.1), constr = T)
    
    indexs <- inla.spde.make.index("s", spde_x$n.spde)
    
    # create stack
    stack_x <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                          data = list(count = df_x$y), 
                          effects=list(s=indexs, # spatial
                                       data.frame(df_x[,-"y"])),
                          A =list(A_x,1))
    
    # Run model
    mod <- inla(form_x, 
                data = inla.stack.data(stack_x),
                family = mod_fam,
                Ntrials = 12,
                control.fixed=list(expand.factor.strategy = 'inla', prec=1),
                control.inla = list(int.strategy='eb', npoints = 21),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_x), link = 1, compute=T),
                control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE))
    print(mod$waic$waic)
    
    # Save model
    save(mod, file = paste("./cv_models/", mod_filename, "M5", reg[i], ".RData", sep = "_"))
    
    ## Model posterior predictive distribution
    s <- 1000
    samples = inla.posterior.sample(1000, mod)
    samples.m <- sapply(samples, function(x) x$latent[1:n])
    samples.m <- exp(samples.m)
    kappa <- as.vector(sapply(samples, function(x) x$hyperpar[1]))
    pred.samples=matrix(NA,nrow=dim(samples.m)[1],ncol=s)
    
    if(j == "P. pipistrellus") {
      
      for (l in 1:dim(samples.m)[1])
      {
        # sample from neg bin
        pred.samples[l,]=rpois(s,samples.m[l,],kappa)
      }
    } else {
      
      for (l in 1:dim(samples.m)[1])
      {
        # sample from neg bin
        pred.samples[l,]=rnegbin(s,samples.m[l,],kappa)
      }
    }
    
    model_pred <- data.frame(Species    = rep(j, nrow(df)),
                             Model      = rep(i, nrow(df)), 
                             CountYear  = df$CountYear,
                             reg        = df$reg,
                             ssites     = df$ssites,
                             easting    = df$easting,
                             northin    = df$northing, 
                             Observed   = df$y,
                             Predicted  = apply(pred.samples,1,mean),
                             lci        = apply(pred.samples,1,quantile,probs=c(0.025),na.rm=T),
                             uci        = apply(pred.samples,1,quantile,probs=c(0.975),na.rm=T)) %>%
      # Subset to predicted regio n
      subset(reg == reg[i])
    
    save(model_pred, file = paste("./cv_models/predictions/", mod_filename, "M5", reg[i], "pred", ".RData", sep = "_"))
    
  }
}


