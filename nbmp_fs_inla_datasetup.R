setwd("C:/Users/ellab/Dropbox/PhD/NBMP_data")

library(dplyr);library(plyr); library(data.table); library(lubridate);library(maptools);library(rgdal);library(raster);library(tidyr)

field <- fread("NBMP_Field_Survey_data.csv")
f_sites <- fread("field_sites.csv")
#merge the two data frames so that field also has long/lat
field <- merge(field, f_sites, by = "ssites")

#############################################################
###!!!###!!!###!!!###!!!######!!!###!!!###!!!###!!!###
## A survey has been input twice - Site ID 120558, countdate 30-Jul-15
# Data from a single spot/walk has also been input 6 times - ID 85712
library(dplyr)
# first remove extra spot/walk rows for ID 85712
fx <- distinct(field, ID, .keep_all = T)

# remove ID column as for the twice input survey the rows have unique IDs, and then do same above 
fx <- fx[, -2]
fx <- distinct(fx)
# check duplicates have been removed 
x = fx[fx$ssites == 120558,]

field <- fx

##############################################################

# Need to work out time survey started after sunset
# make column with date and start/end times together
field$date_starttime = paste(field$CountDate, field$StartTime, sep = " ")
field$date_endtime = paste(field$CountDate, field$EndTime, sep = " ")
#convert to unix time
field$start_unix = as.numeric(as.POSIXct(field$date_starttime, format = '%d-%b-%y %H:%M:%S'))
field$end_unix = as.numeric(as.POSIXct(field$date_endtime, format = '%d-%b-%y %H:%M:%S'))
head(field)

# get sunset time for each day
field$date = as.POSIXct(strptime(field$CountDate, "%d-%b-%y"), tz = "GMT")
head(field)

locs = cbind(field$longitude, field$latitude)
sunset = crepuscule(locs, field$date, direction = 'dusk', solarDep = 1, POSIXct.out = T)
field$sunset = sunset$time
head(field)
#convert to unix
field$sunset_unix = strftime(field$sunset, '%s')

# calculate difference in start time and sunset time
field$mins_after_set = as.integer((field$start_unix - as.numeric(field$sunset_unix))/60)
head(field)


field <- ddply(field, .(ssites, longitude, latitude, GridReference, CountYear, ObserverID, VolExperience, VolSkillSelfAssessment, CountDate, Detector,
                        DetectorID, Temperature, Rainfall, Duration, mins_after_set, WindStrength), summarise,
               pip_count = sum(CommonPipCount, na.rm = T), 
               pip_spot_count = length(CommonPipCount[CommonPipCount > 0]), # number of spots a pipistrellus was observed
               pyg_count = sum(SopranoPipCount, na.rm = T), 
               pyg_spot_count = length(SopranoPipCount[SopranoPipCount > 0]), # number of spots a pygmaeus was observed
               max_spot = max(SpotNo, na.rm = T), 
               noc_count = sum(NoctuleCount, na.rm = T), 
               noc_walk_count = length(NoctuleCount[NoctuleCount >0]), # number of walks a noctule was observed
               ser_count = sum(SerotineCount, na.rm = T), 
               ser_walk_count = length(SerotineCount[SerotineCount >0]), # number of walks a serotine was observed
               leis_count = sum(LeislerCount, na.rm = T), 
               max_walk = max(WalkNo, na.rm = T))
head(field)
nrow(field) # 1576 surveys
# How many surveys had fewer than 12 walks/ spots?
table(field$max_spot)
table(field$max_walk)
#Are these the same surveys?
field[field$max_spot < 12,] # yes
# remove surveys where spot/walk numbers are less than 12 (only two)
field <- field[field$max_spot > 11, ]

# How many surveys have either misrecorded durations or took much longer than the standard 90 mins?
hist(field$Duration)
field[field$Duration > 300, ]
# remove these surveys with erroneous/ overlong  durations
field <- field[field$Duration < 300, ]

# remove NA columns (not sure where these have come from)
field <- field[!is.na(field$ssites), ]

nrow(field) # now 1543 surveys 

## Add country column 
field$country <- ifelse(field$ssites < 200000, "eng", 
                        ifelse((field$ssites > 199999 & field$ssites < 300000), "scot",
                               ifelse(field$ssites > 299999, "wal", NA )))

# Add easting/ northing column 
# use osg_parse function 
library(rnrfa)

east_north <- osg_parse(field$GridReference, coord_system = c("BNG", "WGS84"))

field$easting <- east_north$easting
field$northing <- east_north$northing


#####################################
## Correct the easting/northing coordinates of three sites in scotland
field <- fieldx
head(field)

y = field[!field$ssites == 220016 & !field$ssites ==220072 & !field$ssites ==220102, ]


x <- field[field$ssites == 220016, ]
table(x$ssites)

x$easting = 315000
x$northing = 929000

field2 <- rbind(y, x)

x <- field[field$ssites == 220072, ]
table(x$ssites)

x$easting = 334000
x$northing = 961000

field2 <- rbind(field2, x)

x <- field[field$ssites == 220102, ]
table(x$ssites)

x$easting = 232000
x$northing = 642000

field2 <- rbind(field2, x)

nrow(field2)

##########################################
# European space agency lancd over data 
## 1km rasters created in 'esa_prop_funcs.R'
# read in all esa landcover tifs

esa_prop_tifs <- list.files("C:/Users/ellab/Dropbox/PhD/NBMP_data/esa_lc_1km_aggregate_tifs/", pattern = "esa_prop", full.names = TRUE)

years <- c('X1998', 'X1999', 'X2000', 
           'X2001', 'X2002', 'X2003', 'X2004', 'X2005', 
           'X2006', 'X2007', 'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015')

esa_prop_forest_1km <- stack(esa_prop_tifs[[3]])
names(esa_prop_forest_1km) <- years
esa_prop_broad_1km <- stack(esa_prop_tifs[[2]])
names(esa_prop_broad_1km) <- years
esa_prop_needle_1km <- stack(esa_prop_tifs[[5]])
names(esa_prop_needle_1km) <- years
esa_prop_grass_1km <- stack(esa_prop_tifs[[4]])
names(esa_prop_grass_1km) <- years
esa_prop_agri_1km <- stack(esa_prop_tifs[[1]])
names(esa_prop_agri_1km) <- years
esa_prop_urban_1km <- stack (esa_prop_tifs[[7]])
names(esa_prop_urban_1km) <- years
esa_prop_water_1km <- stack (esa_prop_tifs[[6]])
names(esa_prop_water_1km) <- years
esa_prop_other_1km <- stack (esa_prop_tifs[[8]])
names(esa_prop_other_1km) <- years






## Extract at 1km buffers
f_sp <- f_sites
coordinates(f_sp) <- ~ longitude + latitude
proj4string(f_sp) <- proj4string(esa_prop_agri_1km)
  
  
forest_ex = raster::extract(esa_prop_forest_1km, f_sp)
forest_ex <- as.data.frame(forest_ex)
forest_ex$ssites = f_sites$ssites
forest_ex2 <- gather(forest_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(forest_ex2)[3] <- c('prop_forest')

broad_ex = raster::extract(esa_prop_broad_1km, f_sp)
broad_ex <- as.data.frame(broad_ex)
broad_ex$ssites = f_sites$ssites
broad_ex2 <- gather(broad_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(broad_ex2)[3] <- c('prop_broad')

needle_ex = raster::extract(esa_prop_needle_1km, f_sp)
needle_ex <- as.data.frame(needle_ex)
needle_ex$ssites = f_sites$ssites
needle_ex2 <- gather(needle_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(needle_ex2)[3] <- c('prop_needle')


agri_ex = raster::extract(esa_prop_agri_1km, f_sp)
agri_ex <- as.data.frame(agri_ex)
agri_ex$ssites = f_sites$ssites
agri_ex2 <- gather(agri_ex, key = CountYear, value = value, -ssites) # change rows into columns
head(agri_ex2)
colnames(agri_ex2)[3] <- c('prop_agri')

grass_ex = raster::extract(esa_prop_grass_1km, f_sp)
grass_ex <- as.data.frame(grass_ex)
grass_ex$ssites = f_sites$ssites
grass_ex2 <- gather(grass_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(grass_ex2)[3] <- c('prop_grass')

water_ex = raster::extract(esa_prop_water_1km, f_sp)
water_ex <- as.data.frame(water_ex)
water_ex$ssites = f_sites$ssites
water_ex2 <- gather(water_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(water_ex2)[3] <- c('prop_water')

urb_ex = raster::extract(esa_prop_urban_1km, f_sp)
urb_ex <- as.data.frame(urb_ex)
urb_ex$ssites = f_sites$ssites
urb_ex2 <- gather(urb_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(urb_ex2)[3] <- c('prop_urban')

other_ex = raster::extract(esa_prop_other_1km, f_sp)
other_ex <- as.data.frame(other_ex)
other_ex$ssites = f_sites$ssites
other_ex2 <- gather(other_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(other_ex2)[3] <- c('prop_other')

# bind together into one data frame
prop_1km <- as.data.frame(cbind(agri_ex2$ssites, agri_ex2$CountYear, forest_ex2$prop_forest, broad_ex2$prop_broad, needle_ex2$prop_needle, grass_ex2$prop_grass, agri_ex2$prop_agri, 
                                water_ex2$prop_water, urb_ex2$prop_urban, other_ex2$prop_other))
colnames(prop_1km) <- c("ssites", "CountYear", "prop_forest", "prop_broad", "prop_needle", "prop_grass", "prop_agri", "prop_water", "prop_urban", "prop_other")
prop_1km$CountYear <- substr(prop_1km$CountYear, 2,5) # remove 'X' from each year
prop_1km[, c(3:9)] <- lapply(prop_1km[, c(3:9)], as.numeric) 
prop_1km[,c(1:2)] <-  lapply(prop_1km[, c(1:2)], as.integer)
head(prop_1km)
# merge with field 
field3 <- left_join(field2, prop_1km, by = c('ssites', 'CountYear'))


############################
## Weather variables 

#read in stack with relevant annual weather vars (see XXXX.R for how this was made)
chess <- stack('chess_seas_osgrd/env_covs.grd')
plot(chess)


# subset out specific seasonal var stacks
### precip = precipitation ; sfc = surface wind speed ; tas = air temperature
### 1 = spring ; 2 = summer ; 4 = winter

precip_1 <- chess[[10:28]]
sfc_2 <- chess[[29:47]]
tas_2 <- chess[[(48:66)]]
tas_4 <- chess[[(67:85)]]

## extract for each site, all years

### 1km 
# Precip spring
precip_ex <- raster::extract(precip_1, f_sp, buffer = 1000, fun = mean)
precip_ex <- as.data.frame(precip_ex)
colnames(precip_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
precip_ex$ssites <- f_sites$ssites
precip_ex2 <- gather(precip_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(precip_ex2)[3] <- c('precip_1_1km')

# tas summer 
tas_ex <- raster::extract(tas_2, f_sp, buffer = 1000, fun = mean)
tas_ex <- as.data.frame(tas_ex)
colnames(tas_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
tas_ex$ssites <- f_sites$ssites
tas_ex2 <- gather(tas_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(tas_ex2)[3] <- c('tas_2_1km')

# tas winter 
tas_4_ex <- raster::extract(tas_4, f_sp, buffer = 1000, fun = mean)
tas_4_ex <- as.data.frame(tas_4_ex)
colnames(tas_4_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
tas_4_ex$ssites <- f_sites$ssites
tas_4_ex2 <- gather(tas_4_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(tas_4_ex2)[3] <- c('tas_4_1km')

#sfc wind summer 
sfc_ex <- raster::extract(sfc_2, f_sp, buffer = 1000, fun = mean)
sfc_ex <- as.data.frame(sfc_ex)
colnames(sfc_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
sfc_ex$ssites <- f_sites$ssites
sfc_ex2 <- gather(sfc_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(sfc_ex2)[3] <- c('sfc_2_1km')

# bind prop_x columns together
w_vars_1km <- as.data.frame(cbind(precip_ex2$ssites, precip_ex2$CountYear, precip_ex2$precip_1_1km, sfc_ex2$sfc_2_1km, tas_ex2$tas_2_1km, tas_4_ex2$tas_4_1km))
colnames(w_vars_1km) <- c("ssites", "CountYear", "precip_1_1km", "sfc_2_1km",  "tas_2_1km", "tas_4_1km")
w_vars_1km$CountYear <- substr(w_vars_1km$CountYear, 2,5) # remove 'X' from year
w_vars_1km[, c(2:6)] <- lapply(w_vars_1km[, c(2:6)], as.numeric) 
w_vars_1km$ssites <- as.integer(w_vars_1km$ssites)
head(w_vars_1km)
# merge with year
field3 <- left_join(field3, w_vars_1km, by = c('ssites', 'CountYear'))
(proc.time() - t)

head(field3)

### 3km 
# Precip spring
precip_ex <- raster::extract(precip_1, f_sp, buffer = 3000, fun = mean)
precip_ex <- as.data.frame(precip_ex)
colnames(precip_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
precip_ex$ssites <- f_sites$ssites
precip_ex2 <- gather(precip_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(precip_ex2)[3] <- c('precip_1_3km')

# tas summer 
tas_ex <- raster::extract(tas_2, f_sp, buffer = 3000, fun = mean)
tas_ex <- as.data.frame(tas_ex)
colnames(tas_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
tas_ex$ssites <- f_sites$ssites
tas_ex2 <- gather(tas_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(tas_ex2)[3] <- c('tas_2_3km')

# tas winter 
tas_4_ex <- raster::extract(tas_4, f_sp, buffer = 3000, fun = mean)
tas_4_ex <- as.data.frame(tas_4_ex)
colnames(tas_4_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
tas_4_ex$ssites <- f_sites$ssites
tas_4_ex2 <- gather(tas_4_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(tas_4_ex2)[3] <- c('tas_4_3km')

#sfc wind summer 
sfc_ex <- raster::extract(sfc_2, f_sp, buffer = 3000, fun = mean)
sfc_ex <- as.data.frame(sfc_ex)
colnames(sfc_ex) <- c("X1997", "X1998", "X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006","X2007", "X2008", "X2009", "X2010", "X2011", "X2012", "X2013", "X2014", "X2015")
sfc_ex$ssites <- f_sites$ssites
sfc_ex2 <- gather(sfc_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(sfc_ex2)[3] <- c('sfc_2_3km')

# bind prop_x columns together
w_vars_3km <- as.data.frame(cbind(precip_ex2$ssites, precip_ex2$CountYear, precip_ex2$precip_1_3km, sfc_ex2$sfc_2_3km, tas_ex2$tas_2_3km, tas_4_ex2$tas_4_3km))
colnames(w_vars_3km) <- c("ssites", "CountYear", "precip_1_3km", "sfc_2_3km",  "tas_2_3km", "tas_4_3km")
w_vars_3km$CountYear <- substr(w_vars_3km$CountYear, 2,5) # remove 'X' from year
w_vars_3km[, c(2:6)] <- lapply(w_vars_3km[, c(2:6)], as.numeric) 
w_vars_3km$ssites <- as.integer(w_vars_3km$ssites)
head(w_vars_3km)
# merge with year
field3 <- left_join(field3, w_vars_3km, by = c('ssites', 'CountYear'))
(proc.time() - t)

head(field3)

#####
# Save csv 

write.csv(field3, 'field_chpt3.csv')



#####################
## Detector info
dtr <- fread('field_detectors.csv')
head(dtr)
dtr <- dtr[, c(1,2, 4:8,11)]
str(dtr)

unique(dtr$Microphone_type)
unique(dtr$Detector_type)

field2 <- left_join(field3, dtr, by = 'Detector')
head(field2)

##

########################
## Add a dummy variable column for volunteer self assessment where they have ranked themselves as 'poor' group all the others 

field2$VolSkillSA_Poor <- ifelse(field$VolSkillSelfAssessment == 'Poor', 1, 0)
table(field2$VolSkillSA_Poor, field2$VolSkillSelfAssessment)

write.csv(field2, 'field_chpt3.csv')
# Add region to df

reg <- rgdal::readOGR('Eng_regions/Regions_December_2015_Ultra_Generalised_Clipped_Boundaries_in_England.shp')
plot(reg)

f_sp <- field
coordinates(f_sp) <- ~longitude +latitude
points(f_sp)

eng_reg <- raster::intersect(f_sp, reg)
eng_reg_df <- as.data.frame(eng_reg)

eng_reg_df$reg <- eng_reg_df$rgn15nm
colnames(eng_reg_df)
f_eng <- eng_reg_df[, c(1:66)]

f_n_eng <- field[!field$country == 'eng', ]

field <- rbind(f_eng, f_n_eng)
unique(field$reg)

field$reg <- as.factor(field$reg)
field$reg <- revalue(field$reg, c("scot" = "Scotland", "wal" = "Wales"))
write.csv(field, 'field_chpt3.csv')

###################################
colnames(field2)

fieldx <- na.omit(field2)

# check correlation between evironmental vars at 1km and survey vars 
cor(fieldx[, c(12,14,15,28:33, 40:43)], method = 'spearman')
# check correlation between evironmental vars at 3km and survey vars
cor(fieldx[, c(12,14,15,34:39, 44:47)], method = 'spearman')

################################################################
################################################################

################################################################

################################################################

################################################################
########## -- previous code to get esa proportions -- ########## 


esa_files <- list.files(path = "C:/Users/ellab/Dropbox/PhD/NBMP_data/ESA_LC_1992_2015", pattern=".tif", full.names = TRUE)
esa <- stack(esa_files)

#read in shapefile of UK in order to crop esa rasters
gb <- readOGR("UK_clip_2.shp", "UK_clip_2")
e <- extent(gb)
# crop esa to just the UK
esa_uk <- crop(esa, e)

# rename each layer to the year 
names(esa_uk) <- c('X1992', 'X1993', 'X1994', 'X1995', 'X1996', 'X1997', 'X1998', 'X1999', 'X2000', 
                   'X2001', 'X2002', 'X2003', 'X2004', 'X2005', 
                   'X2006', 'X2007', 'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015')

esa_uk_sub <- subset(esa_uk, 7:24)

esa_uk <- esa_uk_sub

# call functions for proportion of land cover classes
source('esa_prop_funcs.R')

## Extract at 1km and 3km buffers
f_sp <- f_sites
coordinates(f_sp) <- ~ longitude + latitude
proj4string(f_sp) <- proj4string(esa_uk)

# 1km 
t = proc.time()
agri_ex = raster::extract(esa_uk, f_sp, buffer=1000, fun=prop_agri)
agri_ex <- as.data.frame(agri_ex)
agri_ex$ssites = f_sites$ssites
agri_ex2 <- gather(agri_ex, key = CountYear, value = value, -ssites) # change rows into columns
head(agri_ex2)
colnames(agri_ex2)[3] <- c('prop_agri')

forest_ex = raster::extract(esa_uk, f_sp, buffer=1000, fun=prop_forest)
forest_ex <- as.data.frame(forest_ex)
forest_ex$ssites = f_sites$ssites
forest_ex2 <- gather(forest_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(forest_ex2)[3] <- c('prop_forest')

grass_ex = raster::extract(esa_uk, f_sp, buffer=1000, fun=prop_grass)
grass_ex <- as.data.frame(grass_ex)
grass_ex$ssites = f_sites$ssites
grass_ex2 <- gather(grass_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(grass_ex2)[3] <- c('prop_grass')

water_ex = raster::extract(esa_uk, f_sp, buffer=1000, fun=prop_water)
water_ex <- as.data.frame(water_ex)
water_ex$ssites = f_sites$ssites
water_ex2 <- gather(water_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(water_ex2)[3] <- c('prop_water')

urb_ex = raster::extract(esa_uk, f_sp, buffer=1000, fun=prop_urb)
urb_ex <- as.data.frame(urb_ex)
urb_ex$ssites = f_sites$ssites
urb_ex2 <- gather(urb_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(urb_ex2)[3] <- c('prop_urban')

other_ex = raster::extract(esa_uk, f_sp, buffer=1000, fun=prop_other)
other_ex <- as.data.frame(other_ex)
other_ex$ssites = f_sites$ssites
other_ex2 <- gather(other_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(other_ex2)[3] <- c('prop_other')

# bind together into one data frame
prop_1km <- as.data.frame(cbind(agri_ex2$ssites, agri_ex2$CountYear, forest_ex2$prop_forest, grass_ex2$prop_grass, agri_ex2$prop_agri, 
                                water_ex2$prop_water, urb_ex2$prop_urban, other_ex2$prop_other))
colnames(prop_1km) <- c("ssites", "CountYear", "prop_forest", "prop_grass", "prop_agri", "prop_water", "prop_urban", "prop_other")
prop_1km$CountYear <- substr(prop_1km$CountYear, 2,5) # remove 'X' from each year
prop_1km[, c(2:8)] <- lapply(prop_1km[, c(2:8)], as.numeric) 
prop_1km$ssites <- as.integer(prop_1km$ssites)
head(prop_1km)
# merge with field 
field3 <- left_join(field, prop_1km, by = c('ssites', 'CountYear'))



# 3km
t = proc.time()
agri_ex = raster::extract(esa_uk, f_sp, buffer=3000, fun=prop_agri)
agri_ex <- as.data.frame(agri_ex)
agri_ex$ssites = f_sites$ssites
agri_ex2 <- gather(agri_ex, key = CountYear, value = value, -ssites) # change rows into columns
head(agri_ex2)
colnames(agri_ex2)[3] <- c('prop_agri_3km')

forest_ex = raster::extract(esa_uk, f_sp, buffer=3000, fun=prop_forest)
forest_ex <- as.data.frame(forest_ex)
forest_ex$ssites = f_sites$ssites
forest_ex2 <- gather(forest_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(forest_ex2)[3] <- c('prop_forest_3km')

grass_ex = raster::extract(esa_uk, f_sp, buffer=3000, fun=prop_grass)
grass_ex <- as.data.frame(grass_ex)
grass_ex$ssites = f_sites$ssites
grass_ex2 <- gather(grass_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(grass_ex2)[3] <- c('prop_grass_3km')

water_ex = raster::extract(esa_uk, f_sp, buffer=3000, fun=prop_water)
water_ex <- as.data.frame(water_ex)
water_ex$ssites = f_sites$ssites
water_ex2 <- gather(water_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(water_ex2)[3] <- c('prop_water_3km')

urb_ex = raster::extract(esa_uk, f_sp, buffer=3000, fun=prop_urb)
urb_ex <- as.data.frame(urb_ex)
urb_ex$ssites = f_sites$ssites
urb_ex2 <- gather(urb_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(urb_ex2)[3] <- c('prop_urban_3km')

other_ex = raster::extract(esa_uk, f_sp, buffer=3000, fun=prop_other) 
other_ex <- as.data.frame(other_ex)
other_ex$ssites = f_sites$ssites
other_ex2 <- gather(other_ex, key = CountYear, value = value, -ssites) # change rows into columns
colnames(other_ex2)[3] <- c('prop_other_3km')

# bind prop_x columns together
prop_3km <- as.data.frame(cbind(agri_ex2$ssites, agri_ex2$CountYear, forest_ex2$prop_forest_3km, grass_ex2$prop_grass_3km, agri_ex2$prop_agri_3km, 
                                water_ex2$prop_water_3km, urb_ex2$prop_urban_3km, other_ex2$prop_other_3km))
colnames(prop_3km) <- c("ssites", "CountYear", "prop_forest_3km", "prop_grass_3km", "prop_agri_3km", "prop_water_3km", "prop_urban_3km", "prop_other_3km")
prop_3km$CountYear <- substr(prop_3km$CountYear, 2,5) # remove 'X' from year
prop_3km[, c(2:8)] <- lapply(prop_3km[, c(2:8)], as.numeric) 
prop_3km$ssites <- as.integer(prop_3km$ssites)
head(prop_3km)
# merge with year
field3 <- left_join(field3, prop_3km, by = c('ssites', 'CountYear'))
(proc.time() - t)

head(field3)
