# *****************************************************************************************
# Setting up covariate data and running JAGS with real data sets
# 
# use this with functions from
#  PrepCovDataForJags
#  RunJags_DRYModelUpdate_RunMulticov

library(raster)
library(maptools)
library(rgdal)
library("data.table")
library("lattice")
library("gplots")
library(ggplot2)  #fortify
library("plyr")  # join = tools for splitting, applying and merging data
library(sp)
require(moments)  #  alternative for skewness is the package: e1071, psych(?)
require(abind)

# getwd()
# setwd("E:/PhD2016/ResearchExamplesAndTrials/VicDataExercise")


################################################
#
# All the data
#


owl_glider_cov <- read.csv("so_occupancy_cov.csv")
owl_glider_det <- read.csv("so_detection_cov.csv")
owl_glider_sp  <- read.csv("so_species_data.csv")
View(owl_glider_sp)
View(owl_glider_det)
View(owl_glider_cov)
sites <- owl_glider_sp[,1:2]
View(sites)




################################################
#
# Occupancy variables
#

owl_glider_cov <- read.csv("so_occupancy_cov.csv")

#############################
# create a standardised set of data for these occupancy variables
# 1 JanRainFall 
# 2 JulRainFall
# 3 elevation
# 4 terrain_ruggedness_index_sept2012
# 5 wetness_index_saga_sept2012
# 6 topographic_position_index_sept2012
# 7 SouthSummerEVI
# 8 SouthSummerNDVI
ModelOccCov.1 <- owl_glider_cov[,1:8]
colnames(ModelOccCov.1)
# standardise the variables
ModelOccCov.1.std.list <- std_occ_vars(ModelOccCov.1)
ModelOccCov.1.std <- ModelOccCov.1.std.list$hab.std  # matrix: Nsite x NCov

######################
# data exploration

# wetness covariate
# will use SAGA wetness index: wetness_index_saga_sept2012
wetness <- owl_glider_cov$wetness_index_saga_sept2012
hist(wetness,main="Histogram of Saga Wetness Index",xlab="wetness index")
skewness(owl_glider_cov$wetness_index_saga_sept2012)
kurtosis(owl_glider_cov$wetness_index_saga_sept2012)

# vegetation lushness
# nvdi: Normalised Difference vegetation index  - less influenced by topography
# evi:  Enhanced Vegetation Index  -less influenced by cloud/terrain
# considered complementary
par(mfrow=c(2,2))
hist(owl_glider_cov$SouthSummerEVI ,main="Histogram of EVI", xlab="Enhanced vegetation index")
hist(owl_glider_cov$SouthSummerNDVI ,main="Histogram of NDVI", xlab="Normalised Difference vegetation index")
#hist_vegidx_2panel
plot(owl_glider_cov$SouthSummerEVI,owl_glider_cov$SouthSummerNDVI,
     main="vegetation index readings: EVI versus NVDI",
     xlab="Enhanced vegetation index",ylab="Normalised Difference vegetation index")
skewness(owl_glider_cov$SouthSummerNDVI)
skewness(owl_glider_cov$SouthSummerEVI)
kurtosis(owl_glider_cov$SouthSummerNDVI)
kurtosis(owl_glider_cov$SouthSummerEVI)

# terrain ruggedness
rugged_shape <- owl_glider_cov$topographic_position_index_sept2012
rugged_slope <- owl_glider_cov$terrain_ruggedness_index_sept2012
par(mfrow=c(2,2))
hist(rugged_slope ,main=list("Histogram of Terrain Ruggedness (Change)",cex=0.9),xlab="Terrain Ruggedness Index")
hist(rugged_slope ,main=list("Histogram of Terrain Ruggedness",cex=0.9),xlab="Log Terrain Ruggedness Index")
plot(rugged_shape,rugged_slope,
     main="Terrain ruggedness index readings: TPI versus Ruggedness",
     xlab="Terrain Ruggedness Index (Average slope)",ylab="Terrain Position Index (shape)")
skewness(rugged_shape)
skewness(rugged_slope)
kurtosis(rugged_shape)
kurtosis(rugged_slope)

# rainfall 
par(mfrow=c(2,2))
hist(owl_glider_cov$RainfallJan75m ,main=list("Histogram of Summer Rainfall (January)",cex=0.9),xlab="Rainfall (mm)")
hist(owl_glider_cov$RainfallJul75m ,main=list("Histogram of Winter Rainfall (July)",cex=0.9),xlab="Rainfall (mm)")
skewness(owl_glider_cov$RainfallJan75m)
skewness(owl_glider_cov$RainfallJul75m)
kurtosis(owl_glider_cov$RainfallJan75m)
kurtosis(owl_glider_cov$RainfallJul75m)

# elevation (possible covaraite only one speces effected)
hist(owl_glider_cov$Dem75mInteger  ,main="Histogram of Elevation",xlab="Elevation")
skewness(owl_glider_cov$Dem75mInteger)
kurtosis(owl_glider_cov$Dem75mInteger)



################################################
#
# detection variables
#
 
owl_glider_det <- read.csv("so_detection_cov.csv")
# View(owl_glider_det)
Model1DecCov.wind <- cbind(owl_glider_det$wind1,owl_glider_det$wind2,owl_glider_det$wind3)
colnames(Model1DecCov.wind) <- colnames(owl_glider_det[,1:3])
Model1DecCov.yday <- cbind(owl_glider_det$yday1,owl_glider_det$yday2,owl_glider_det$yday3)
colnames(Model1DecCov.yday) <- colnames(owl_glider_det[,4:6])
Model1DecCov.tpi <- rep(owl_glider_cov[,6],3)
dim(Model1DecCov.tpi) <- c(length(owl_glider_cov[,6]),3)
colnames(Model1DecCov.tpi) <- c("TPI1", "TPI2", "TPI3" ) 
Model1DecCov.season <- ifelse ( Model1DecCov.yday>=0, 1,0)  # 1=winter, 0=Autumn
colnames(Model1DecCov.season) <- c("Season1", "Season2", "Season3" ) 


################
# Data exploration
par(mfrow=c(2,2))
for ( i in 1:3) {
  hist(Model1DecCov.wind[,i],main=paste0("Wind Day ",i),xlab="Wind")
}
par(mfrow=c(2,2))
for ( i in 1:3) {
  hist(Model1DecCov.yday[,i],main=paste0("Histogram of yday",i),xlab="yday")
}


################################
# groupings of variables tried: 

# group 1
var_list <- list(wind=Model1DecCov.wind,yday=Model1DecCov.yday,tpi=Model1DecCov.tpi)
ModelDecCov.1.array <- array( dim=c(wind.dim[1],wind.dim[2],length(var_list)) ) # sites, reps, no_vars
for ( i in 1:length(var_list)) {
  ModelDecCov.1.array[,,i] <- var_list[[i]]
}
dim(ModelDecCov.1.array)
dimnames(ModelDecCov.1.array) <- list(NULL,c("Observation_1","Observation_2","Observation_3"),c("Wind","YDay","TPI"))
dimnames(ModelDecCov.1.array)

# group 1 with standardisation
# standardize dec vars yday and tpi, but not wind
ModelDecCov.1.std.list <-  std_det_vars(ModelDecCov.1.array[,,c(2,3)])
ModelDecCov.1.std <- ModelDecCov.1.std.list$det.std
ModelDecCov.1.all <- array( dim=dim(ModelDecCov.1.array))
ModelDecCov.1.all[,,1] <- ModelDecCov.1.array[,,1]
ModelDecCov.1.all[,,2] <- ModelDecCov.1.std[,,1]
ModelDecCov.1.all[,,3] <- ModelDecCov.1.std[,,2]
# set the same NA patter for the TPI (constant terrain vars)
ModelDecCov.1.all[is.na(ModelDecCov.1.std[,2,1]),2,3] <- NA
ModelDecCov.1.all[is.na(ModelDecCov.1.std[,3,1]),3,3] <- NA


#  group 2
var_list <- list(wind=Model1DecCov.wind,season=Model1DecCov.season)
wind.dim <- dim(Model1DecCov.wind)
ModelDecCov.2.all <- array( dim=c(wind.dim[1],wind.dim[2],length(var_list)) )
for ( i in 1:length(var_list)) {
  ModelDecCov.2.all[,,i] <- var_list[[i]]
}
dimnames(ModelDecCov.2.all) <- list(NULL,c("Observation_1","Observation_2","Observation_3"),c("Wind","Winter"))
dim(ModelDecCov.2.all)
dimnames(ModelDecCov.2.all)
ModelDecCov.2.std.list <- NULL  # this is null  because no values are standardised
# we have wind unstandardised and yday a categorical variable



# set the survey replications for each site
sreps <- apply(ModelDecCov.1.all[,,1],1, function(x) {sum(!is.na(x))} )



#############################################
#
# Species Data
#



# now looking at the y data for each site, each rep, eac species
# note: dim 1: nsite, dim 2: nrep, dim 3: nspec
owl_glider_sp  <- read.csv("so_species_data.csv")
# View(owl_glider_sp)
head(owl_glider_sp[,3:14])

# data we want is 12:23
MultiSpec.Svy  <- owl_glider_sp[,3:14]
MultiPA <- array(dim=dim((MultiSpec.Svy)))
for (i in 1: dim(MultiSpec.Svy)[2] ) {
  MultiPA[,i] <- ifelse( !is.na(MultiSpec.Svy[,i]),ifelse(MultiSpec.Svy[,i]>0,1,0), NA)
}
colnames(MultiPA) <- colnames(owl_glider_sp[,3:14])
# View(MultiPA)
dim(MultiPA)

# now put data into the format: site, rep, spec
n.reps <- 3
MultiPA.array <- array( dim=c( dim(MultiPA)[1],n.reps,(dim(MultiPA)[2]/n.reps) )  )
for (j in 1:(dim(MultiPA)[2]/n.reps) ) {
  for ( i in 1:n.reps ) {
    MultiPA.array[,i,j] <- MultiPA[,((j-1)*n.reps+i)]
  }
}
dim(MultiPA.array)
head(MultiPA.array)



#############################################
#
# Get training and validation ata set
#


# taking a random sample (no blocking)
# training sample
#round(202*2/3,0)
#n.site = dim(MultiPA)[1]

# set any directory preferences here
savedir_res <- ""
savedir <- ""
readdir <- ""


training <- sample.int(dim(MultiPA)[1],size=round(dim(MultiPA)[1]*2/3),replace=FALSE)

MultiPA.array.tr <- MultiPA.array[training,,]    # species data
MultiPA.array.va <- MultiPA.array[-training,,]
ModelOccCov.1.std.tr <- ModelOccCov.1.std[training,]  # occupancy data
ModelOccCov.1.std.va <- ModelOccCov.1.std[-training,]
ModelDecCov.1.all.tr <- ModelDecCov.1.all[training,,]  # detection data (2 possible sets)
ModelDecCov.1.all.va <- ModelDecCov.1.all[-training,,]
ModelDecCov.2.all.tr  <- ModelDecCov.2.all[training,,]
ModelDecCov.2.all.va  <- ModelDecCov.2.all[-training,,]
no_train = dim(MultiPA.array.tr)[1]
no_valid = dim(MultiPA.array.va)[1]


# positioning data (not used directly)
svy_pos.tr <- owl_glider_sp[training,1:3]
svy_pos.va <- owl_glider_sp[-training,1:3]

# put together the training and validation data
# make first list with the non corrected habitation covariate data
ModelData.1.std <- list(obs.tr=MultiPA.array.tr, obs.va=MultiPA.array.va,
                        dec.tr=ModelDecCov.1.all.tr, dec.va=ModelDecCov.1.all.va,
                        hab.tr =ModelOccCov.1.std.tr, std.hab.va=ModelOccCov.1.std.va,
                        svy_pos.tr=svy_pos.tr, svy_pos.va=svy_pos.va,
                        cov_pos.tr=svy_pos.tr,cov_pos.va=svy_pos.va,
                        ModelDecCov.1.std.list,ModelOccCov.1.std.list   # need lists to get the mean and sd data
)
ModelData.2.std <- list(obs.tr=MultiPA.array.tr, obs.va=MultiPA.array.va,
                        dec.tr=ModelDecCov.2.all.tr, dec.va=ModelDecCov.2.all.va,
                        hab.tr =ModelOccCov.1.std.tr, std.hab.va=ModelOccCov.1.std.va,
                        svy_pos.tr=svy_pos.tr, svy_pos.va=svy_pos.va,
                        cov_pos.tr=svy_pos.tr,cov_pos.va=svy_pos.va,
                        ModelDecCov.2.std.list,ModelOccCov.1.std.list   # need lists to get the mean and sd data
)




# Full data lists for detection groups 1 and 2: training and validation groups
realdat.1.std.train1 <- list(obs.tr=ModelData.1.std$obs.tr,obs.va=ModelData.1.std$obs.va,
                             hab.tr=ModelData.1.std$hab.tr,hab.va=ModelData.1.std$std.hab.va,
                             det.tr=ModelData.1.std$dec.tr,det.va=ModelData.1.std$dec.va,
                             n.train = dim(ModelData.1.std$hab.tr)[1],
                             n.valid = dim(ModelData.1.std$std.hab.va)[1],
                             n.spec = dim(ModelData.1.std$obs.tr)[3],
                             n.rep = dim(ModelData.1.std$obs.tr)[2]
                             )
realdat.2.std.train1 <- list(obs.tr=ModelData.2.std$obs.tr,obs.va=ModelData.2.std$obs.va,
                             hab.tr=ModelData.2.std$hab.tr,hab.va=ModelData.2.std$std.hab.va,
                             det.tr=ModelData.2.std$dec.tr,det.va=ModelData.2.std$dec.va,
                             n.train = dim(ModelData.2.std$hab.tr)[1],
                             n.valid = dim(ModelData.2.std$std.hab.va)[1],
                             n.spec = dim(ModelData.2.std$obs.tr)[3],
                             n.rep = dim(ModelData.2.std$obs.tr)[2]
)


# for all 4 species
windats.1.std.train1.4sp <- set.win.cov.real(realdat.1.std.train1)
windats.2.std.train1.4sp <- set.win.cov.real(realdat.2.std.train1)
saveRDS(windats.1.std.train1.4sp,paste0(savedir,"windats.1.std.train1.4sp.RDS"))
saveRDS(windats.2.std.train1.4sp,paste0(savedir,"windats.2.std.train1.4sp.RDS"))



#####################################
# example of cut down data, can do this for any data set

# species combinations:
PO_v_SO = c(1,2)
PO_v_GG = c(1,3)
PO_v_YB = c(1,4)
SO_v_GG = c(2,3)
SO_v_YB = c(2,4)
GG_v_YB = c(3,4)
PO_v_SO_GG = c(1,2,3)
PO_v_SO_YB = c(1,2,4)
PO_v_GG_YB = c(1,3,4)
SO_v_GG_YB = c(2,3,4)

# illustrating two subsets with 2 species
# for 2 species: SO versus GG, correlation, preditor/prey
SO_v_GG = c(2,3)
realdat.1.std.train1.SO_GG <- list(obs.tr=ModelData.1.std$obs.tr[,,SO_v_GG],obs.va=ModelData.1.std$obs.va[,,SO_v_GG],
                                   hab.tr=ModelData.1.std$hab.tr,hab.va=ModelData.1.std$std.hab.va,
                                   det.tr=ModelData.1.std$dec.tr,det.va=ModelData.1.std$dec.va,
                                   n.train = dim(ModelData.1.std$hab.tr)[1],
                                   n.valid = dim(ModelData.1.std$std.hab.va)[1],
                                   n.spec = dim(ModelData.1.std$obs.tr)[3],
                                   n.rep = dim(ModelData.1.std$obs.tr)[2]
)
windats.1.std.train1.SO_GG <- set.win.cov.real(realdat.1.std.train1.SO_GG)
saveRDS(realdat.1.std.train1.SO_GG,paste0(savedir,"realdat.1.std.train1.SO_GG",".RDS"))
saveRDS(windats.1.std.train1.SO_GG,paste0(savedir,"windats.1.std.train1.SO_GG",".RDS"))


realdat.2.std.train1.SO_GG <- list(obs.tr=ModelData.2.std$obs.tr[,,SO_v_GG],obs.va=ModelData.2.std$obs.va[,,SO_v_GG],
                                   hab.tr=ModelData.2.std$hab.tr,hab.va=ModelData.2.std$std.hab.va,
                                   det.tr=ModelData.2.std$dec.tr,det.va=ModelData.2.std$dec.va,
                                   n.train = dim(ModelData.2.std$hab.tr)[1],
                                   n.valid = dim(ModelData.2.std$std.hab.va)[1],
                                   n.spec = dim(ModelData.2.std$obs.tr[,,SO_v_GG])[3],
                                   n.rep = dim(ModelData.2.std$obs.tr)[2]
)
windats.2.std.train1.SO_GG <- set.win.cov.real(realdat.2.std.train1.SO_GG)
saveRDS(realdat.2.std.train1.SO_GG,paste0(savedir,"realdat.2.std.train1.SO_GG",".RDS"))
saveRDS(windats.2.std.train1.SO_GG,paste0(savedir,"windats.2.std.train1.SO_GG",".RDS"))


# for 2 species: PO versus GG, correlation, preditor/prey
PO_v_GG = c(1,3)
realdat.1.std.train1.PO_GG <- list(obs.tr=ModelData.1.std$obs.tr[,,PO_v_GG],obs.va=ModelData.1.std$obs.va[,,PO_v_GG],
                             hab.tr=ModelData.1.std$hab.tr,hab.va=ModelData.1.std$std.hab.va,
                             det.tr=ModelData.1.std$dec.tr,det.va=ModelData.1.std$dec.va,
                             n.train = dim(ModelData.1.std$hab.tr)[1],
                             n.valid = dim(ModelData.1.std$std.hab.va)[1],
                             n.spec = dim(ModelData.1.std$obs.tr)[3],
                             n.rep = dim(ModelData.1.std$obs.tr)[2]
)
windats.1.std.train1.PO_GG <- set.win.cov.real(realdat.1.std.train1.PO_GG)
saveRDS(realdat.1.std.train1.PO_GG,paste0(savedir,"realdat.1.std.train1.PO_GG",".RDS"))
saveRDS(windats.1.std.train1.PO_GG,paste0(savedir,"windats.1.std.train1.PO_GG",".RDS"))


realdat.2.std.train1.PO_GG <- list(obs.tr=ModelData.2.std$obs.tr[,,PO_v_GG],obs.va=ModelData.2.std$obs.va[,,PO_v_GG],
                                   hab.tr=ModelData.2.std$hab.tr,hab.va=ModelData.2.std$std.hab.va,
                                   det.tr=ModelData.2.std$dec.tr,det.va=ModelData.2.std$dec.va,
                                   n.train = dim(ModelData.2.std$hab.tr)[1],
                                   n.valid = dim(ModelData.2.std$std.hab.va)[1],
                                   n.spec = dim(ModelData.2.std$obs.tr)[3],
                                   n.rep = dim(ModelData.2.std$obs.tr)[2]
)

windats.2.std.train1.PO_GG <- set.win.cov.real(realdat.2.std.train1.PO_GG)
saveRDS(realdat.2.std.train1.PO_GG,paste0(savedir,"realdat.2.std.train1.PO_GG",".RDS"))
saveRDS(windats.2.std.train1.PO_GG,paste0(savedir,"windats.2.std.train1.PO_GG",".RDS"))



#############################################
# list of data

# data sets
std_list.tr1 <- list(windats.1.std.train1.4sp,
                 windats.1.std.train1.SO_GG,
                 windats.1.std.train1.PO_GG,
                 windats.1.std.train1.PO_YB)
# id list
id.list.tr1 <- c("All_4sp","SO_GG","PO_GG")
std_list_id.tr1 <- windats_id.list(std_list.tr1,id.list.tr1)

saveRDS(std_list.tr1,paste0(savedir,"std_list.tr1.RDS"))
saveRDS(std_list_id.tr1,paste0(savedir,"std_list_id.tr1.RDS"))

# retrieve
std_list.tr1 <- readRDS(paste0(savedir,"std_list.tr1.RDS"))
std_list_id.tr1 <- readRDS(paste0(savedir,"std_list_id.tr1.RDS"))




####################################################
# reduced variables
#
# many variables were not signficant after analysis for
# any species.  This set of covariates resulted in more significant results


# data with reduced data set:
# reduce
dim(ModelData.1.std$hab.tr)
# reduced occ var set
hab_red.set = c(1,3,5,6,8)   # summer rain, elevation, wetness (saga index), TPI, NDVI (south summer)
det_red.set = c(1,2)   # wind, season (winter or autumn)
dim(ModelData.1.std$hab.tr[,hab_red.set])

realdat.3.std.train1 <- list(obs.tr=ModelData.2.std$obs.tr,obs.va=ModelData.2.std$obs.va,
                             hab.tr=ModelData.2.std$hab.tr[,hab_red.set],hab.va=ModelData.2.std$std.hab.va[,hab_red.set],
                             det.tr=ModelData.2.std$dec.tr[,,det_red.set],det.va=ModelData.2.std$dec.va[,,det_red.set],
                             n.train = dim(ModelData.2.std$hab.tr)[1],
                             n.valid = dim(ModelData.2.std$std.hab.va)[1],
                             n.spec = dim(ModelData.2.std$obs.tr)[3],
                             n.rep = dim(ModelData.2.std$obs.tr)[2]
)

# for all 4 species
windats.3.std.train1.4sp <- set.win.cov.real(realdat.3.std.train1)
saveRDS(realdat.3.std.train1,paste0(savedir,"realdat.3.std.train1",".RDS"))
saveRDS(windats.3.std.train1.4sp,paste0(savedir,"windats.3.std.train1.4sp",".RDS"))

# for 2 species: SO versus GG, correlation, preditor/prey
SO_v_GG = c(2,3)
realdat.3.std.train1.SO_GG <- list(obs.tr=ModelData.2.std$obs.tr[,,SO_v_GG],obs.va=ModelData.2.std$obs.va[,,SO_v_GG],
                                   hab.tr=ModelData.2.std$hab.tr[,hab_red.set],hab.va=ModelData.2.std$std.hab.va[,hab_red.set],
                                   det.tr=ModelData.2.std$dec.tr[,,det_red.set],det.va=ModelData.2.std$dec.va[,,det_red.set],
                                   n.train = dim(ModelData.2.std$hab.tr)[1],
                                   n.valid = dim(ModelData.2.std$std.hab.va)[1],
                                   n.spec = dim(ModelData.2.std$obs.tr)[3],
                                   n.rep = dim(ModelData.2.std$obs.tr)[2]
)
windats.3.std.train1.SO_GG <- set.win.cov.real(realdat.3.std.train1.SO_GG)
saveRDS(realdat.3.std.train1.SO_GG,paste0(savedir,"realdat.3.std.train1.SO_GG",".RDS"))
saveRDS(windats.3.std.train1.SO_GG,paste0(savedir,"windats.3.std.train1.SO_GG",".RDS"))


PO_v_GG = c(1,3)
realdat.3.std.train1.PO_GG <- list(obs.tr=ModelData.2.std$obs.tr[,,PO_v_GG],obs.va=ModelData.2.std$obs.va[,,PO_v_GG],
                                   hab.tr=ModelData.2.std$hab.tr[,hab_red.set],hab.va=ModelData.2.std$std.hab.va[,hab_red.set],
                                   det.tr=ModelData.2.std$dec.tr[,,det_red.set],det.va=ModelData.2.std$dec.va[,,det_red.set],
                                   n.train = dim(ModelData.2.std$hab.tr)[1],
                                   n.valid = dim(ModelData.2.std$std.hab.va)[1],
                                   n.spec = dim(ModelData.2.std$obs.tr)[3],
                                   n.rep = dim(ModelData.2.std$obs.tr)[2]
)
windats.3.std.train1.PO_GG <- set.win.cov.real(realdat.3.std.train1.PO_GG)
saveRDS(realdat.3.std.train1.PO_GG,paste0(savedir,"realdat.3.std.train1.PO_GG",".RDS"))
saveRDS(windats.3.std.train1.PO_GG,paste0(savedir,"windats.3.std.train1.PO_GG",".RDS"))


# data lists
std_list.3.tr1 <- list(windats.3.std.train1.4sp,
                       windats.3.std.train1.SO_GG,
                       windats.3.std.train1.PO_GG)

# id list
id.list.tr1 <- c("All_4sp","SO_GG","PO_GG")
std_list_id.3.tr1 <- windats_id.list(std_list.3.tr1,id.list.tr1)
saveRDS(std_list.3.tr1,paste0(savedir,"std_list.3.tr1",".RDS"))
saveRDS(std_list_id.3.tr1,paste0(savedir,"std_list_id.3.tr1",".RDS"))

# retrieve list later
std_list.3.tr1 <- readRDS(paste0(savedir,"std_list.3.tr1",".RDS"))
std_list_id.3.tr1 <- readRDS(paste0(savedir,"std_list_id.3.tr1",".RDS"))



#################################################################
#
# Run the MVP Jags Model
#
###################################################################


# model with detection
# model with collapsed y
# model no detection, from 1st ob data (as missing data in others)

# model types to run 
# det models used with data with replications and explicit detection
model_name.dep.det <- "model11_10_sehII_dep_mvp_cov_v05b.txt"
model_name.dep.det.Kcov <- "model11_10_sehII_dep_mvp_Kcov_v05b.txt"
model_name.dep.det.Kcov.rl <- "model11_10_sehII_dep_mvp_Kcov_v06b.txt"

# nodet models used with collapsed data or single observation data
model_name.dep.nodet <- "model11_10_sehII_dep_nodet_mvp_cov_v05c.txt"
model_name.dep.nodet.Kcov <- "model11_10_sehII_dep_nodet_mvp_Kcov_v05c.txt"

# set the JAGS initialisation and burn in parameters
# params.dep.det  # from RunJAgs_DRY_Models_seh_idep_v07_0.  .ch.R
# quick check that functions are running
ni=800;nt=10;nb=400;nc=3;na=200


# recommended
ni=150000;nt=10;nb=100000;nc=3;na=1000

# you must set a fid string.  This could be a meaningful string or a random string
# to identify different runs
fid=""
fid="test1"

#**********************************************************************
# follow up JAGS studies 
#
# Concentrating on SO and GG which have best correlation
# want to eliminate unused/insignificant parameters to improve estimates
# start with set



#**********************************************************************
# Now run jags - study 2
#
# run for all variables - using standardisation only

# 
# windats: data for the MVP Jags model in correct format:
# windats.3.std.train1.4sp
# windats.3.std.train1.SO_GG
# windats.3.std.train1.PO_GG

# create a list of data sets to examine
# note every data set in the list will run
windat.list <- list()
windat.list[[1]] <-windats.3.std.train1.4sp
windat.list[[2]] <-windats.3.std.train1.SO_GG
windat.list[[3]] <-windats.3.std.train1.PO_GG

# if want a single data set to run
windat.list <- list()
windat.list[[1]] <-windats.3.std.train1.4sp



###################################
#
# running the model with explicit detection and multiple replications

real_std.stdy3.full.det <- Run_jags.clist.cov.rl(caseList=windat.list,type=1,wnd.idx=1,
                                          params=params.dep.det,model_name=model_name.dep.det.Kcov.rl,
                                          ni=ni,nt=nt,nb=nb,nc=nc,na=na,
                                          fid=fid,savedir=savedir_res)

real_std.stdy3.full.det[[1]]$resFile <- paste0(real_std.stdy3.full.det[[1]]$resFile,".RDS")
saveRDS(real_std.stdy3.full.det,paste0(savedir_res,"real_std.stdy3.full.det",fid,".RDS"))



###################################
#
# running the model with no explicit detection and collapsed data

# collpased Y
real_std.stdy3.full.clpy <- Run_jags.clist.cov.rl(caseList=windat.list,type=2,wnd.idx=1,
                                                 params=params.dep.nodet,model_name=model_name.dep.nodet.Kcov,
                                                 ni=ni,nt=nt,nb=nb,nc=nc,na=na,
                                                 fid=fid,savedir=savedir_res)

real_std.stdy3.full.clpy[[1]]$resFile <- paste0(real_std.stdy3.full.clpy[[1]]$resFile,".RDS")
saveRDS(real_std.stdy3.full.clpy,paste0(savedir_res,"real_std.stdy3.full.clpy",fid,".RDS"))




###################################
#
# running the model with no explicit detection and a single observation

# first obs only
real_std.stdy3.full.fst <- Run_jags.clist.cov.rl(caseList=windat.list,type=3,wnd.idx=1,
                                                  params=params.dep.nodet,model_name=model_name.dep.nodet.Kcov,
                                                  ni=ni,nt=nt,nb=nb,nc=nc,na=na,
                                                  fid=fid,savedir=savedir_res)

real_std.stdy3.full.fst[[1]]$resFile <- paste0(real_std.stdy3.full.fst[[1]]$resFile,".RDS")
saveRDS(real_std.stdy3.full.fst,paste0(savedir_res,"real_std.stdy3.full.fst",fid,".RDS"))

#*************************************
# get the print out variables



pr_real_std.stdy3.full <- list()
# put data into easily accessible format
pr_real_std.stdy3.full[[1]] <- printouts_rl2(list_jagfiles=real_std.stdy3.full.det,windat.list=windat.list,n.train=realdat.2.std.train1.SO_GG$n.train,wnd.idx=1,savedir=savedir_res)
pr_real_std.stdy3.full[[2]] <- printouts_rl2(list_jagfiles=real_std.stdy3.full.clpy,windat.list=windat.list,n.train=realdat.2.std.train1.SO_GG$n.train,wnd.idx=1,savedir=savedir_res)
pr_real_std.stdy3.full[[3]] <- printouts_rl2(list_jagfiles=real_std.stdy3.full.fst,windat.list=windat.list,n.train=realdat.2.std.train1.SO_GG$n.train,wnd.idx=1,savedir=savedir_res)
saveRDS(pr_real_std.stdy3.full,paste0(savedir_res,"pr_real_std.stdy3.full",fid,".RDS"))

# Model with detection
str(pr_real_std.stdy3.full[[1]][[1]])
pr_real_std.stdy3.full[[1]][[1]]$beta_sum
pr_real_std.stdy3.full[[1]][[1]]$alpha_sum
pr_real_std.stdy3.full[[1]][[1]]$rho.est[1,2,]
pr_real_std.stdy3.full[[1]][[1]]$DIC

# collapsed data
pr_real_std.stdy3.full[[2]][[1]]$beta_sum
pr_real_std.stdy3.full[[2]][[1]]$alpha_sum
pr_real_std.stdy3.full[[2]][[1]]$rho.est[1,2,]
pr_real_std.stdy3.full[[2]][[1]]$DIC  


#*************************************

# getting ROC/AUC
# this wil return ROC values but also
# plot graphs of ROC/AIC

# getOccDecProb from BasicsimJagsPrepV05b
ROC_AUC_std.stdy3.list <- list()


idx=1;type=1 
ROC_AUC_std.stdy3.list[[type]] <- getOccDecProb(
              list_jagfiles=real_std.stdy3.full.det,
              list_prresp=pr_real_std.stdy3.full[[type]],
              Valid_set=windat.list[[1]][[idx]]$other$y_valid,
              savedir=savedir_res)
  
idx=1;type=2 
ROC_AUC_std.stdy3.list[[type]] <- getOccDecProb(
  list_jagfiles=real_std.stdy3.full.clpy,
  list_prresp=pr_real_std.stdy3.full[[type]],
  Valid_set=windat.list[[1]][[idx]]$other$y_valid,
  savedir=savedir_res)


#**********************************************************************************
# Process data for plotting using jagcovrespsforplot
# jagcovrespsforplot returns the correlation data
# call: jagcovrespsforplot(list_jagfiles,traceplot=FALSE,savedircov,testFlag=FALSE)

# if just interested in the correlation data and raw intercept estimates
# note correlation values calculated over all bayesian posterior estimates 
# and for all bayesian posterior parameter sets
pl_real_std.stdy3_det_corr <- jagrespsforplot(real_std.stdy3.full.det,savedircov=savedir,traceplot=FALSE,testFlag=TRUE,test_params=1:450)
pl_real_std.stdy3_clp_corr <- jagrespsforplot(real_std.stdy3.full.clpy,savedircov=savedir,traceplot=FALSE,testFlag=TRUE,test_params=1:450)

# jagcovrespsforplot as above but returns correlation estimates, plus raw estimates
# for all alpha and beta parameters based on summary data
# this is more suitable if have covariate data
pl_real_std.stdy3.full <- list()
pl_real_std.stdy3.full[[1]] <- jagcovrespsforplot(real_std.stdy3.full.det,savedircov=savedir,traceplot=FALSE,testFlag=TRUE,test_params=1:50)
saveRDS(pl_real_std.stdy3.full[[1]],paste0(savedir,"pl_","real_std.stdy3.full.det",".RDS"))

pl_real_std.stdy3.full[[2]] <- jagcovrespsforplot(real_std.stdy3.full.clpy,savedircov=savedir,traceplot=FALSE,testFlag=FALSE)
saveRDS(pl_real_std.stdy3.full[[2]],paste0(savedir,"pl_","real_std.stdy3.full.clpy",".RDS"))

pl_real_std.stdy3.full[[3]] <- jagcovrespsforplot(real_std.stdy3.full.fst,savedircov=savedir,traceplot=FALSE,testFlag=FALSE)
saveRDS(pl_real_std.stdy3.full[[2]],paste0(savedir,"pl_","real_std.stdy3.full.fst",".RDS"))



# Process data to get Beta estimates using jagcovresps_beta_est
# jagcovresps_beta_est returns standardised estimates calculated
# from the full distribution of bayesian posterior estimates
# and for all bayesian posterior parameter sets
# call: jagcovresps_beta_est(list_jagfiles, savedircov)
pl_beta_real_std.stdy3.full <- list()
pl_beta_real_std.stdy3.full[[1]] <- jagcovresps_beta_est(real_std.stdy3.full.det,savedircov=savedir)
saveRDS(pl_beta_real_std.stdy3.full[[1]],paste0(savedir,"pl_","beta_real_std.stdy3.full.det",".RDS"))

pl_beta_real_std.stdy3.full[[2]] <- jagcovresps_beta_est(real_std.stdy3.full.clpy,savedircov=savedir)
saveRDS(pl_beta_real_std.stdy3.full[[2]],paste0(savedir,"pl_","beta_real_std.stdy3.full.clpy",".RDS"))

pl_beta_real_std.stdy3.full[[3]] <- jagcovresps_beta_est(real_std.stdy3.full.fst,savedircov=savedir)
saveRDS(pl_beta_real_std.stdy3.full[[2]],paste0(savedir,"pl_","beta_real_std.stdy3.full.fst",".RDS"))


# Combine beta occupancy parameters and alpha detection parameters
# for plots that show alpha and beta parameters

pl_alphabeta_real_std.stdy3.full <- combineAlphaBetaJagsParamEstLists(pl_beta_real_std.stdy3.full,pl_real_std.stdy3.full)
pl_alphabeta_real_std.stdy3.full.det <- combineAlphaBetaJagsParamEsts(pl_beta_real_std.stdy3.full[[1]],pl_real_std.stdy3.full[[1]])


#******************************************************************
# species 1: PO
# species 2: SO
# species 3: GG
# species 4: YB
#

pr_resp <- pr_real_std.stdy3.full[[1]] # det
pr_resp <- pr_real_std.stdy3.full[[2]] # clpy
pr_resp <- pr_real_std.stdy3.full[[3]] # fst

par(mfrow=c(2,3))

title_cex=1.5
axis_cex=1.3
hist(pr_resp[[1]]$rho.list[,1,2],main=list("Correlation between PO and SO",cex=title_cex),xlab= list("Correlation Estimates",cex=axis_cex),ylab= list("Frequency",cex=axis_cex) )
abline(v=mean(pr_resp[[1]]$rho.list[,1,2]),col="red",lty="solid" )
abline(v=median(pr_resp[[1]]$rho.list[,1,2]),col="blue",lty="dashed" )
abline(v=quantile(pr_resp[[1]]$rho.list[,1,2],0.025),col="grey",lty="dotted" )
abline(v=quantile(pr_resp[[1]]$rho.list[,1,2],0.975),col="grey",lty="dotted" )


hist(pr_resp[[1]]$rho.list[,1,3],main=list("Correlation between PO and GG",cex=title_cex),xlab= list("Correlation Estimates",cex=axis_cex),ylab= list("Frequency",cex=axis_cex) )
abline(v=mean(pr_resp[[1]]$rho.list[,1,3]),col="red",lty="solid" )
abline(v=median(pr_resp[[1]]$rho.list[,1,3]),col="blue",lty="dashed" )
abline(v=quantile(pr_resp[[1]]$rho.list[,1,3],0.025),col="grey",lty="dotted" )
abline(v=quantile(pr_resp[[1]]$rho.list[,1,3],0.975),col="grey",lty="dotted" )

hist(pr_resp[[1]]$rho.list[,1,4],main=list("Correlation between PO and YB",cex=title_cex),xlab= list("Correlation Estimates",cex=axis_cex),ylab= list("Frequency",cex=axis_cex) )
abline(v=mean(pr_resp[[1]]$rho.list[,1,4]),col="red",lty="solid" )
abline(v=median(pr_resp[[1]]$rho.list[,1,4]),col="blue",lty="dashed" )
abline(v=quantile(pr_resp[[1]]$rho.list[,1,4],0.025),col="grey",lty="dotted" )
abline(v=quantile(pr_resp[[1]]$rho.list[,1,4],0.975),col="grey",lty="dotted" )

hist(pr_resp[[1]]$rho.list[,2,3],main=list("Correlation between SO and GG",cex=title_cex),xlab= list("Correlation Estimates",cex=axis_cex),ylab= list("Frequency",cex=axis_cex) )
abline(v=mean(pr_resp[[1]]$rho.list[,2,3]),col="red",lty="solid" )
abline(v=median(pr_resp[[1]]$rho.list[,2,3]),col="blue",lty="dashed" )
abline(v=quantile(pr_resp[[1]]$rho.list[,2,3],0.025),col="grey",lty="dotted" )
abline(v=quantile(pr_resp[[1]]$rho.list[,2,3],0.975),col="grey",lty="dotted" )

hist(pr_resp[[1]]$rho.list[,2,4],main=list("Correlation between SO and YB",cex=title_cex),xlab= list("Correlation Estimates",cex=axis_cex),ylab= list("Frequency",cex=axis_cex) )
abline(v=mean(pr_resp[[1]]$rho.list[,2,4]),col="red",lty="solid" )
abline(v=median(pr_resp[[1]]$rho.list[,2,4]),col="blue",lty="dashed" )
abline(v=quantile(pr_resp[[1]]$rho.list[,2,4],0.025),col="grey",lty="dotted" )
abline(v=quantile(pr_resp[[1]]$rho.list[,2,4],0.975),col="grey",lty="dotted" )

hist(pr_resp[[1]]$rho.list[,3,4],main=list("Correlation between GG and YB",cex=title_cex),xlab= list("Correlation Estimates",cex=axis_cex),ylab= list("Frequency",cex=axis_cex) )
abline(v=mean(pr_resp[[1]]$rho.list[,3,4]),col="red",lty="solid" )
abline(v=median(pr_resp[[1]]$rho.list[,3,4]),col="blue",lty="dashed" )
abline(v=quantile(pr_resp[[1]]$rho.list[,3,4],0.025),col="grey",lty="dotted" )
abline(v=quantile(pr_resp[[1]]$rho.list[,3,4],0.975),col="grey",lty="dotted" )



#******************************************************************


