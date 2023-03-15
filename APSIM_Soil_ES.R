
#Housekeeping

library(GGally)
library(apsimx)
library(ggplot2)
library(collapse)
library(tidyverse)
library(terra)
library(nasapower)
library(svMisc)
library(stringr)
library(corrplot)
library(network)
library(sna)
library(scales)
library(ggpubr)
library(gridExtra)
library(rstatix)

setwd("C:/Users/hugow/Documents/ESModelling")

##########################################
#-------------DATA COLLECTION------------#
##########################################

Model_location <- "C:/Users/hugow/Documents/ESModelling/"

#Load UK wheat raster and reduce in size
wheat_rast <- rast("Soil data/Land use 1km 2015 RASTER/ww_uk_eucropmap.tif") %>%
  project("epsg:27700") %>%
  clamp(lower= 1000, value=FALSE)


#Get coords for later use
wheat_crds <- wheat_rast %>% 
  as.points(na.rm = FALSE) %>% 
  project("epsg:4326") %>% 
  crds() %>% as.data.frame()

#Load soil data and rasterise
phbd_vect <- vect("Soil data/Bulk Density/5dd624a9-55c9-4cc0-b366-d335991073c7/data/CS_topsoil_pH_bulkDensity.shp")
CN_vect <- vect("Soil data/Topsoil C-N/7055965b-7fe5-442b-902d-63193cbe001c/data/7055965b-7fe5-442b-902d-63193cbe001c/CS_topsoil_nutrients.shp")%>%
  project("epsg:27700")
C_vect <- vect("Soil data/Topsoil C concentration/9e4451f8-23d3-40dc-9302-73e30ad3dd76/data/9e4451f8-23d3-40dc-9302-73e30ad3dd76/CS_topsoil_carbon.shp")%>%
  project("epsg:27700")
sat_rast <- rast("Soil data/FS_water_2016_07/ths_fao_octop.tif")%>%
  project("epsg:27700")
dul_rast <- rast("Soil data/FS_water_2016_07/fc_fao.tif")%>%
  project("epsg:27700")
ll15_rast <- rast("Soil data/FS_water_2016_07/wp_fao.tif")%>%
  project("epsg:27700")


#rasterise to wheat

ph_rast <- rasterize(phbd_vect, wheat_rast, field = "PH_07")
bd_rast <- rasterize(phbd_vect, wheat_rast, field  = "BULKD_07")
CN_rast <- rasterize(CN_vect, wheat_rast, field = "CN_07")
N_rast <- rasterize(CN_vect, wheat_rast, field = "NCONC_07")
C_rast <- rasterize(C_vect, wheat_rast, field = "CCONC_07")
sat_rast <- crop(sat_rast[[1]], wheat_rast)%>%
  resample(wheat_rast)
dul_rast <- crop(dul_rast[[1]], wheat_rast)%>%
  resample(wheat_rast)
ll15_rast <- crop(ll15_rast[[1]], wheat_rast)%>%
  resample(wheat_rast)

#Extract and and add to dataframe

Soildata <- spatSample(wheat_rast, 1300, method = "weights",
                   replace = FALSE, na.rm = TRUE,
                   cells = TRUE, xy = TRUE, values = TRUE)

Soildata$BD <- unlist(bd_rast[Soildata$cell])
Soildata$CN <- unlist(CN_rast[Soildata$cell])
Soildata$C <- unlist(C_rast[Soildata$cell])
Soildata$pH <- unlist(ph_rast[Soildata$cell])
Soildata$sat <- unlist(sat_rast[Soildata$cell])
Soildata$dul <- unlist(dul_rast[Soildata$cell])
Soildata$ll15 <- unlist(ll15_rast[Soildata$cell])
Soildata$Lon <- wheat_crds$x[Soildata$cell]
Soildata$Lat <- wheat_crds$y[Soildata$cell]

Soildata <- Soildata[Soildata$C >= 0, ]
Soildata <- Soildata %>% na.omit() %>% sample_n(1100, replace = FALSE)

write.csv(Soildata, file = "Outputs/Soildata_noES.csv", row.names = FALSE)

##########################################
#---------------PROCESSING---------------#
##########################################

Soildata <- read.csv("Outputs/Soildata_noES.csv")

apsim_file <- "Ecosystem_Services_Conventional.apsimx"

for(i in 1:nrow(Soildata)){
  

#Make Soil Profile with data
  
Soil_profile <- apsimx_soil_profile(
                nlayers = 4,
                Thickness = c(150, 250, 300, 300),
                BD = c(Soildata$BD[i], Soildata$BD[i], Soildata$BD[i], Soildata$BD[i]),
                AirDry = c(Soildata$ll15[i] * 0.35, Soildata$ll15[i] * 0.80, Soildata$ll15[i] * 0.80, Soildata$ll15[i] * 0.80),
                LL15 = c(Soildata$ll15[i], Soildata$ll15[i], Soildata$ll15[i], Soildata$ll15[i]),
                DUL = c(min(Soildata$dul[i], 0.552), min(Soildata$dul[i], 0.53), min(Soildata$dul[i], 0.53), min(Soildata$dul[i], 0.53)),
                SAT = c(min(Soildata$sat[i], 0.562), min(Soildata$sat[i], 0.54), min(Soildata$sat[i], 0.54), min(Soildata$sat[i], 0.54)),
                KS = c(20, 20, 20, 20),
                crop.LL = c(Soildata$ll15[i], Soildata$ll15[i], Soildata$ll15[i], Soildata$ll15[i]),
                crop.KL = c(0.06, 0.06, 0.06, 0.06),
                crop.XF = c(1, 1, 1, 1),
                Carbon = c(Soildata$C[i] * 0.1, (Soildata$C[i] * 0.09), Soildata$C[i] * 0.06, Soildata$C[i] * 0.01),
                SoilCNRatio = c(Soildata$CN[i], Soildata$CN[i], Soildata$CN[i], Soildata$CN[i]),
                FOM = c(350, 200, 1, 1),
                FOM.CN = 35,
                FBiom = c(0.05, 0.01, 0.01, 0.01),
                FInert = c(0.4, 0.6, 0.8, 1),
                NO3N = c(20, 5, 0.1, 0.1),
                NH4N = c(0.5, 0.1, 0.1, 0.1),
                PH = c(Soildata$pH[i], Soildata$pH[i], Soildata$pH[i],  Soildata$pH[i]),
                soil.bottom = 1000,
                water.table = 2000,
                soil.type = 0,
                crops = c("Wheat", "Chickpea", "WhiteClover", "AGPRyegrass"),
                metadata = NULL,
                soilwat = NA,
                swim = NA,
                soilorganicmatter = NA,
                dist.parms = list(a = 0, b = 0.2)
                )

edit_apsimx_replace_soil_profile(
  file = apsim_file,
  src.dir = Model_location,
  wrt.dir = Model_location,
  soil.profile = Soil_profile,
  edit.tag = "",
  overwrite = TRUE,
  verbose = TRUE,
)

#check soil profile has been replaced

inspect_apsimx(apsim_file, Model_location,
                node = "Soil",
                soil.child = "Physical"
               )

#Generate met file for grid location and save

cat("Getting weather data for", i)
get_power_apsim_met(
  c(Soildata$Lon[i], Soildata$Lat[i]),
  c("1995-09-01", "2020-12-31"),
  wrt.dir = "C:/Users/hugow/Documents/ESModelling/Met data/NasaPower met data/",
  filename = "ES_weather.met"
)

#Edit apsim model for new met file

edit_apsimx(apsim_file,
           node = "Weather", 
           value = "Met data/NasaPower met data/ES_weather.met",
           edit.tag = ""
           )

#Run apsim file and tidy data

cat("Running model for", i)

apsim_output_25yr <- apsimx(apsim_file)
apsim_output_25yr$N_Uptake <- apsim_output_25yr$"N_Uptake(1)" + apsim_output_25yr$"N_Uptake(2)" + apsim_output_25yr$"N_Uptake(3)" + apsim_output_25yr$"N_Uptake(4)"

apsim_output_12yr <- apsim_output_25yr[1:4383, ]
apsim_output_4yr <- apsim_output_25yr[1:1461, ]

##################################
#Nitrogen Provided by soil to wheat crop (ratio)
##################################

#25yr
N_output <- apsim_output_25yr[complete.cases(apsim_output_25yr[ ,"N_Uptake"]) & apsim_output_25yr$N_Uptake != 0, ]
N_Supply_day <- pmax(0, pmin(1, (N_output$Nmin_som + N_output$Nmin_res) / (N_output$N_Uptake * 10)))
N_Supply_avg <- mean(N_Supply_day)

Soildata[i, "NSupply25yr"] <- N_Supply_avg

#12yr
N_output <- apsim_output_12yr[complete.cases(apsim_output_12yr[ ,"N_Uptake"]) & apsim_output_12yr$N_Uptake != 0, ]
N_Supply_day <- pmax(0, pmin(1, (N_output$Nmin_som + N_output$Nmin_res) / (N_output$N_Uptake * 10)))
N_Supply_avg <- mean(N_Supply_day)

Soildata$NSupply12yr[i] <- N_Supply_avg

#4yr
N_output <- apsim_output_4yr[complete.cases(apsim_output_4yr[ ,"N_Uptake"]) & apsim_output_4yr$N_Uptake != 0, ]
N_Supply_day <- pmax(0, pmin(1, (N_output$Nmin_som + N_output$Nmin_res) / (N_output$N_Uptake * 10)))
N_Supply_avg <- mean(N_Supply_day)

Soildata$NSupply4yr[i] <- N_Supply_avg

#############################
#Water supply to plants (mm per yr)
#############################

#25yr
Water_supply <- apsim_output_25yr %>%
  group_by(Year)%>%
  summarize(Year_watersupply = sum(Transpiration))
Water_supply <- Water_supply[-1, ]
Water_supply <- mean(Water_supply$Year_watersupply, na.rm = T)

Soildata$Water_Supply25yr[i] <- Water_supply

#12yr
Water_supply <- apsim_output_12yr %>%
  group_by(Year)%>%
  summarize(Year_watersupply = sum(Transpiration))
Water_supply <- Water_supply[-1, ]
Water_supply <- mean(Water_supply$Year_watersupply, na.rm = T)

Soildata$Water_Supply12yr[i] <- Water_supply

#4yr
Water_supply <- apsim_output_4yr %>%
  group_by(Year)%>%
  summarize(Year_watersupply = sum(Transpiration))
Water_supply <- Water_supply[-1, ]
Water_supply <- mean(Water_supply$Year_watersupply, na.rm = T)

Soildata$Water_Supply4yr[i] <- Water_supply

######################################
#Carbon Sequestration (T Carbon per yr)
######################################

#25yr
Carbon_stocks <- apsim_output_25yr %>%
  group_by(Year)%>%
  summarize(Year_carbon = (OC_Kgha[1] - OC_Kgha[365]) * 0.001)
CO2e_seq <- mean(Carbon_stocks$Year_carbon, na.rm = T)
Soildata$CO2e_seq25yr[i] <- CO2e_seq

#12yr
Carbon_stocks <- apsim_output_12yr %>%
  group_by(Year)%>%
  summarize(Year_carbon = (OC_Kgha[1] - OC_Kgha[365]) * 0.001)
CO2e_seq <- mean(Carbon_stocks$Year_carbon, na.rm = T)
Soildata$CO2e_seq12yr[i] <- CO2e_seq

#4yr
Carbon_stocks <- apsim_output_4yr %>%
  group_by(Year)%>%
  summarize(Year_carbon = (OC_Kgha[1] - OC_Kgha[365]) * 0.001)
CO2e_seq <- mean(Carbon_stocks$Year_carbon, na.rm = T)
Soildata$CO2e_seq4yr[i] <- CO2e_seq

###############################
#N2O emissions (Kg CO2e per yr)
###############################

#25yr
N2O_CO2e <- sum(apsim_output_25yr$N2O) * 273 / 25000
Soildata$N2O_CO2e25yr[i] <- N2O_CO2e

#12yr
N2O_CO2e <- sum(apsim_output_12yr$N2O) * 273 / 25000
Soildata$N2O_CO2e12yr[i] <- N2O_CO2e

#4yr
N2O_CO2e <- sum(apsim_output_4yr$N2O) * 273 / 25000
Soildata$N2O_CO2e4yr[i] <- N2O_CO2e

###############################
#Water quality (ratio per yr)
###############################

#25yr
Water_quality_day <- ((1- apsim_output_25yr$N_leaching) / apsim_output_25yr$N_tot)
Water_quality_avg <- mean(Water_quality_day)

Soildata$Water_quality25yr[i] <- Water_quality_avg

#12yr
Water_quality_day <- ((1- apsim_output_12yr$N_leaching) / apsim_output_12yr$N_tot)
Water_quality_avg <- mean(Water_quality_day)

Soildata$Water_quality12yr[i] <- Water_quality_avg

#4yr
Water_quality_day <- ((1- apsim_output_4yr$N_leaching) / apsim_output_4yr$N_tot)
Water_quality_avg <- mean(Water_quality_day)

Soildata$Water_quality4yr[i] <- Water_quality_avg

####################################
#Groundwater recharge (mm per yr)
####################################

#25yr
GW_recharge <- apsim_output_25yr%>%
  group_by(Year)%>%
  summarize(Drainage_year = sum(Drainage, na.rm = TRUE))

Soildata$GW_recharge25yr[i] <- mean(GW_recharge$Drainage_year)

#12yr
GW_recharge <- apsim_output_12yr%>%
  group_by(Year)%>%
  summarize(Drainage_year = sum(Drainage, na.rm = TRUE))

Soildata$GW_recharge12yr[i] <- mean(GW_recharge$Drainage_year)

#4yr
GW_recharge <- apsim_output_4yr%>%
  group_by(Year)%>%
  summarize(Drainage_year = sum(Drainage, na.rm = TRUE))

Soildata$GW_recharge4yr[i] <- mean(GW_recharge$Drainage_year)

##############################
#Flood Mitigation (mm per yr)
##############################

#25yr
Flood_mitigation <- apsim_output_25yr %>%
  group_by(Year)%>%
  summarize(Year_Flood = sum(Infiltration, na.rm = TRUE))
Soildata$Flood_mitigation25yr[i] <- mean(Flood_mitigation$Year_Flood)

#12yr
Flood_mitigation <- apsim_output_12yr %>%
  group_by(Year)%>%
  summarize(Year_Flood = sum(Infiltration, na.rm = TRUE))
Soildata$Flood_mitigation12yr[i] <- mean(Flood_mitigation$Year_Flood)

#25yr
Flood_mitigation <- apsim_output_4yr %>%
  group_by(Year)%>%
  summarize(Year_Flood = sum(Infiltration, na.rm = TRUE))
Soildata$Flood_mitigation4yr[i] <- mean(Flood_mitigation$Year_Flood)

##################
#Wheat Yield
##################

#25yr
Wheat_yield <- apsim_output_25yr[apsim_output_25yr$Wheat_stage == "EndGrainFill", ]
Wheat_yield_tot <- sum(Wheat_yield$Yield)
Soildata$Yield25yr[i] <- Wheat_yield_tot

#12yr
Wheat_yield <- apsim_output_12yr[apsim_output_12yr$Wheat_stage == "EndGrainFill", ]
Wheat_yield_tot <- sum(Wheat_yield$Yield)
Soildata$Yield12yr[i] <- Wheat_yield_tot

#4yr
Wheat_yield <- apsim_output_4yr[apsim_output_4yr$Wheat_stage == "EndGrainFill", ]
Wheat_yield_tot <- sum(Wheat_yield$Yield)
Soildata$Yield4yr[i] <- Wheat_yield_tot

}

#Save for type of intervention

write.csv(Soildata, file = "Outputs/Conventional_ES.csv", row.names = FALSE)

##########################################
#---------------ANALYSIS-----------------#
##########################################

#Read data from disk *AND CORRECT FOR CARBON AND YIELD MISCALCULATIONS...*

CON <- as.data.frame(read.csv("Outputs/Conventional_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)
LEY <- as.data.frame(read.csv("Outputs/Ley-arable_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)
WCV <- as.data.frame(read.csv("Outputs/Winter-cover_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)
STW <- as.data.frame(read.csv("Outputs/Straw_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)
FYM <- as.data.frame(read.csv("Outputs/Farmyard-manure_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)
LEG <- as.data.frame(read.csv("Outputs/Legume_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)
COM <- as.data.frame(read.csv("Outputs/Combined_ES.csv")) %>%
  mutate(CO2e_seq25yr = CO2e_seq25yr * -1, CO2e_seq12yr = CO2e_seq12yr * -1, CO2e_seq4yr = CO2e_seq4yr * -1,
         Yield25yr = Yield25yr / 25000, Yield12yr = Yield12yr / 12000, Yield4yr = Yield4yr / 4000)

#determine all cells with na values in any dataset
CON_na <- CON[rowSums(is.na(CON)) > 0 , ] 
LEY_na <- LEY[rowSums(is.na(LEY)) > 0 , ]
WCV_na <- WCV[rowSums(is.na(WCV)) > 0 , ] 
STW_na <- STW[rowSums(is.na(STW)) > 0 , ] 
FYM_na <- FYM[rowSums(is.na(FYM)) > 0 , ] 
LEG_na <- LEG[rowSums(is.na(LEG)) > 0 , ]
COM_na <- COM[rowSums(is.na(COM)) > 0 , ] 

To_remove <- merge(CON_na, LEY_na, all = TRUE)%>%
  merge(WCV_na, all = TRUE)%>%
  merge(STW_na, all = TRUE)%>%
  merge(FYM_na, all = TRUE)%>%
  merge(LEG_na, all = TRUE)%>%
  merge(COM_na, all = TRUE)

#remove na values and reduce to 1000
CON <- CON[!CON$cell %in% To_remove$cell, ]%>%
  slice(1:1000)
LEY <- LEY[!LEY$cell %in% To_remove$cell, ]%>%
  slice(1:1000)
WCV <- WCV[!WCV$cell %in% To_remove$cell, ]%>%
  slice(1:1000)
STW <- STW[!STW$cell %in% To_remove$cell, ]%>%
  slice(1:1000)
FYM <- FYM[!FYM$cell %in% To_remove$cell, ]%>%
  slice(1:1000)
LEG <- LEG[!LEG$cell %in% To_remove$cell, ]%>%
  slice(1:1000)
COM <- COM[!COM$cell %in% To_remove$cell, ]%>%
  slice(1:1000)

################
#data processing
################

#normalize by min/max scaling 
all_data <- rbind(CON = CON,
                  LEY = LEY,
                  WCV = WCV,
                  STW = STW,
                  FYM = FYM,
                  LEG = LEG,
                  COM = COM
                  )

normalise <- function(x) {(x - min(x)) / (max(x) - min(x))}

all_data[] <- all_data %>% lapply(normalise)

CON_norm <- all_data[grepl("CON", rownames(all_data)), ]
LEY_norm <- all_data[grepl("LEY", rownames(all_data)), ]
WCV_norm <- all_data[grepl("WCV", rownames(all_data)), ]
STW_norm <- all_data[grepl("STW", rownames(all_data)), ]
FYM_norm <- all_data[grepl("FYM", rownames(all_data)), ]
LEG_norm <- all_data[grepl("LEG", rownames(all_data)), ]
COM_norm <- all_data[grepl("COM", rownames(all_data)), ]


# Calculate relative differences for absolute values
rel_diff_abs <- data.frame(
  LEY = (abs(LEY[,14:ncol(CON)]) - abs(CON[,14:ncol(CON)])) / abs(CON[,14:ncol(CON)]),
  WCV = (abs(WCV[,14:ncol(CON)]) - abs(CON[,14:ncol(CON)])) / abs(CON[,14:ncol(CON)]),
  STW = (abs(STW[,14:ncol(CON)]) - abs(CON[,14:ncol(CON)])) / abs(CON[,14:ncol(CON)]),
  FYM = (abs(FYM[,14:ncol(CON)]) - abs(CON[,14:ncol(CON)])) / abs(CON[,14:ncol(CON)]),
  LEG = (abs(LEG[,14:ncol(CON)]) - abs(CON[,14:ncol(CON)])) / abs(CON[,14:ncol(CON)]),
  COM = (abs(COM[,14:ncol(CON)]) - abs(CON[,14:ncol(CON)])) / abs(CON[,14:ncol(CON)])
)

#summarise medians
rel_diff_abs_medians <- data.frame(
  LEY = apply(rel_diff_abs[, grepl("LEY", colnames(rel_diff_abs))], 2, median), 
  WCV = apply(rel_diff_abs[, grepl("WCV", colnames(rel_diff_abs))], 2, median),
  STW = apply(rel_diff_abs[, grepl("STW", colnames(rel_diff_abs))], 2, median),
  FYM = apply(rel_diff_abs[, grepl("FYM", colnames(rel_diff_abs))], 2, median),
  LEG = apply(rel_diff_abs[, grepl("LEG", colnames(rel_diff_abs))], 2, median),
  COM = apply(rel_diff_abs[, grepl("COM", colnames(rel_diff_abs))], 2, median)
) %>% rownames_to_column("Ecosystem_Service")
rel_diff_abs_medians$Ecosystem_Service <- sub('^LEY.', '', rel_diff_abs_medians$Ecosystem_Service) 


#Standard error
rel_diff_abs_se <- data.frame(
  LEY = apply(LEY [,14:ncol(CON )], 2, sd) / sqrt(nrow(LEY )),
  WCV = apply(WCV [,14:ncol(CON )], 2, sd) / sqrt(nrow(WCV )),
  STW = apply(STW [,14:ncol(CON )], 2, sd) / sqrt(nrow(STW )),
  FYM = apply(FYM [,14:ncol(CON )], 2, sd) / sqrt(nrow(FYM )),
  LEG = apply(LEG [,14:ncol(CON )], 2, sd) / sqrt(nrow(LEG )),
  COM = apply(COM [,14:ncol(CON )], 2, sd) / sqrt(nrow(COM ))
)%>% rownames_to_column("Ecosystem_Service")

# Perform t-tests and calculate p-values
pvals_abs <- data.frame(
  LEY = apply(LEY [,14:ncol(CON )], 2, function(x) t.test(x, CON [,14:ncol(CON )])$p.value),
  WCV = apply(WCV [,14:ncol(CON )], 2, function(x) t.test(x, CON [,14:ncol(CON )])$p.value),
  STW = apply(STW [,14:ncol(CON )], 2, function(x) t.test(x, CON [,14:ncol(CON )])$p.value),
  FYM = apply(FYM [,14:ncol(CON )], 2, function(x) t.test(x, CON [,14:ncol(CON )])$p.value),
  LEG = apply(LEG [,14:ncol(CON )], 2, function(x) t.test(x, CON [,14:ncol(CON )])$p.value),
  COM = apply(COM [,14:ncol(CON )], 2, function(x) t.test(x, CON [,14:ncol(CON )])$p.value)
)%>% rownames_to_column("Ecosystem_Service")

#wide to long for absolute and normalised
rel_diff_abs_medians_long <- rel_diff_abs_medians %>% 
  gather(key = "Intervention", value = "ES_change_median", 2:ncol(rel_diff_abs_medians), factor_key = TRUE)

rel_diff_abs_se_long <- rel_diff_norm_se %>% 
  gather(key = "Intervention", value = "se", 2:ncol(rel_diff_abs_se), factor_key = TRUE)

pvals_abs_long <- pvals_abs %>% 
  gather(key = "Intervention", value = "pval", 2:ncol(pvals_abs), factor_key = TRUE)
pvals_abs_long <- pvals_abs_long %>% mutate(sig = case_when(pvals_abs_long$pval > 0.05 ~ "*", pvals_abs_long$pval <= 0.05 ~ ""))

#summarise
summary_abs <- merge(rel_diff_abs_medians_long, rel_diff_abs_se_long, all = TRUE)%>%
  merge(pvals_abs_long, all = TRUE) %>% 
  mutate(ES_change_median = if_else(Ecosystem_Service == "N2O_CO2e25yr", ES_change_median * -1, ES_change_median)) %>%
  mutate(ES_change_median = if_else(Ecosystem_Service == "N2O_CO2e12yr", ES_change_median * -1, ES_change_median)) %>%
  mutate(ES_change_median = if_else(Ecosystem_Service == "N2O_CO2e4yr", ES_change_median * -1, ES_change_median))

summary_abs$duration_yr <-ifelse(grepl("25yr", summary_abs$Ecosystem_Service), "25",
                          ifelse(grepl("12yr", summary_abs$Ecosystem_Service), "12",
                          ifelse(grepl("4yr", summary_abs$Ecosystem_Service), "4", NA)))

summary_abs$Ecosystem_Service<-ifelse(grepl("CO2e_seq", summary_abs$Ecosystem_Service), "CO2e_seq",
                                ifelse(grepl("Flood_mitigation", summary_abs$Ecosystem_Service), "Flood_mitigation",
                                ifelse(grepl("GW_recharge", summary_abs$Ecosystem_Service), "GW_recharge",
                                ifelse(grepl("N2O_CO2e", summary_abs$Ecosystem_Service), "N2O_CO2e",
                                ifelse(grepl("NSupply", summary_abs$Ecosystem_Service), "NSupply",
                                ifelse(grepl("Water_quality", summary_abs$Ecosystem_Service), "Water_quality",
                                ifelse(grepl("Water_Supply", summary_abs$Ecosystem_Service), "Water_Supply",
                                ifelse(grepl("Yield", summary_abs$Ecosystem_Service), "Yield", NA))))))))


summary_abs_25 <- summary_abs[summary_abs$duration_yr == "25",]

#####################################
#Correlation matrix
#####################################

#make names readable
corr_data <- all_data[, c(5:8, 10, 14, 17, 20, 23, 26, 29, 32, 35)]
colnames(corr_data) <- c("Bulk Density",
             "Carbon: Nitrogen",
             "SOC",
             "Soil pH",
             "Field Capacity",
             "N Supply",
             "Water Supply",
             "C Stocks",
             "Nitrous Oxide Emissions",
             "Nitrate Retention",
             "Drainage",
             "Infiltration",
             "Total Yield"
             )

#make matrix and p values
cor <- cor_mat(corr_data) %>% column_to_rownames(var = "rowname") %>% as.matrix()
cor_p <- cor_pmat(corr_data) %>% column_to_rownames(var = "rowname") %>% as.matrix()

#plot
png(height=800, width=1000, file="Outputs/ES_correlations_median.png", type = "cairo")
corr_plot <- corrplot(cor,
                  method="color",
                  type="upper",
                  is.corr=TRUE,
                  col.lim=NULL,
                  diag=FALSE,
                  outline=FALSE,
                  bg="white",
                  p.mat=cor_p,
                  tl.pos="td",
                  tl.col="black",
                  tl.srt=45,
                  tl.cex=1.8,
                  cl.pos="r",
                  cl.ratio=0.2,
                  cl.offset=0,
                  sig.level=0.05,
                  cl.length=11,
                  cl.cex=1.5,
                  pch.col="darkgrey",
                  na.label="NA",
                  na.label.col="grey",
                  fin=c(1,1),
                  new=TRUE
                  )
dev.off()

corr_plot

####################################
#All ES comparison
#####################################

#filter for only 25yrs simulation
comp_matrix <- summary_abs_25 %>% 
  select(Ecosystem_Service, Intervention, ES_change_median) %>%
  spread(key = "Intervention", value = "ES_change_median") %>%
  select(!Ecosystem_Service)%>%
  mutate(AVG = rowMeans(.)) %>%
  as.matrix()
rownames(comp_matrix) <- c(
                     "Climate regulation (Carbon Stocks Changes)",
                     "Flood Mitigation",
                     "Groundwater Supply",
                     "Climate regulation (Nitrous oxide retained)",
                     "Crop provisioning (Nitrogen supply)",
                     "Abiotic filtration and sequestration",
                     "Crop provisioning (Water Supply)",
                     "Total Yield"
)

comp_matrix <- comp_matrix[order(row.names(comp_matrix)), ]

#p-values
comp_matrix_pvals <- pvals_abs %>% 
  filter(grepl("25", Ecosystem_Service)) %>% 
  column_to_rownames(var = "Ecosystem_Service") %>%
  mutate(AVG = rep(0, 8)) %>%
  as.matrix()
rownames(comp_matrix_pvals) <- c(
                     "Crop provisioning (Nitrogen supply)",
                     "Crop provisioning (Water Supply)",
                     "Climate regulation (Carbon Stocks Changes)",
                     "Climate regulation (Nitrous oxide retained)",
                     "Abiotic filtration and sequestration",
                     "Groundwater Supply",
                     "Flood Mitigation",
                     "Total Yield"
)
comp_matrix_pvals <- comp_matrix_pvals[order(row.names(comp_matrix_pvals)), ]

#standard deviations
comp_matrix_se <- rel_diff_abs_se %>%
  filter(grepl("25", rel_diff_abs_se$Ecosystem_Service)) %>%
  select(!"Ecosystem_Service")
rownames(comp_matrix_se) <- c(
  "Crop provisioning (Nitrogen supply)",
  "Crop provisioning (Water Supply)",
  "Climate regulation (Carbon Stocks Changes)",
  "Climate regulation (Nitrous oxide retained)",
  "Abiotic filtration and sequestration",
  "Groundwater Supply",
  "Flood Mitigation",
  "Total Yield"
)

comp_matrix_sd[] <- comp_matrix_se %>%
  lapply(function(x) x * sqrt(1000)) %>%
  as.data.frame()

#colour
col <- colorRampPalette(c("#EA8300", "#FF8F00", "#FF9A19", "#FFA532", "#FFB14C",
                          "#f1f1f1",
                          "#4CC9FF", "#32C1FF", "#19B9FF", "#00B2FF", "#00A1E7"))

#plot
png(height=1000, width=1400, file="Outputs/comp_matrix_median.png", type = "cairo")
rel_plot <- corrplot(
                    comp_matrix,
                    method="color",
                    type="full",
                    is.corr=FALSE,
                    col.lim=c(-2, 2),
                    col=col(300),
                    outline=FALSE,
                    bg="white",
                    addgrid.col=NA,
                    p.mat=comp_matrix_pvals,
                    tl.col="black",
                    tl.srt=45,
                    tl.cex=2.5,
                    tl.pos="lt",
                    cl.pos="r",
                    cl.ratio=0.25,
                    cl.offset=0.5,
                    sig.level=0.001,
                    cl.length=5,
                    cl.cex=2,
                    pch.col="black",
                    na.label="NA",
                    na.label.col="grey",
                    cl.align.text="l",
                    addCoef.col="black",
                    addCoefasPercent=TRUE,
                    number.cex=2.5,
                    mar=c(0,0,0,0),
                    number.font = 1,
                    number.digits = 1
)
dev.off()

##############################
#absolute box plots
############################

#get data in a readable format
box_dataCONLEY <- merge(CON, LEY, by = c("cell", "x", "y", "ww_uk_eucropmap", "BD", "CN", "C", "pH", "sat", "dul", "ll15", "Lon", "Lat"),
                  all = TRUE,
                  suffixes = c(".CON", ".LEY"))
box_dataWCVFYM <- merge(WCV, FYM, by = c("cell", "x", "y", "ww_uk_eucropmap", "BD", "CN", "C", "pH", "sat", "dul", "ll15", "Lon", "Lat"),
                        all = TRUE,
                        suffixes = c(".WCV", ".FYM"))
box_dataSTWLEG <- merge(STW, LEG, by = c("cell", "x", "y", "ww_uk_eucropmap", "BD", "CN", "C", "pH", "sat", "dul", "ll15", "Lon", "Lat"),
                        all = TRUE,
                        suffixes = c(".STW", ".LEG"))
box_dataSTWLEGCOM <- merge(box_dataSTWLEG, COM, by = c("cell", "x", "y", "ww_uk_eucropmap", "BD", "CN", "C", "pH", "sat", "dul", "ll15", "Lon", "Lat"),
                          all = TRUE)
box_data <- cbind(box_dataCONLEY, box_dataWCVFYM, box_dataSTWLEG, COM)
colnames(box_data) <- make.unique(names(box_data))
box_data <- select(box_data, -contains("cell"))
box_data <-   select(box_data, -contains("x"))
box_data <-   select(box_data, -contains("y."), -"y")
box_data <-   select(box_data, -contains("ww_uk_eucropmap"))
box_data <-   select(box_data, -contains("BD"))
box_data <-   select(box_data, -contains("CN"))
box_data <-   select(box_data, -contains("C."), -"C")
box_data <-   select(box_data, -contains("pH"))
box_data <-   select(box_data, -contains("sat"))
box_data <-   select(box_data, -contains("dul"))
box_data <-   select(box_data, -contains("ll15"))
box_data <-   select(box_data, -contains("Lon"))
box_data <-   select(box_data, -contains("Lat"))

box_data <- box_data %>%
  gather(key = "Name", value = "value")

box_data$Ecosystem_Service <- ifelse(grepl("CO2e_seq", box_data$Name), "CO2e_seq",
                              ifelse(grepl("Flood_mitigation", box_data$Name), "Flood_mitigation",
                              ifelse(grepl("GW_recharge", box_data$Name), "GW_recharge",
                              ifelse(grepl("N2O_CO2e", box_data$Name), "N2O_CO2e",
                              ifelse(grepl("NSupply", box_data$Name), "NSupply",
                              ifelse(grepl("Water_quality", box_data$Name), "Water_quality",
                              ifelse(grepl("Water_Supply", box_data$Name), "Water_Supply",
                              ifelse(grepl("Yield", box_data$Name), "Yield",
                              NA))))))))

box_data$Intervention <- ifelse(grepl("CON", box_data$Name), "CON",
                         ifelse(grepl("LEY", box_data$Name), "LEY",
                         ifelse(grepl("WCV", box_data$Name), "WCV",
                         ifelse(grepl("STW", box_data$Name), "STW",
                         ifelse(grepl("FYM", box_data$Name), "FYM",
                         ifelse(grepl("LEG", box_data$Name), "LEG",
                         "COM"))))))

box_data$Year <- ifelse(grepl("25yr", box_data$Name), "25",
                        ifelse(grepl("12yr", box_data$Name), "12",
                               ifelse(grepl("4yr", box_data$Name), "4",
                                      NA)))

box_data_25 <- box_data %>% filter(Year == "25") %>% select(!"Name")

#Separate out ESs for plotting and order nicely
custom_order <- c("CON", "LEY", "WCV", "STW", "FYM", "LEG", "COM")

box_data_25_CStocks <- box_data_25 %>% filter(Ecosystem_Service == "CO2e_seq")
box_data_25_CStocks$Intervention_ordered <- factor(box_data_25_CStocks$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_CStocks <- arrange(box_data_25_CStocks, Intervention_ordered)

box_data_25_Flood <- box_data_25 %>% filter(Ecosystem_Service == "Flood_mitigation")
box_data_25_Flood$Intervention_ordered <- factor(box_data_25_Flood$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_Flood <- arrange(box_data_25_Flood, Intervention_ordered)

box_data_25_GW <- box_data_25 %>% filter(Ecosystem_Service == "GW_recharge")
box_data_25_GW$Intervention_ordered <- factor(box_data_25_GW$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_GW <- arrange(box_data_25_GW, Intervention_ordered)

box_data_25_N2O <- box_data_25 %>% filter(Ecosystem_Service == "N2O_CO2e")
box_data_25_N2O$Intervention_ordered <- factor(box_data_25_N2O$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_N2O <- arrange(box_data_25_N2O, Intervention_ordered)

box_data_25_NSupply <- box_data_25 %>% filter(Ecosystem_Service == "NSupply")
box_data_25_NSupply$Intervention_ordered <- factor(box_data_25_NSupply$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_NSupply <- arrange(box_data_25_NSupply, Intervention_ordered)

box_data_25_WQuality <- box_data_25 %>% filter(Ecosystem_Service == "Water_quality")
box_data_25_WQuality$Intervention_ordered <- factor(box_data_25_WQuality$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_WQuality <- arrange(box_data_25_WQuality, Intervention_ordered)

box_data_25_WSupply <- box_data_25 %>% filter(Ecosystem_Service == "Water_Supply")
box_data_25_WSupply$Intervention_ordered <- factor(box_data_25_WSupply$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_WSupply <- arrange(box_data_25_WSupply, Intervention_ordered)

box_data_25_Yield <- box_data_25 %>% filter(Ecosystem_Service == "Yield")
box_data_25_Yield$Intervention_ordered <- factor(box_data_25_Yield$Intervention, levels = custom_order, ordered = TRUE)
box_data_25_Yield <- arrange(box_data_25_Yield, Intervention_ordered)


#plot box plots for each ES
CStocks <- ggplot(box_data_25_CStocks, aes(x=Intervention_ordered, y=value)) +
        geom_boxplot(fill = "#907247", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
        geom_boxplot(color = "#907247",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
        ylim(-0.6, 0.6) +
        labs(x = "Intervention", y = expression("Change in C Stocks (Mg C yr"^-1*" ha"^-1*")")) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
          panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
          panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
          panel.background = element_rect(fill = "#F5F5F5"),
          axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5,),
          plot.margin = margin(1, 1, 1, 1, "cm")
          )

CStocks

Flood <- ggplot(box_data_25_Flood, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "lightblue", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "lightblue",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  ylim(600, 1000) +
  labs(x = "Intervention", y = expression("Infiltration (mm yr"^-1*")")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
Flood

GW <- ggplot(box_data_25_GW, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "#878ED6", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "#878ED6",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  ylim(0, 700) +
  labs(x = "Intervention", y = expression("Bottom layer drainage (mm yr"^-1*")")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
GW

N2O <- ggplot(box_data_25_N2O, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "#DA71D2", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "#DA71D2",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  ylim(0, 2.5) +
  labs(x = "Intervention", y = expression("N"[2]*"O emissions (Mg CO"[2]*"e yr"^-1*" ha"^-1*")")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
N2O

NSupply <- ggplot(box_data_25_NSupply, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "#6CEB71", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "#6CEB71",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  ylim(0, 0.4) +
  labs(x = "Intervention", y = "Mineral nitrogen supply to wheat") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
NSupply

WQuality <- ggplot(box_data_25_WQuality, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "#D0EB6C", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "#D0EB6C",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  ylim(0, 0.075) +
  labs(x = "Intervention", y = expression("NO"[3]*" retention capacity")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
WQuality

WSupply <- ggplot(box_data_25_WSupply, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "#EBBF6C", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "#EBBF6C",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  ylim(0, 400) +
  labs(y = expression("Wheat evapotranspiration (mm yr "^-1*")")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
WSupply

Yield <- ggplot(box_data_25_Yield, aes(x=Intervention_ordered, y=value)) +
  geom_boxplot(fill = "#F0D840", outlier.alpha = 0.7, outlier.colour = "black", outlier.stroke = 0, outlier.size = 0.9, lwd = 0.9, fatten = 1.5) + 
  geom_boxplot(color = "#F0D840",fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F, lwd = 0.905) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  labs(y = expression("Wheat yield (Mg ha "^-1*" yr"^-1*")")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", face = "bold", vjust = 0),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "#F5F5F5"),
    axis.ticks = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
Yield

#custom legend
legend_data <- box_data_25
legend <- ggplot(box_data_25, aes(x=value, y=Intervention, fill = Ecosystem_Service)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#907247", "lightblue", "#878ED6", "#DA71D2", "#6CEB71", "#D0EB6C", "#EBBF6C", "#F0D840"),
                       labels = c("Climate Mitigation (C Stocks)", "Flood Mitigation", "Groundwater Recharge", "Climate Regulation (Nitrous Oxide Emissions)", "Crop Provisioning (N Supply)", "Abiotic Filtration, Sequestration, and Storage of Waste", "Crop Provisioning (Water Supply)", "Yield"),
                       name = "Key") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm'),
    legend.box.margin = margin(1, 1, 1, 4, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20, face = "bold"),
    legend.key = element_blank(),
    legend.key.width = unit(2, "cm")
  ) +
  geom_boxplot(color = "white",fatten = NULL, fill = "white", coef = 0, outlier.alpha = 0, show.legend = FALSE)
legend

#Put into one figure
all_boxplots <- ggarrange(Yield, CStocks, N2O, Flood, GW, NSupply, WSupply, WQuality,
                labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                ncol = 3,
                nrow = 3,
                font.label = list(size = 24)
                )

ggsave(
  "Outputs/all_boxplots_ES_3x3.png",
  plot = all_boxplots,
  device = "png",
  scale = 1.8,
  width = 28,
  height = 24,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = "white"
  )

#######################
#ES Change over time
#######################
time_data <- summary_abs %>%
  filter(Intervention == "WCV", Ecosystem_Service != "CO2e_seq")

time <- ggplot(time_data, aes(x=duration_yr, y=ES_change_median, group = Ecosystem_Service)) +
  geom_line(aes(colour = Ecosystem_Service)) +
  scale_colour_manual(values = c("#907247", "lightblue", "#878ED6", "#DA71D2", "#6CEB71", "#D0EB6C", "#EBBF6C", "#F0D840"),
                    labels = c("Climate Mitigation (C Stocks)", "Flood Mitigation", "Groundwater Recharge", "Climate Regulation (Nitrous Oxide Retention", "Crop Provisioning (N Supply)", "Abiotic Filtration, Sequestration, and Storage of Waste", "Crop Provisioning (Water Supply", "Yield"),
                    name = "Key") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
time


#######################
#Wheat raster
#######################

wheat_rast_export <- rast("Soil data/Land use 1km 2015 RASTER/ww_uk_eucropmap.tif") %>%
  project("epsg:27700")


