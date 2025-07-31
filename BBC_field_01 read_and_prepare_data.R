library(neonUtilities) # functions to download portal data
library(tidyverse)

# load data from NEON portal using loadByProduct function
# bbcDatField <- try(neonUtilities::loadByProduct(
#   dpID = 'DP1.10067.001',
#   check.size = F,
#   site = "all",
#   startdate = "2012-01",
#   enddate = "2024-06", ,
#   release = "current",
#   package = "basic",
#   include.provisional = TRUE,
#   token = Sys.getenv('NEON_PAT')
# ),silent=T)
#saveRDS(bbcDatField, "bbcDatField.rds") 
# originally loaded and saved provisional data on April 19, 2024
# re-loaded and re-saved provisional data on July 26, 2025 and re-ran through all field analysis code-chunks

# read in saved data
bbcDatField <- readRDS("bbcDatField.rds")


list2env(bbcDatField, envir=.GlobalEnv) 

# rootSampleArea unit is squareMeter
# dryMass unit is gram

# prepare data
bbc_percore_mod <- bbc_percore

bbc_percore_mod$year <- substr(bbc_percore_mod$eventID,10,13)

bbc_percore_mod <- bbc_percore_mod %>% filter(sampleID != "BBC.YELL055671.20200813.NORTH") # remove the one exraneousrecord from this bout

stature_class <- bbc_percore_mod %>% distinct(eventID, plotID) %>% group_by(eventID) %>% summarise(n = n() ) %>% filter(n > 20) %>%
   mutate(siteID = substr(eventID, 5,8)) %>% distinct(siteID) %>% mutate(stature = "small")
year_class <- bbc_percore_mod %>% distinct(eventID, siteID) %>% group_by(siteID) %>% summarise(years = n() )
stature_year <- merge(year_class, stature_class, by = c("siteID"), all.x = T)
  stature_year$stature <- if_else(stature_year$stature == "small", stature_year$stature, "large", "large")
stature_year$stature_yr <- paste0(stature_year$stature, "_", stature_year$years,"times")

bbc_percore_mod$nlcdClass <- ifelse(bbc_percore_mod$plotID == "SERC_061", "deciduousForest", bbc_percore_mod$nlcdClass) # nlcdClass missing but other plots at site had decidiousForest

bbc_percore_plotCount <- bbc_percore_mod %>% distinct(eventID, siteID, plotID) %>% group_by(eventID, siteID) %>% summarise(plot_count = n()) %>% select(-siteID)
bbc_percore_clipIDCount <- bbc_percore_mod %>% filter(samplingImpractical == "OK") %>% distinct(eventID, siteID, plotID, clipID) %>% group_by(eventID,plotID) %>% summarise(clip_count = n())
bbc_percore_coreIDCount <- bbc_percore_mod %>% filter(samplingImpractical == "OK") %>% group_by(eventID, siteID, clipID) %>% summarise(core_count = n()) %>% ungroup() %>% select(-siteID)
bbc_percore_nlcdCount <- bbc_percore_mod %>% filter(samplingImpractical == "OK") %>% distinct(eventID, siteID, nlcdClass) %>% group_by(eventID, siteID) %>% summarise(nlcd_count = n()) %>% ungroup() %>% select(-siteID)

bbc_percore_mod <- bbc_percore_mod %>% filter(samplingImpractical == "OK") %>% left_join(bbc_percore_plotCount, by = "eventID") %>% left_join(bbc_percore_clipIDCount, by = c("eventID","plotID")) %>% 
  left_join(bbc_percore_coreIDCount, by = c("eventID","clipID")) %>% left_join(bbc_percore_nlcdCount, by = "eventID") 

# get the average of all cores and sizes within a clipID on a particular datebbc_rootChemistry_mod <- bbc_rootChemistry 
bbc_rootChemistry_mod$cnSampleID <- gsub("0-05.POOL.", "0-1.POOL.", bbc_rootChemistry_mod$cnSampleID)
bbc_rootChemistry_mod$cnSampleID <- gsub("05-1.POOL.", "0-1.POOL.", bbc_rootChemistry_mod$cnSampleID)
bbc_rootChemistry_mod <- bbc_rootChemistry_mod %>% group_by(cnSampleID) %>% summarize(nitrogenPercent = mean(nitrogenPercent)) %>% ungroup() # get the average of analyticalReps
bbc_rootChemistry_mod$sampleID_simple <- substr(bbc_rootChemistry_mod$cnSampleID, 1,23)
bbc_rootChemistry_clipMean <- bbc_rootChemistry_mod %>% filter(!is.na(nitrogenPercent)) %>% group_by(sampleID_simple) %>% summarize(nitrogenPercent_ave = mean(nitrogenPercent)) 


# get the average of all cores and sizes within a clipID on a particular datebbc_rootChemistry_lumped <- bbc_rootChemistry 
bbc_rootChemistry_lumped$cnSampleID <- gsub("0-05.POOL.", "0-2.POOL.", bbc_rootChemistry_lumped$cnSampleID)
bbc_rootChemistry_lumped$cnSampleID <- gsub("05-1.POOL.", "0-2.POOL.", bbc_rootChemistry_lumped$cnSampleID)
bbc_rootChemistry_lumped$cnSampleID <- gsub("1-2.POOL.", "0-2.POOL.", bbc_rootChemistry_lumped$cnSampleID)
bbc_rootChemistry_lumped <- bbc_rootChemistry_lumped %>% group_by(cnSampleID) %>% summarize(nitrogenPercent = mean(nitrogenPercent)) %>% ungroup() # get the average of analyticalReps
bbc_rootChemistry_lumped$sampleID_simple <- substr(bbc_rootChemistry_lumped$cnSampleID, 1,23)
bbc_rootChemistry_clipMean_lumped <- bbc_rootChemistry_lumped %>% filter(!is.na(nitrogenPercent)) %>% group_by(sampleID_simple) %>% summarize(nitrogenPercent_ave = mean(nitrogenPercent)) 


# add dead and live root masses, and combine old 0-05 and 05-1 sizeCategories
bbc_rootMass_grouped <- bbc_rootmass %>% filter(qaDryMass == "N", !is.na(dryMass))
bbc_rootMass_grouped$sampleID_simple <- substr(bbc_rootMass_grouped$sampleID, 1, 23)
bbc_rootMass_grouped$sizeCategory <- ifelse(bbc_rootMass_grouped$sizeCategory == "0-05" | bbc_rootMass_grouped$sizeCategory == "05-1", "0-1", bbc_rootMass_grouped$sizeCategory)
bbc_rootMass_grouped <- bbc_rootMass_grouped %>% group_by(sampleID, sampleID_simple, sizeCategory) %>% summarise(dryMass = sum(dryMass))  
bbc_rootMass_grouped$cnSampleID <- paste0(bbc_rootMass_grouped$sampleID_simple,".",bbc_rootMass_grouped$sizeCategory,".POOL.CN")

bbc_rootMass_grouped_lumped <- bbc_rootMass_grouped
bbc_rootMass_grouped_lumped$sizeCategory <- ifelse(bbc_rootMass_grouped_lumped$sizeCategory == "0-1" | bbc_rootMass_grouped_lumped$sizeCategory == "1-2", "0-2", bbc_rootMass_grouped_lumped$sizeCategory)
bbc_rootMass_grouped_lumped <- bbc_rootMass_grouped_lumped %>% group_by(sampleID, sampleID_simple, sizeCategory) %>% summarise(dryMass = sum(dryMass))  # combine 0-1 and 1-2 sizeCategories
bbc_rootMass_grouped_lumped$cnSampleID <- paste0(bbc_rootMass_grouped_lumped$sampleID_simple,".",bbc_rootMass_grouped_lumped$sizeCategory,".POOL.CN")


# merge mass and core and calculate g per m2
data_merge <- merge(bbc_rootMass_grouped, bbc_percore_mod, by = "sampleID", all.x = TRUE) %>% mutate(mass_g_m2 = dryMass / rootSampleArea) %>% 
   filter(!is.na(mass_g_m2)) %>% 
   select(-c(dryMass, uid, namedLocation, decimalLatitude, decimalLongitude, geodeticDatum, coordinateUncertainty, elevation, samplingProtocolVersion, coreDiameter, monolithLength, monolithWidth, toxicodendronPossible, measuredBy, recordedBy, publicationDate))

# add summed mass as a new column and calc each size category's fraction of total mass for core
data_total_small <- data_merge %>% group_by(sampleID) %>% reframe(mass_g_m2_tot = sum(mass_g_m2)) # sum the masses of separate size categories
data_merge <- data_merge %>% left_join(data_total_small, by = "sampleID") %>% mutate(mass_fract = mass_g_m2/mass_g_m2_tot) 

# add nitrogen percent and then weight that nitrogen value by the fraction of core mass contributed by size category
data_merge <- data_merge %>% left_join(bbc_rootChemistry_mod, by = c("sampleID_simple","cnSampleID")) %>% left_join(bbc_rootChemistry_clipMean, by = "sampleID_simple")
 data_merge$nitrogenPercent <- ifelse(is.na(data_merge$nitrogenPercent), data_merge$nitrogenPercent_ave, data_merge$nitrogenPercent)
data_merge <- data_merge %>% mutate(nitrogen_weight = nitrogenPercent * mass_fract) 

# Sum the masses of sizeCategories within core, and also sum the nitrogen_weights, giving a total mass and appropriately mass weighted average N per for each core
data_total <- data_merge %>% group_by(sampleID, sampleID_simple, siteID, plotID, domainID, subplotID, clipID, coreID, collectDate, eventID, samplingImpractical, sampleCondition, wst10cmDist, wst1cmDist, bareGround,
        litterDepth, rootSamplingMethod, rootSampleArea, rootSampleDepth, remarks, release, nlcdClass, year, plot_count, clip_count, core_count, nlcd_count) %>% # years, stature, stature_yr, 
        reframe(mass_g_m2 = sum(mass_g_m2), nitrogenPercent = sum(nitrogen_weight))
data_total$sizeCategory <- "total"
data_total$cnSampleID <- NA

data_merge$nitrogenPercent_ave <- data_merge$mass_g_m2_tot <- data_merge$mass_fract <- data_merge$nitrogen_weight <- data_merge$elevationUncertainty <- NULL
data_merge <- rbind(data_merge, data_total) # combine the aggregate core mass and N with size-category specific mass and N
data_merge <- data_merge %>% left_join(stature_year, by = "siteID")
data_merge <- data_merge %>% rename(response = mass_g_m2)

#### finish lumping 0-1 cm and 1-2 cm diameter size classes into single 0-2 cm diameter size class for use in subsequent analyses
# merge mass and core and calculate g per m2
data_merge_lumped <- merge(bbc_rootMass_grouped_lumped, bbc_percore_mod, by = "sampleID", all.x = TRUE) %>% mutate(mass_g_m2 = dryMass / rootSampleArea) %>% 
   filter(!is.na(mass_g_m2)) %>% 
   select(-c(dryMass, uid, namedLocation, decimalLatitude, decimalLongitude, geodeticDatum, coordinateUncertainty, elevation, samplingProtocolVersion, coreDiameter, monolithLength, monolithWidth, toxicodendronPossible, measuredBy, recordedBy, publicationDate))
    
# add summed mass as a new column and calc each size categorie's fraction of total mass for core
data_total_small_lumped <- data_merge_lumped %>% group_by(sampleID) %>% reframe(mass_g_m2_tot = sum(mass_g_m2)) # sum the masses of separate size categories
data_merge_lumped <- data_merge_lumped %>% left_join(data_total_small_lumped, by = "sampleID") %>% mutate(mass_fract = mass_g_m2/mass_g_m2_tot) 

# add nitrogen percent and then weight that nitrogen value by the fraction of core mass contributed by size category
data_merge_lumped <- data_merge_lumped %>% left_join(bbc_rootChemistry_lumped, by = c("sampleID_simple","cnSampleID")) %>% left_join(bbc_rootChemistry_clipMean_lumped, by = "sampleID_simple")
 data_merge_lumped$nitrogenPercent <- ifelse(is.na(data_merge_lumped$nitrogenPercent), data_merge_lumped$nitrogenPercent_ave, data_merge_lumped$nitrogenPercent)
data_total$cnSampleID <- NA
 

data_merge_lumped$nitrogenPercent_ave <- data_merge_lumped$mass_g_m2_tot <- data_merge_lumped$mass_fract <- data_merge_lumped$elevationUncertainty <- NULL
data_merge_lumped <- rbind(data_merge_lumped, data_total) # combine the aggregate core mass and N with size-category specific mass and N
data_merge_lumped <- data_merge_lumped %>% left_join(stature_year, by = "siteID")
data_merge_lumped <- data_merge_lumped %>% rename(response = mass_g_m2)