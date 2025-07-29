#' #######################################################################################
#' Analysis Script for BBC (root) lab dilution analysis
#' Sam Simkin (ssimkin@batelleecology) originally authored in 2023, with bootstrapping approach based on Courtney Meier code for stem map sampling simulations
#' #######################################################################################

library(neonUtilities)
library(tidyverse)
library(broom) # includes broom function
library(infer) # includes rep_slice_sample function
library(ggplot2)


bbcDat <- loadByProduct(dpID="DP1.10067.001", package = "basic", check.size = FALSE, token = Sys.getenv('NEON_TOKEN')) # pulled from portal 2023-05-02
 list2env(bbcDat ,.GlobalEnv) # unlist all data frames

bbc_rootmass_lt2 <- bbc_rootmass %>%  filter(sizeCategory != "2-10")
rootmass_by_sampleID <- bbc_rootmass_lt2 %>% group_by(sampleID) %>% summarise(root_dryMass = sum(dryMass))
 
bbc_percore <- bbc_percore %>% select(subplotID, clipID, coreID, sampleID, wst10cmDist, wst1cmDist, litterDepth, rootSampleArea, rootSampleDepth)

rootmass_per_m2 <- merge(rootmass_by_sampleID, bbc_percore, by = "sampleID", all.x = TRUE)
rootmass_per_m2 <- rootmass_per_m2 %>% mutate(root_g_m2 = round(root_dryMass / rootSampleArea , digits = 4) ) %>% select(sampleID, root_dryMass, rootSampleArea, root_g_m2)

bbc_dilution <- merge(bbc_dilution, bbc_percore, by = "sampleID", all.x = TRUE)
bbc_dilution <- bbc_dilution %>% filter(!is.na(rootSampleArea))

bbc_dilution$year <- substr(bbc_dilution$collectDate,1,4)
bbc_dilution$eventID <- paste0(bbc_dilution$siteID, "_", bbc_dilution$year)
bbc_dil_unique_eventID <- bbc_dilution %>% select(siteID, year, eventID) %>% distinct()  
bbc_dilution_yr_count <- bbc_dil_unique_eventID %>% group_by(siteID) %>% tally()
 
bbc_dilution$year <- substr(bbc_dilution$collectDate,1,4)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "BART" & bbc_dilution$year == "2016", 1, bbc_dilution$year)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "BART" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "GUAN" & bbc_dilution$year == "2018", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "GUAN" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "JORN" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "JORN" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "KONZ" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "KONZ" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "MOAB" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "MOAB" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ONAQ" & bbc_dilution$year == "2016", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ONAQ" & bbc_dilution$year == "2021", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ORNL" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ORNL" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "SCBI" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "SCBI" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STEI" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STEI" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STER" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STER" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TALL" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TALL" & bbc_dilution$year == "2021", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TOOL" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TOOL" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "UNDE" & bbc_dilution$year == "2016", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "UNDE" & bbc_dilution$year == "2019", 2, bbc_dilution$year2)
bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "WOOD" & bbc_dilution$year == "2016", 1, bbc_dilution$year2)
 bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "WOOD" & bbc_dilution$year == "2021", 2, bbc_dilution$year2)  

# look at all 47 sites and for the 14 sites with 2 bouts look at both bouts
# filter to remove outliers > 99 percentile for whole dataset and filter to remove fragMass < 0 g
bbc_dilution <- bbc_dilution %>% filter(dryMass >= 0) %>%  filter(dryMass < quantile(dryMass, probs = 0.99, na.rm =TRUE))

# if there are < 10 subsamples per core then drop the entire core
#dil_is_10 <- bbc_dilution %>% group_by(sampleID) %>% summarise(nDM = n()) %>%  filter(nDM == 10)
#  bbc_dilution <- merge(bbc_dilution, dil_is_10, all.y=TRUE)
dil_is_gt7 <- bbc_dilution %>% group_by(sampleID) %>% summarise(nDM = n()) %>%  filter(nDM > 7)
  bbc_dilution <- merge(bbc_dilution, dil_is_gt7, all.y=TRUE)
  bbc_dilution$nDM <- NULL

### Calculate dryMass for fragments < 1 cm length for each dilutionSubsampleID (calculated mass per unit area but could do on a soil volume basis instead)
dil_calc <- bbc_dilution %>%
  mutate(fragMass = round(dryMass*(sampleVolume/dilutionSubsampleVolume), digits = 4), frag_g_m2 = round(fragMass / rootSampleArea , digits = 4),
         somMass = round(somDryMass*(sampleVolume/dilutionSubsampleVolume), digits = 4),  som_g_m2 = round(somMass / rootSampleArea , digits = 4) )

#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_per_m2 <- dil_calc %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth) %>% 
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass) %>% rename(frag_g_m2 = dryMass) 
bbc_per_m2 <- merge(dil_per_m2, rootmass_per_m2, by = "sampleID")
  bbc_per_m2 <- bbc_per_m2 %>%  mutate(per_frag = round(100* frag_g_m2 / (frag_g_m2 + root_g_m2), digits = 1) )
bbc_per_m2_multiyear <- bbc_per_m2 %>% filter(year2 ==1 | year2 ==2) %>% filter(!is.na(rootSampleArea)) %>%  mutate(bbc_g_m2 = frag_g_m2 + root_g_m2) 
frag_year_glm <- bbc_per_m2_multiyear %>% group_by(siteID) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) ))
  write.table(frag_year_glm, paste0("year_frag_glm_",Sys.Date(),".txt"), sep = "\t", row.names=F)
bbc_year_glm <- bbc_per_m2_multiyear %>% group_by(siteID) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) ))
  write.table(bbc_year_glm, paste0("year_bbc_glm_",Sys.Date(),".txt"), sep = "\t", row.names=F)

bbc_per_m2_multiyear_for_plotting <- bbc_per_m2_multiyear %>% group_by(domainID, siteID, year2) %>% 
    summarise(bbc_sd = sd(bbc_g_m2), core_count = n(), bbc_se = sd(bbc_g_m2)/sqrt(core_count), bbc_g_m2 = mean(bbc_g_m2))

bbc_per_m2_multiyear_plotID_count <- bbc_per_m2_multiyear %>% group_by(domainID, siteID, plotID, year2) %>% summarise(plot_count = n())

#bbc_per_m2_multiyear_wide <- bbc_per_m2_multiyear  %>% pivot_wider(id_cols = c(domainID, siteID, plotID), names_from = year2, values_from = c(frag_g_m2) ) %>% 
#     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1, bbc_change_per = 100*bbc_change/bbc_g_m2_1)

## Hack to add asterisks for the 4 sites with significant (p < 0.05) bout effect
## If sites or statistical results change then would need to change.
dat_text <- data.frame(label = c("","","","","*","*","*","","","*","","","",""),
  siteID   = c("BART","GUAN","JORN","KONZ","MOAB","ONAQ","ORNL","SCBI","STEI","STER","TALL","TOOL","UNDE","WOOD") )

pdf(paste0("bbc_2_bouts_",Sys.Date(),".pdf"), width =8, height = 5)
ggplot(data = bbc_per_m2_multiyear_for_plotting, aes(x=year2, y=bbc_g_m2)) + geom_point(aes(colour=siteID))  +
  geom_errorbar(aes(ymin=bbc_g_m2-bbc_se, ymax=bbc_g_m2+bbc_se), width=.2,  position=position_dodge(.9))  +
   facet_wrap(~siteID, ncol = 5, scales = "free") +
   geom_text(data = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -3,   vjust   = -1, size =9)
#  annotate("text",x=-1,y=-3.1,label="Scatterplot Display")+coord_cartesian(ylim=c(-2.5,3),clip="off")
dev.off()

frag_site <- bbc_per_m2 %>% group_by(year, year2, domainID, siteID) %>% summarise(frag_g_m2_true = mean(frag_g_m2, na.rm = TRUE), root_g_m2_true = mean(root_g_m2, na.rm = TRUE)) %>% 
  mutate(bbc_g_m2_true = frag_g_m2_true + root_g_m2_true, per_frag = 100 * frag_g_m2_true / (frag_g_m2_true + root_g_m2_true), "bbc-10%" = bbc_g_m2_true * 0.9, "bbc+10%" = bbc_g_m2_true * 1.1, "frag-10%" = frag_g_m2_true * 0.9, "frag+10%" = frag_g_m2_true * 1.1)
frag_site_wide <- frag_site %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID), names_from = year2, values_from = c(frag_g_m2_true, bbc_g_m2_true) ) %>% 
     mutate(frag_change_true = frag_g_m2_true_2 - frag_g_m2_true_1, bbc_change_true = bbc_g_m2_true_2 - bbc_g_m2_true_1, bbc_change_per_true = 100*bbc_change_true/bbc_g_m2_true_1,
            "frag_change-10%" = frag_change_true - (abs(frag_change_true) * 0.1), "frag_change+10%" = frag_change_true + (abs(frag_change_true) * 0.1),
            "bbc_change-10%" = bbc_change_true - (abs(bbc_change_true) * 0.1), "bbc_change+10%" = bbc_change_true  + (abs(bbc_change_true) * 0.1),
            "bbc_change_per-5%" = bbc_change_per_true - (abs(bbc_change_per_true) * 0.05), "bbc_change_per+5%" = bbc_change_per_true  + (abs(bbc_change_per_true) * 0.05))
frag_per_summary <- frag_site %>% select(siteID, year, per_frag)

core_effort <- dil_per_m2 %>% group_by(domainID, siteID, year, year2) %>% tally() %>% rename("cores_per_site" = "n")
 cores_per_plot <- dil_per_m2 %>% group_by(domainID, siteID,plotID, year, year2) %>% tally()
plot_effort <- cores_per_plot %>% group_by(domainID, siteID, year, year2) %>% tally() %>% rename("plots_per_site" = "n")
effort <- merge(plot_effort, core_effort, by = c("domainID", "siteID","year","year2"))
 effort$year2 <- if_else(effort$year2 !="1" & effort$year2 !="2" , "2", effort$year2, effort$year2)
effort_wide <- effort  %>% pivot_wider(id_cols = c(domainID, siteID), names_from = year2, values_from = c(year, plots_per_site, cores_per_site) )
 View(effort_wide)
write.table(effort_wide, paste0("effort_",Sys.Date(),".txt"), sep = "\t", row.names=F)

groups <- unique(dil_calc$dilutionSampleID)
repetitions <- 1000

########## 10 ###############
dil_calc_10sub_rep = data.frame()
sample_set <- 10 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_10sub_rep = dplyr::bind_rows(dil_calc_10sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_10sub <- dil_calc_10sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_10sub <- merge(dil_var_10sub, rootmass_per_m2, by = "sampleID")
dil_var_10sub <- dil_var_10sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_10sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_10sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_10sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_10sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_10sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_10sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_10sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_10sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_10sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_10sub <- dil_var_10sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_10sub_wide <- frag_site_10sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1, bbc_change_per = 100*bbc_change/bbc_g_m2_1)
    frag_site_10sub_wide <- merge(frag_site_10sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_10sub_wide$fragFlag <- ifelse(frag_site_10sub_wide$frag_change > frag_site_10sub_wide$`frag_change-10%` & frag_site_10sub_wide$frag_change < frag_site_10sub_wide$`frag_change+10%`, 1, 0)
    frag_site_10sub_wide$bbcFlag <- ifelse(frag_site_10sub_wide$bbc_change > frag_site_10sub_wide$`bbc_change-10%` & frag_site_10sub_wide$bbc_change < frag_site_10sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_10sub_wide$bbcFlag_per <- ifelse(frag_site_10sub_wide$bbc_change_per > frag_site_10sub_wide$`bbc_change_per-5%` & frag_site_10sub_wide$bbc_change_per < frag_site_10sub_wide$`bbc_change_per+5%`, 1, 0)
  frag_site_10sub_wide_summary <- frag_site_10sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_10sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_10sub = round((frag_reps_inRange/reps_Number)*100, digits=1),
      bbc_reps_change_inRange = sum(bbcFlag_per), reps_Number = n(), bbc_per_change_inRange_10sub = round((bbc_reps_change_inRange/reps_Number)*100, digits=1)) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number, -bbc_reps_change_inRange) )
  
frag_site_10sub <- merge(frag_site_10sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_10sub$bbcFlag <- ifelse(frag_site_10sub$bbc_g_m2 > frag_site_10sub$`bbc-10%` & frag_site_10sub$bbc_g_m2 < frag_site_10sub$`bbc+10%`, 1, 0)
frag_site_10sub$fragFlag <- ifelse(frag_site_10sub$frag_g_m2 > frag_site_10sub$`frag-10%` & frag_site_10sub$frag_g_m2 < frag_site_10sub$`frag+10%`, 1, 0)

frag_site_10sub_summary <- frag_site_10sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_10sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_10sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 9 ###############
dil_calc_9sub_rep = data.frame()
sample_set <- 9 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_9sub_rep = dplyr::bind_rows(dil_calc_9sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_9sub <- dil_calc_9sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_9sub <- merge(dil_var_9sub, rootmass_per_m2, by = "sampleID")
dil_var_9sub <- dil_var_9sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_9sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_9sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_9sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_9sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_9sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_9sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_9sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_9sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_9sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_9sub <- dil_var_9sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_9sub_wide <- frag_site_9sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_9sub_wide <- merge(frag_site_9sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_9sub_wide$bbcFlag <- ifelse(frag_site_9sub_wide$bbc_change > frag_site_9sub_wide$`bbc_change-10%` & frag_site_9sub_wide$bbc_change < frag_site_9sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_9sub_wide$fragFlag <- ifelse(frag_site_9sub_wide$frag_change > frag_site_9sub_wide$`frag_change-10%` & frag_site_9sub_wide$frag_change < frag_site_9sub_wide$`frag_change+10%`, 1, 0)
  frag_site_9sub_wide_summary <- frag_site_9sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_9sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_9sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_9sub <- merge(frag_site_9sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_9sub$bbcFlag <- ifelse(frag_site_9sub$bbc_g_m2 > frag_site_9sub$`bbc-10%` & frag_site_9sub$bbc_g_m2 < frag_site_9sub$`bbc+10%`, 1, 0)
frag_site_9sub$fragFlag <- ifelse(frag_site_9sub$frag_g_m2 > frag_site_9sub$`frag-10%` & frag_site_9sub$frag_g_m2 < frag_site_9sub$`frag+10%`, 1, 0)

frag_site_9sub_summary <- frag_site_9sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_9sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_9sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 8 ###############
dil_calc_8sub_rep = data.frame()
sample_set <- 8 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_8sub_rep = dplyr::bind_rows(dil_calc_8sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_8sub <- dil_calc_8sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_8sub <- merge(dil_var_8sub, rootmass_per_m2, by = "sampleID")
dil_var_8sub <- dil_var_8sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_8sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_8sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_8sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_8sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_8sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_8sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_8sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_8sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_8sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_8sub <- dil_var_8sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_8sub_wide <- frag_site_8sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_8sub_wide <- merge(frag_site_8sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_8sub_wide$bbcFlag <- ifelse(frag_site_8sub_wide$bbc_change > frag_site_8sub_wide$`bbc_change-10%` & frag_site_8sub_wide$bbc_change < frag_site_8sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_8sub_wide$fragFlag <- ifelse(frag_site_8sub_wide$frag_change > frag_site_8sub_wide$`frag_change-10%` & frag_site_8sub_wide$frag_change < frag_site_8sub_wide$`frag_change+10%`, 1, 0)
  frag_site_8sub_wide_summary <- frag_site_8sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_8sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_8sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_8sub <- merge(frag_site_8sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_8sub$bbcFlag <- ifelse(frag_site_8sub$bbc_g_m2 > frag_site_8sub$`bbc-10%` & frag_site_8sub$bbc_g_m2 < frag_site_8sub$`bbc+10%`, 1, 0)
frag_site_8sub$fragFlag <- ifelse(frag_site_8sub$frag_g_m2 > frag_site_8sub$`frag-10%` & frag_site_8sub$frag_g_m2 < frag_site_8sub$`frag+10%`, 1, 0)

frag_site_8sub_summary <- frag_site_8sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_8sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_8sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 7 ###############
dil_calc_7sub_rep = data.frame()
sample_set <- 7 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_7sub_rep = dplyr::bind_rows(dil_calc_7sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_7sub <- dil_calc_7sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_7sub <- merge(dil_var_7sub, rootmass_per_m2, by = "sampleID")
dil_var_7sub <- dil_var_7sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_7sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_7sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_7sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_7sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_7sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_7sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_7sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_7sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_7sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_7sub <- dil_var_7sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_7sub_wide <- frag_site_7sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_7sub_wide <- merge(frag_site_7sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_7sub_wide$bbcFlag <- ifelse(frag_site_7sub_wide$bbc_change > frag_site_7sub_wide$`bbc_change-10%` & frag_site_7sub_wide$bbc_change < frag_site_7sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_7sub_wide$fragFlag <- ifelse(frag_site_7sub_wide$frag_change > frag_site_7sub_wide$`frag_change-10%` & frag_site_7sub_wide$frag_change < frag_site_7sub_wide$`frag_change+10%`, 1, 0)
  frag_site_7sub_wide_summary <- frag_site_7sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_7sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_7sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_7sub <- merge(frag_site_7sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_7sub$bbcFlag <- ifelse(frag_site_7sub$bbc_g_m2 > frag_site_7sub$`bbc-10%` & frag_site_7sub$bbc_g_m2 < frag_site_7sub$`bbc+10%`, 1, 0)
frag_site_7sub$fragFlag <- ifelse(frag_site_7sub$frag_g_m2 > frag_site_7sub$`frag-10%` & frag_site_7sub$frag_g_m2 < frag_site_7sub$`frag+10%`, 1, 0)

frag_site_7sub_summary <- frag_site_7sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_7sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_7sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 6 ###############
dil_calc_6sub_rep = data.frame()
sample_set <- 6 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_6sub_rep = dplyr::bind_rows(dil_calc_6sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_6sub <- dil_calc_6sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_6sub <- merge(dil_var_6sub, rootmass_per_m2, by = "sampleID")
dil_var_6sub <- dil_var_6sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_6sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_6sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_6sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_6sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_6sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_6sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_6sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_6sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_6sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_6sub <- dil_var_6sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_6sub_wide <- frag_site_6sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_6sub_wide <- merge(frag_site_6sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_6sub_wide$bbcFlag <- ifelse(frag_site_6sub_wide$bbc_change > frag_site_6sub_wide$`bbc_change-10%` & frag_site_6sub_wide$bbc_change < frag_site_6sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_6sub_wide$fragFlag <- ifelse(frag_site_6sub_wide$frag_change > frag_site_6sub_wide$`frag_change-10%` & frag_site_6sub_wide$frag_change < frag_site_6sub_wide$`frag_change+10%`, 1, 0)
  frag_site_6sub_wide_summary <- frag_site_6sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_6sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_6sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_6sub <- merge(frag_site_6sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_6sub$bbcFlag <- ifelse(frag_site_6sub$bbc_g_m2 > frag_site_6sub$`bbc-10%` & frag_site_6sub$bbc_g_m2 < frag_site_6sub$`bbc+10%`, 1, 0)
frag_site_6sub$fragFlag <- ifelse(frag_site_6sub$frag_g_m2 > frag_site_6sub$`frag-10%` & frag_site_6sub$frag_g_m2 < frag_site_6sub$`frag+10%`, 1, 0)

frag_site_6sub_summary <- frag_site_6sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_6sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_6sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 5 ###############
dil_calc_5sub_rep = data.frame()
sample_set <- 5 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_5sub_rep = dplyr::bind_rows(dil_calc_5sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_5sub <- dil_calc_5sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_5sub <- merge(dil_var_5sub, rootmass_per_m2, by = "sampleID")
dil_var_5sub <- dil_var_5sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_5sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_5sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_5sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_5sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_5sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_5sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_5sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_5sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_5sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_5sub <- dil_var_5sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_5sub_wide <- frag_site_5sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_5sub_wide <- merge(frag_site_5sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_5sub_wide$bbcFlag <- ifelse(frag_site_5sub_wide$bbc_change > frag_site_5sub_wide$`bbc_change-10%` & frag_site_5sub_wide$bbc_change < frag_site_5sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_5sub_wide$fragFlag <- ifelse(frag_site_5sub_wide$frag_change > frag_site_5sub_wide$`frag_change-10%` & frag_site_5sub_wide$frag_change < frag_site_5sub_wide$`frag_change+10%`, 1, 0)
  frag_site_5sub_wide_summary <- frag_site_5sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_5sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_5sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_5sub <- merge(frag_site_5sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_5sub$bbcFlag <- ifelse(frag_site_5sub$bbc_g_m2 > frag_site_5sub$`bbc-10%` & frag_site_5sub$bbc_g_m2 < frag_site_5sub$`bbc+10%`, 1, 0)
frag_site_5sub$fragFlag <- ifelse(frag_site_5sub$frag_g_m2 > frag_site_5sub$`frag-10%` & frag_site_5sub$frag_g_m2 < frag_site_5sub$`frag+10%`, 1, 0)

frag_site_5sub_summary <- frag_site_5sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_5sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_5sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 4 ###############
dil_calc_4sub_rep = data.frame()
sample_set <- 4 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_4sub_rep = dplyr::bind_rows(dil_calc_4sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_4sub <- dil_calc_4sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_4sub <- merge(dil_var_4sub, rootmass_per_m2, by = "sampleID")
dil_var_4sub <- dil_var_4sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_4sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_4sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_4sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_4sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_4sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_4sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_4sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_4sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_4sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_4sub <- dil_var_4sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_4sub_wide <- frag_site_4sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_4sub_wide <- merge(frag_site_4sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_4sub_wide$bbcFlag <- ifelse(frag_site_4sub_wide$bbc_change > frag_site_4sub_wide$`bbc_change-10%` & frag_site_4sub_wide$bbc_change < frag_site_4sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_4sub_wide$fragFlag <- ifelse(frag_site_4sub_wide$frag_change > frag_site_4sub_wide$`frag_change-10%` & frag_site_4sub_wide$frag_change < frag_site_4sub_wide$`frag_change+10%`, 1, 0)
  frag_site_4sub_wide_summary <- frag_site_4sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_4sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_4sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_4sub <- merge(frag_site_4sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_4sub$bbcFlag <- ifelse(frag_site_4sub$bbc_g_m2 > frag_site_4sub$`bbc-10%` & frag_site_4sub$bbc_g_m2 < frag_site_4sub$`bbc+10%`, 1, 0)
frag_site_4sub$fragFlag <- ifelse(frag_site_4sub$frag_g_m2 > frag_site_4sub$`frag-10%` & frag_site_4sub$frag_g_m2 < frag_site_4sub$`frag+10%`, 1, 0)

frag_site_4sub_summary <- frag_site_4sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_4sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_4sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

########## 3 ###############
dil_calc_3sub_rep = data.frame()
sample_set <- 3 # full set is 10 dilutions, try subsets of 2, 4, 6, and 8
for(i in 1:length(groups)){
reps <- dil_calc  %>% filter(dilutionSampleID == groups[i]) %>% rep_slice_sample(n = sample_set, replace=TRUE, reps = repetitions)
reps$sample<-rep(c(1:sample_set),times=repetitions)
dil_calc_3sub_rep = dplyr::bind_rows(dil_calc_3sub_rep, reps)
}

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub); rm(frag_year_sub)
#   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
dil_var_3sub <- dil_calc_3sub_rep %>%
  group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth, replicate) %>%
  summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
            somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass)  %>% rename(frag_g_m2 = dryMass)

dil_var_3sub <- merge(dil_var_3sub, rootmass_per_m2, by = "sampleID")
dil_var_3sub <- dil_var_3sub %>% mutate(bbc_g_m2 = frag_g_m2 + root_g_m2)
bbc_per_m2_multiyear_sub <- dil_var_3sub %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(bbc_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_3sub_summary <- bbc_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_3sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_3sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_3sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(domainID,siteID,replicate) %>% do(tidy(glm(frag_g_m2 ~ year2, data = .) )) %>% filter(term == "year22")
frag_year_sub$yearFlag <- ifelse(frag_year_sub$p.value < 0.05, 1, 0)
frag_year_sub$yearDir <- ifelse(frag_year_sub$estimate > 0, 1, 0)
frag_year_3sub_summary <- frag_year_sub %>% group_by(domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_3sub = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_3sub = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_3sub = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )

frag_site_3sub <- dil_var_3sub %>% group_by(year, year2, domainID, siteID, replicate) %>% summarise(bbc_g_m2 = mean(bbc_g_m2, na.rm = TRUE), frag_g_m2 = mean(frag_g_m2, na.rm = TRUE)) 

  frag_site_3sub_wide <- frag_site_3sub %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(domainID, siteID, replicate), names_from = year2, values_from = c(frag_g_m2, bbc_g_m2) ) %>% 
     mutate(frag_change = frag_g_m2_2 - frag_g_m2_1, bbc_change = bbc_g_m2_2 - bbc_g_m2_1)
    frag_site_3sub_wide <- merge(frag_site_3sub_wide, frag_site_wide, by = c("domainID", "siteID"))
    frag_site_3sub_wide$bbcFlag <- ifelse(frag_site_3sub_wide$bbc_change > frag_site_3sub_wide$`bbc_change-10%` & frag_site_3sub_wide$bbc_change < frag_site_3sub_wide$`bbc_change+10%`, 1, 0)
    frag_site_3sub_wide$fragFlag <- ifelse(frag_site_3sub_wide$frag_change > frag_site_3sub_wide$`frag_change-10%` & frag_site_3sub_wide$frag_change < frag_site_3sub_wide$`frag_change+10%`, 1, 0)
  frag_site_3sub_wide_summary <- frag_site_3sub_wide %>% group_by(domainID, siteID) %>% 
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_3sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_3sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
      select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )

frag_site_3sub <- merge(frag_site_3sub, frag_site, by = c("year", "year2", "domainID", "siteID"))
frag_site_3sub$bbcFlag <- ifelse(frag_site_3sub$bbc_g_m2 > frag_site_3sub$`bbc-10%` & frag_site_3sub$bbc_g_m2 < frag_site_3sub$`bbc+10%`, 1, 0)
frag_site_3sub$fragFlag <- ifelse(frag_site_3sub$frag_g_m2 > frag_site_3sub$`frag-10%` & frag_site_3sub$frag_g_m2 < frag_site_3sub$`frag+10%`, 1, 0)

frag_site_3sub_summary <- frag_site_3sub %>% group_by(year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_3sub = round((bbc_reps_inRange/reps_Number)*100, digits=1),
            frag_reps_inRange = sum(fragFlag), reps_Number = n(), frag_per_inRange_3sub = round((frag_reps_inRange/reps_Number)*100, digits=1) ) %>%
  select(c(-bbc_reps_inRange, -frag_reps_inRange, -reps_Number) )



save(dil_calc_10sub_rep, dil_calc_9sub_rep, dil_calc_8sub_rep, dil_calc_7sub_rep, dil_calc_6sub_rep, dil_calc_5sub_rep, 
     dil_calc_4sub_rep, dil_calc_3sub_rep, file = "dilReps.RData")
#load("dilReps.RData")


############# merge the different dilution subsets #################################################

df_list <- list(frag_site_10sub_summary, frag_site_9sub_summary, frag_site_8sub_summary, frag_site_7sub_summary, frag_site_6sub_summary, 
                           frag_site_5sub_summary, frag_site_4sub_summary, frag_site_3sub_summary)      
frag_site_summary <- df_list %>% reduce(full_join, by=c("year", "year2", "domainID", "siteID") )
 View(frag_site_summary)
frag_site_summary <- merge(frag_site_summary, frag_per_summary, by = c("year", "year2", "domainID", "siteID") )
write.table(frag_site_summary, paste0("frag_site_summary_",Sys.Date(),".txt"), sep = "\t", row.names=F)


bbc_year_list <- list(bbc_year_10sub_summary, bbc_year_9sub_summary, bbc_year_8sub_summary, bbc_year_7sub_summary, bbc_year_6sub_summary, 
                           bbc_year_5sub_summary, bbc_year_4sub_summary, bbc_year_3sub_summary)      
bbc_year_summary <- bbc_year_list %>% reduce(full_join, by=c("domainID", "siteID") )
 View(bbc_year_summary)
write.table(bbc_year_summary, paste0("year_summary_bbc_",Sys.Date(),".txt"), sep = "\t", row.names=F)

bbc_year_summary_subset <- bbc_year_summary %>% select(domainID, siteID, estimate_10sub, year_per_sig_10sub, 
                     estimate_7sub, year_per_sig_7sub, estimate_3sub, year_per_sig_3sub) %>% 
    filter(siteID == "ORNL" | siteID == "STER" | siteID == "MOAB" | siteID == "ONAQ")
View(bbc_year_summary_subset)
write.table(bbc_year_summary_subset, paste0("year_summary_bbc_subset_",Sys.Date(),".txt"), sep = "\t", row.names=F)



frag_year_list <- list(frag_year_10sub_summary, frag_year_9sub_summary, frag_year_8sub_summary, frag_year_7sub_summary, frag_year_6sub_summary, 
                           frag_year_5sub_summary, frag_year_4sub_summary, frag_year_3sub_summary)      
frag_year_summary <- frag_year_list %>% reduce(full_join, by=c("domainID", "siteID") )
 View(frag_year_summary)
write.table(frag_year_summary, paste0("year_summary_frag_",Sys.Date(),".txt"), sep = "\t", row.names=F)

frag_year_summary_subset <- frag_year_summary %>% select(domainID, siteID, estimate_10sub, year_per_sig_10sub, 
                     estimate_7sub, year_per_sig_7sub, estimate_3sub, year_per_sig_3sub) %>% 
    filter(siteID == "WOOD" | siteID == "STER" | siteID == "ONAQ" | siteID == "GUAN" | siteID == "ORNL" | siteID == "JORN" | siteID == "UNDE")
View(frag_year_summary_subset)
frag_year_summary_subset <- frag_year_summary_subset %>% arrange(year_per_sig_7sub,year_per_sig_3sub)
write.table(frag_year_summary_subset, paste0("year_summary_frag_subset_",Sys.Date(),".txt"), sep = "\t", row.names=F)

############ look at frags PLUS long roots #######################
bbc_site_summary <- frag_site_summary %>% select(year, year2, domainID, siteID, contains("bbc_"))
#View(bbc_site_summary)
write.table(bbc_site_summary, paste0("bbc_lt2mm_site_summary_",Sys.Date(),".txt"), sep = "\t", row.names=F)

bbc_site_summary_long <- bbc_site_summary %>% 
  pivot_longer(
    cols = !c(year, year2, domainID, siteID), 
    names_to = "subsamples", 
    values_to = "percent"
  )
bbc_site_summary_long$subsamples <- gsub("bbc_per_inRange_", "", bbc_site_summary_long$subsamples)
bbc_site_summary_long$subsamples <- as.numeric (gsub("sub", "", bbc_site_summary_long$subsamples) )
#write.table(bbc_site_summary_long, paste0("bbc_lt2mm_site_summary_long_",Sys.Date(),".txt"), sep = "\t", row.names=F)

pdf(paste0("bbc_lt2mm_site_facets_",Sys.Date(),".pdf"), width =8, height = 5)
#png("bbc_lt2mm_site_facets.png")
ggplot(data = bbc_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(colour=siteID, linetype = year2)) + 
  scale_linetype_manual(values=c("twodash","solid","solid","solid","solid","solid","solid")) + scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) +  ylim(75, 100) + 
  facet_wrap(~siteID, ncol = 10) + 
  labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
    theme_minimal() + theme(legend.position="none")
dev.off()

pdf(paste0("bbc_lt2mm_sites_combined_",Sys.Date(),".pdf"), width =8, height = 5)
#png("bbc_lt2mm_sites_combined.png")
ggplot(data = bbc_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(colour=siteID, linetype = year2))  + 
  scale_linetype_manual(values=c("twodash","solid","solid","solid","solid","solid","solid")) +  scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) + ylim(75, 100) +
   labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
    theme_minimal() + guides(linetype = "none")
dev.off()

############ look at ONLY frags #######################
frag_only_site_summary <- frag_site_summary %>% select(year, year2, domainID, siteID, contains("frag_"))
View(frag_only_site_summary)
write.table(frag_only_site_summary, paste0("frag_only_site_summary_",Sys.Date(),".txt"), sep = "\t", row.names=F)

frag_only_site_summary_long <- frag_only_site_summary %>% 
  pivot_longer(
    cols = !c(year, year2, domainID, siteID), 
    names_to = "subsamples", 
    values_to = "percent"
  )
frag_only_site_summary_long$subsamples <- gsub("frag_per_inRange_", "", frag_only_site_summary_long$subsamples)
frag_only_site_summary_long$subsamples <- as.numeric (gsub("sub", "", frag_only_site_summary_long$subsamples) )
View(frag_only_site_summary_long)
write.table(frag_only_site_summary_long, paste0("frag_only_site_summary_long_",Sys.Date(),".txt"), sep = "\t", row.names=F)

pdf(paste0("frag_only_site_facets_",Sys.Date(),".pdf"), width =8, height = 5)
#png("frag_only_site_facets.png")
ggplot(data = frag_only_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(colour=siteID, linetype = year2))  + 
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid")) + scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) + facet_wrap(~siteID, ncol = 10) + 
   labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
   theme_minimal() + theme(legend.position="none")
dev.off()

pdf(paste0("frag_only_sites_combined_",Sys.Date(),".pdf"), width =8, height = 5)
#png("frag_only_sites_combined.png")
ggplot(data = frag_only_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(colour=siteID, linetype = year2)) + 
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid"))+ scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) + 
   labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))   + 
   theme_minimal() + guides(linetype = "none") # + scale_linetype_discrete(guide = "none")
dev.off()



