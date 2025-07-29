
############ look at frags PLUS long roots #######################
#bbc_site_summary <- frag_site_summary %>% select(year, year2, domainID, siteID, contains("bbc_"))
#View(bbc_site_summary)
#write.table(bbc_site_summary, paste0("bbc_lt2mm_site_summary_",Sys.Date(),".txt"), sep = "\t", row.names=F)

bbc_site_summary <- read.delim("bbc_lt2mm_site_summary_2023-05-03.txt", header=T, sep = "\t", stringsAsFactors = F)
bbc_site_summary$year <- as.character(bbc_site_summary$year)
bbc_site_summary$year2 <- as.character(bbc_site_summary$year2)

bbc_site_summary_long <- bbc_site_summary %>% 
  pivot_longer(
    cols = !c(year, year2, domainID, siteID), 
    names_to = "subsamples", 
    values_to = "percent"
  )
bbc_site_summary_long$subsamples <- gsub("bbc_per_inRange_", "", bbc_site_summary_long$subsamples)
bbc_site_summary_long$subsamples <- as.numeric (gsub("sub", "", bbc_site_summary_long$subsamples) )


pdf(paste0("bbc_lt2mm_site_facets_",Sys.Date(),".pdf"), width =8, height = 5)
#png("bbc_lt2mm_site_facets.png")
ggplot(data = bbc_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(linetype = year2)) + 
  scale_linetype_manual(values=c("twodash","solid","solid","solid","solid","solid","solid")) + scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) +  ylim(85, 100) + 
  facet_wrap(~siteID, ncol = 10) + 
  labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations \n within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
    theme_minimal() + theme(legend.position="none") +
  theme(axis.title = element_text(size = 14))  
dev.off()

pdf(paste0("bbc_lt2mm_sites_combined_",Sys.Date(),".pdf"), width =8, height = 5)
#png("bbc_lt2mm_sites_combined.png")
ggplot(data = bbc_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(colour=siteID, linetype = year2))  + 
  scale_linetype_manual(values=c("twodash","solid","solid","solid","solid","solid","solid")) +  scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) + ylim(85, 100) +
   labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations \n within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
    theme_minimal() + guides(linetype = "none") +
  theme(axis.title = element_text(size = 14))  
dev.off()

############ look at ONLY frags #######################
# frag_only_site_summary <- frag_site_summary %>% select(year, year2, domainID, siteID, contains("frag_"))
# View(frag_only_site_summary)
#write.table(frag_only_site_summary, paste0("frag_only_site_summary_",Sys.Date(),".txt"), sep = "\t", row.names=F)
frag_only_site_summary <- read.delim("frag_only_site_summary_2023-05-03.txt", header=T, sep = "\t", stringsAsFactors = F)
frag_only_site_summary$year <- as.character(frag_only_site_summary$year)
frag_only_site_summary$year2 <- as.character(frag_only_site_summary$year2)

frag_only_site_summary_long <- frag_only_site_summary %>% 
  pivot_longer(
    cols = !c(year, year2, domainID, siteID), 
    names_to = "subsamples", 
    values_to = "percent"
  )
frag_only_site_summary_long$subsamples <- gsub("frag_per_inRange_", "", frag_only_site_summary_long$subsamples)
frag_only_site_summary_long$subsamples <- as.numeric (gsub("sub", "", frag_only_site_summary_long$subsamples) )
View(frag_only_site_summary_long)


pdf(paste0("frag_only_site_facets_",Sys.Date(),".pdf"), width =8, height = 5)
#png("frag_only_site_facets.png")
ggplot(data = frag_only_site_summary_long, aes(x=subsamples, y=percent)) + 
  geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) + facet_wrap(~siteID, ncol = 10) + 
  geom_line(aes(linetype = year2))  + 
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid")) + scale_x_reverse() + 
   labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations \n within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
   theme_minimal() + theme(legend.position="none") +
  theme(axis.title = element_text(size = 14))  
dev.off()

pdf(paste0("frag_only_sites_combined_",Sys.Date(),".pdf"), width =8, height = 5)
#png("frag_only_sites_combined.png")
ggplot(data = frag_only_site_summary_long, aes(x=subsamples, y=percent)) + geom_line(aes(colour=siteID, linetype = year2)) + 
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid"))+ scale_x_reverse() + geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) + 
   labs(x = "Dilution subsamples") + labs(y = "Percent of resample iterations \n within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))   + 
   theme_minimal() + guides(linetype = "none")  +
  theme(axis.title = element_text(size = 14))  # + scale_linetype_discrete(guide = "none")
dev.off()



## multi-year effect (don't have the needed .txt files so tring to reproduce them here)
library(neonUtilities)
library(tidyverse)
library(broom) # includes broom function
library(infer) # includes rep_slice_sample function
library(ggplot2)

# 
# bbcDat <- loadByProduct(dpID="DP1.10067.001", enddate = "2022-12", package = "basic", check.size = FALSE, token = Sys.getenv('NEON_PAT')) # pulled from portal 2025-06-25
#  list2env(bbcDat ,.GlobalEnv) # unlist all data frames  
# 
# bbc_rootmass_lt2 <- bbc_rootmass %>%  filter(sizeCategory != "2-10")
# rootmass_by_sampleID <- bbc_rootmass_lt2 %>% group_by(sampleID) %>% summarise(root_dryMass = sum(dryMass))
#  
# bbc_percore <- bbc_percore %>% select(subplotID, clipID, coreID, sampleID, wst10cmDist, wst1cmDist, litterDepth, rootSampleArea, rootSampleDepth)
# 
# rootmass_per_m2 <- merge(rootmass_by_sampleID, bbc_percore, by = "sampleID", all.x = TRUE)
# rootmass_per_m2 <- rootmass_per_m2 %>% mutate(root_g_m2 = round(root_dryMass / rootSampleArea , digits = 4) ) %>% select(sampleID, root_dryMass, rootSampleArea, root_g_m2)
# 
# bbc_dilution <- merge(bbc_dilution, bbc_percore, by = "sampleID", all.x = TRUE)
# bbc_dilution <- bbc_dilution %>% filter(!is.na(rootSampleArea))
# 
# bbc_dilution$year <- substr(bbc_dilution$collectDate,1,4)
# bbc_dilution$eventID <- paste0(bbc_dilution$siteID, "_", bbc_dilution$year)
# bbc_dil_unique_eventID <- bbc_dilution %>% select(siteID, year, eventID) %>% distinct()  
# bbc_dilution_yr_count <- bbc_dil_unique_eventID %>% group_by(siteID) %>% tally()
#  
# bbc_dilution$year <- substr(bbc_dilution$collectDate,1,4)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "BART" & bbc_dilution$year == "2016", 1, bbc_dilution$year)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "BART" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "GUAN" & bbc_dilution$year == "2018", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "GUAN" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "JORN" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "JORN" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "KONZ" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "KONZ" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "MOAB" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "MOAB" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ONAQ" & bbc_dilution$year == "2016", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ONAQ" & bbc_dilution$year == "2021", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ORNL" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "ORNL" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "SCBI" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "SCBI" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STEI" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STEI" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STER" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "STER" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TALL" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TALL" & bbc_dilution$year == "2021", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TOOL" & bbc_dilution$year == "2017", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "TOOL" & bbc_dilution$year == "2022", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "UNDE" & bbc_dilution$year == "2016", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "UNDE" & bbc_dilution$year == "2019", 2, bbc_dilution$year2)
# bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "WOOD" & bbc_dilution$year == "2016", 1, bbc_dilution$year2)
#  bbc_dilution$year2 <- ifelse(bbc_dilution$siteID == "WOOD" & bbc_dilution$year == "2021", 2, bbc_dilution$year2)  
# 
# # look at all 47 sites and for the 14 sites with 2 bouts look at both bouts
# # filter to remove outliers > 99 percentile for whole dataset and filter to remove fragMass < 0 g
# bbc_dilution <- bbc_dilution %>% filter(dryMass >= 0) %>%  filter(dryMass < quantile(dryMass, probs = 0.99, na.rm =TRUE))
# 
# # if there are < 10 subsamples per core then drop the entire core
# #dil_is_10 <- bbc_dilution %>% group_by(sampleID) %>% summarise(nDM = n()) %>%  filter(nDM == 10)
# #  bbc_dilution <- merge(bbc_dilution, dil_is_10, all.y=TRUE)
# dil_is_gt7 <- bbc_dilution %>% group_by(sampleID) %>% summarise(nDM = n()) %>%  filter(nDM > 7)
#   bbc_dilution <- merge(bbc_dilution, dil_is_gt7, all.y=TRUE)
#   bbc_dilution$nDM <- NULL
# 
# ### Calculate dryMass for fragments < 1 cm length for each dilutionSubsampleID (calculated mass per unit area but could do on a soil volume basis instead)
# dil_calc <- bbc_dilution %>%
#   mutate(fragMass = round(dryMass*(sampleVolume/dilutionSubsampleVolume), digits = 4), frag_g_m2 = round(fragMass / rootSampleArea , digits = 4),
#          somMass = round(somDryMass*(sampleVolume/dilutionSubsampleVolume), digits = 4),  som_g_m2 = round(somMass / rootSampleArea , digits = 4) )
# 
# #   Calculate mean and variability of mass of SOM and fragments < 1 cm length for dilutionSubsampleIDs within each sampleID
# dil_per_m2 <- dil_calc %>%
#   group_by(year, year2, domainID, siteID, plotID, sampleID, wst10cmDist, wst1cmDist, litterDepth) %>% 
#   summarise(dryMass = mean(frag_g_m2), sdDM = sd(frag_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass,
#             somDryMass = mean(som_g_m2), sdSOM = sd(som_g_m2), nSOM = n(), seSOM = sdSOM/sqrt(nSOM), relseSOM = 100 * seSOM / somDryMass) %>% rename(frag_g_m2 = dryMass) 
# bbc_per_m2 <- merge(dil_per_m2, rootmass_per_m2, by = "sampleID")
#   bbc_per_m2 <- bbc_per_m2 %>%  mutate(per_frag = round(100* frag_g_m2 / (frag_g_m2 + root_g_m2), digits = 1) )
# bbc_per_m2_multiyear <- bbc_per_m2 %>% filter(year2 ==1 | year2 ==2) %>% filter(!is.na(rootSampleArea)) %>%  mutate(bbc_g_m2 = frag_g_m2 + root_g_m2) 
#  write.table(bbc_per_m2_multiyear, paste0("bbc_per_m2_multiyear_",Sys.Date(),".txt"), sep = "\t", row.names=F)
bbc_per_m2_multiyear <- read.delim("bbc_per_m2_multiyear_2025-06-25.txt", header=T, sep = "\t", stringsAsFactors = F)

bbc_per_m2_multiyear_for_plotting <- bbc_per_m2_multiyear %>% group_by(domainID, siteID, year2) %>% 
    summarise(bbc_sd = sd(bbc_g_m2), core_count = n(), bbc_se = sd(bbc_g_m2)/sqrt(core_count), bbc_g_m2 = mean(bbc_g_m2))
bbc_per_m2_multiyear_for_plotting$year2 <- as.character(bbc_per_m2_multiyear_for_plotting$year2)

## Hack to add asterisks for the 4 sites with significant (p < 0.05) bout effect
## If sites or statistical results change then would need to change.
dat_text <- data.frame(label = c("","","","","*","*","*","","","*","","","",""),
  siteID   = c("BART","GUAN","JORN","KONZ","MOAB","ONAQ","ORNL","SCBI","STEI","STER","TALL","TOOL","UNDE","WOOD") )

pdf(paste0("bbc_2_bouts_",Sys.Date(),".pdf"), width =8, height = 5)
ggplot(data = bbc_per_m2_multiyear_for_plotting, aes(x=year2, y=bbc_g_m2)) + geom_point(aes())  +
  geom_errorbar(aes(ymin=bbc_g_m2-bbc_se, ymax=bbc_g_m2+bbc_se), width=.2,  position=position_dodge(.9))  +
   facet_wrap(~siteID, ncol = 5, scales = "free") +
  xlab("Bout") +   
  ylab(bquote(Roots ~ (g/m^2))) +   
  geom_text(data = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -3,   vjust   = -1, size =9) +
  theme(axis.title = element_text(size = 14))  
dev.off()



