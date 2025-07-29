library(tidyverse)
library(lme4)
library(performance) #icc partitions variance for random effects
library(broom) # needed for tidy function
library(infer) # includes rep_slice_sample function


data_merge <- data_merge_lumped
#stature_year <- read.delim("stature_year.txt", header = TRUE, sep = "\t")
#data_merge <- read.delim("data_merge_lumped.txt", header = TRUE, sep = "\t")

stature_class <- data_merge %>% distinct(siteID, stature)

data_merge <- data_merge %>% filter(core_count == 2) # removes all of BARR, HEAL, YELL, plus some records in other sites
data_merge <- data_merge %>% filter(!(clip_count == 1 & stature == "large")) # removes ALL of STER and PUUM, plus some records in other sites

data_merge <- arrange(data_merge, eventID) 
boutNumber <- data_merge %>% select(siteID, eventID) %>% distinct() %>%
    dplyr::group_by(siteID) %>% 
    dplyr::mutate(year2 = row_number()) %>% ungroup() %>%
    dplyr::select(-siteID)
data_merge <- data_merge %>% merge(boutNumber, by = "eventID", all.x = TRUE) 

site_counts <- data_merge %>% select(siteID, year2, plot_count,  nlcd_count) %>% distinct()
 site_counts$nlcd_multi <- ifelse(site_counts$nlcd_count > 1, ">1 nlcd", "1 nlcd")


options(stringsAsFactors = FALSE)

###################################################
#check for normality
data_merge <- data_merge %>%
  mutate(response_log = log1p(response))


windows(9,5)
hist(data_merge$response_log)
site_list <- unique(data_merge$siteID)
site_list

size_list <- unique(data_merge$sizeCategory)
size_list


################### Variance partitioning #########################################################################
###################################################################################################################
# variable list for variance partitioning
my_varpart_list <- c('clipID','plotID','nlcdClass')

y_var <- 'response_log'

data_in <- data_merge %>% filter(year2 == 1) # %>% filter(sizeCategory == "total") ########## for analyses just focusing on the more complete first year of sampling per site
data_in$year <- as.numeric(data_in$year)

#hist(data_in$response)
#hist(data_in$response_log)

###################################################
# modify random vars

# change to factors for mixed modeling
for(i in my_varpart_list){
  data_in[,i] <- as_factor(data_in[,i])
}

###################################################
data_single_var_r2 <- data.frame()
for(site_id in (site_list)){
 for(size_id in (size_list)){
  for(x_var in my_varpart_list){
    mod_result_i <- NULL
    data_single_var_r2_i <- data.frame(
      siteID = site_id,
      sizeCategory = size_id,
      y_var = y_var,
      x_var = x_var,
       r2_cond = NA,
      anova_pval = NA)
    
    data_i <- data_in %>% filter(siteID == site_id, sizeCategory == size_id)
    
    # random effects model
    form_i <- paste0(y_var, " ~ 1 + (1|", x_var,")")
    try({mod_result_i <- lme4::lmer(formula = form_i, data = data_i)})
    
    if(!is.null(mod_result_i) &&
       !performance::check_singularity(mod_result_i)){
      r2_i <- performance::r2(mod_result_i)
      data_single_var_r2_i$r2_cond <- r2_i[[1]]["Conditional R2"]
    }
    
    # ANOVA
    form_i_aov <- paste0(y_var, " ~ ", x_var)
    try({mod_result_aov_i <- anova(lm(formula = form_i_aov, data = data_i))})
    
    if(!is.null(mod_result_aov_i)){
      data_single_var_r2_i$anova_pval <- round(
        mod_result_aov_i$`Pr(>F)`[1],
        5)
    }
    
    # add to table
    data_single_var_r2 <- data_single_var_r2 %>%
      bind_rows(data_single_var_r2_i)
    
    message(site_id, " | ", x_var)
  }
 }
}

site_counts_mod <- site_counts %>% filter(year2 == 1)

data_single_var_r2 <- data_single_var_r2 %>% left_join(site_counts_mod, by = "siteID")

data_single_var_r2 <- data_single_var_r2 %>% left_join (stature_year, by = "siteID") %>% filter(!(nlcd_count == 1 & x_var == "nlcdClass")) 

# manual removal of variables that have variance near zero in mixed model that causes model to fail (same when there is no filtering based on number of cores or clips AND where there is filtering based on number of cores)
data_single_var_r2 <- data_single_var_r2 %>% filter(!((siteID == "BONA" | siteID == "DCFS" | siteID == "DELA" | siteID == "GRSM" | siteID == "HARV" | siteID == "JERC" | siteID == "JORN" | 
    siteID == "KONZ" | siteID == "LENO" | siteID == "MOAB" | siteID == "NOGP" | siteID == "OAES" | siteID == "ORNL" | siteID == "SCBI" | siteID == "SJER" | siteID == "TALL" | siteID == "UNDE") & x_var == "nlcdClass"))
data_single_var_r2 <- data_single_var_r2 %>% filter(!((siteID == "DELA" | siteID == "DEJU" | siteID == "MLBS"  | siteID == "NIWO" | siteID == "OSBS"  | siteID == "SCBI" | 
    siteID == "SOAP" | siteID == "SRER" | siteID == "TEAK") & x_var == "clipID"))
data_single_var_r2 <- data_single_var_r2 %>% filter(!((siteID == "BONA" | siteID == "LENO" | siteID == "NIWO" |siteID == "SJER" |  siteID == "SOAP" | siteID == "STEI"  | siteID == "WREF") & x_var == "plotID"))

data_single_var_r2 <- data_single_var_r2 %>% filter(!((stature == "small") & x_var == "plotID"))

######################################################
# create formulas based off of prelim analysis above

site_formulas <- data.frame()
for(j in 1:length(size_list)){
for(i in 1:length(site_list)){
  site_i <- site_list[i]
  size_j <- size_list[j]
  
  # filter by significant anova results
  data_i <- data_single_var_r2 %>%
     filter(siteID == site_i, sizeCategory == size_j #,
#            anova_pval < 0.10,
#            !is.na(r2_cond)
)
  
  y_var_i <- data_i$y_var[1]
  
  if(nrow(data_i) > 0){
    form_i <- paste0(
      data_i$y_var[1],
      " ~ ",
      paste("(1|",data_i$x_var,")", sep = "", collapse = " + "))
  }else{
    form_i = NA_character_
  }
  
  site_formulas_i <- data.frame(
    siteID = site_i,
    sizeCategory = size_j,
    y_var = y_var_i,
    form = form_i)
  
  site_formulas <- bind_rows(
    site_formulas, site_formulas_i)
}
}

site_formulas$form <- ifelse(grepl("nlcd", site_formulas$form), gsub('1\\|plotID','1\\|plotID:nlcdClass',site_formulas$form), site_formulas$form)
site_formulas$form <- ifelse(grepl("nlcd", site_formulas$form), gsub('1\\|clipID','1\\|clipID:plotID:nlcdClass',site_formulas$form), site_formulas$form)
site_formulas$form <- ifelse(!grepl("nlcd", site_formulas$form) & grepl("plotID", site_formulas$form), gsub('1\\|clipID','1\\|clipID:plotID',site_formulas$form), site_formulas$form)
site_formulas$form <- ifelse(!grepl("nlcd", site_formulas$form) & !grepl("plotID", site_formulas$form), gsub('1\\|clipID','1\\|clipID',site_formulas$form), site_formulas$form)

write_delim(site_formulas, paste0("site_formulas.txt"), delim = "\t")


##############################################
# nest data for model fitting
################################

data_nested <- data_in %>%
  group_by(siteID, sizeCategory) %>%
  nest() %>%
  left_join(site_formulas, by = c("siteID","sizeCategory"))


# get counts of unique vals for each random factor
for(i in my_varpart_list){
  my_col_name <- paste0("n_",i)
  data_nested <- data_nested %>%
    mutate(
      !!my_col_name := map_int(data, ~length(unique(unlist(select(.,!!i))))))
}

# get total number of observations or samples at each site
# get a formula for each site
data_nested <- data_nested %>% # filter(sizeCategory == "total") %>% 
  mutate(
    n_obs = map_int(data, ~nrow(.)), n_per_plot = n_obs/n_plotID) #, clips_per_plot = n_clipID / n_plotID)
  

#############################################
# model fitting function

fn_fit_models <- function(data_x, form_x) {
  cat(as.data.frame(data_x)$siteID[1], "\n")
  tryCatch(
    {
      mod <- lmer(
        formula = as.formula(form_x),
        data = as.data.frame(data_x))
      
      return(mod)
    },
    error = function(cond) NA
  )
}


############################################
# fit models by looping through rows in nested data

# copy over nested data into a results nested data object
data_nested <- data_nested 

# initialize new column for mod results
data_nested$mods <- NA
data_nested$r2_cond <- NA
data_nested$r2_marg <- NA
data_nested$AIC <- NA
data_nested$comps <- NA
data_nested$iccs <- NA

# loop though and fit models, save results in list
for(i in 1:nrow(data_nested)){
  message(i)
  message(data_nested$siteID[[i]])
  message(data_nested$sizeCategory[[i]])

  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_nested$data[[i]],
                                form_x = data_nested$form[[i]])
  print(mod_result_i)
  
  if(class(mod_result_i)=="lmerMod"){

   data_nested$AIC[i] <- extractAIC(mod_result_i)[2] # 
    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_nested$mods[i] <- list(mod_result_i)
      
      # get r2 stats
      r2_i <- performance::r2(mod_result_i,tolerance = 1e-7)
      data_nested$r2_cond[i] <- r2_i[[1]]["Conditional R2"] #total var explained
      data_nested$r2_marg[i] <- r2_i[[1]]["Marginal R2"] #var explained by fixed effects
      
      # get ICCs
      icc_i <- performance::icc(mod_result_i,by_group = TRUE,tolerance = 1e-7)
      data_nested$comps[i] <- list(icc_i$Group)
      data_nested$iccs[i] <- list(icc_i$ICC)
      
    }
  }
}

# make a table that can write out as a csv
data_var_comps <- data_nested %>%
  select(-c(data,mods)) %>%
  arrange(siteID) %>% 
  unnest(c(comps,iccs))

data_var_comps$comps <- gsub("clipID:plotID:nlcdClass", "clipID", data_var_comps$comps)
data_var_comps$comps <- gsub("clipID:plotID", "clipID", data_var_comps$comps)
data_var_comps$comps <- gsub("plotID:nlcdClass", "plotID", data_var_comps$comps)
data_var_comps$comps <- gsub("nlcdClass", "NLCD", data_var_comps$comps)

data_var_comps$sizeCategoryDisplay <- gsub("0-1", "0-1 mm", data_var_comps$sizeCategory)
data_var_comps$sizeCategoryDisplay <- gsub("1-2", "1-2 mm", data_var_comps$sizeCategoryDisplay)
data_var_comps$sizeCategoryDisplay <- gsub("0-2", "0-2 mm", data_var_comps$sizeCategoryDisplay)
data_var_comps$sizeCategoryDisplay <- gsub("2-10", "2-10 mm", data_var_comps$sizeCategoryDisplay)

data_mods_nofixed_AIC <- data_nested %>%
  select(c(siteID,form, AIC)) 

#save(data_nested, file = "data_nested.rda")


################################
# Plot variance partitioning results (first available year only)
################
# data_var_comps has 87 rows, of which only 47 have iccs > 0
data_var_comps <- data_var_comps %>% left_join(stature_year, by = "siteID") %>% left_join(site_counts_mod, by = "siteID")

data_var_comps <- data_var_comps %>%  # currently each of the three sizeCategories and the total are plotted separately
  select(sizeCategoryDisplay, stature_yr, stature, years, plot_count, nlcd_count, nlcd_multi, siteID, y_var, comps, n_obs, iccs) %>%
  mutate(`Total Variance` = 'r2_cond') 

# Fig -- boxplot -- variance components by site
#graphics.off()

windows(18, 25, pointsize = 10)
data_var_comps %>% filter(stature == "large") %>% 
  ggplot(aes(siteID, iccs, color=comps, fill=comps)) +  
   geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + facet_wrap(vars(sizeCategoryDisplay, nlcd_multi), nrow = 4, scales="free_x") +
  geom_col(size = 0.2) +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_manual(values=c("deepskyblue", "#009E73","#E69F00","white")) + 
  xlab('Site ID') + 
  ylab(bquote(Total~variance~(R^2~conditional))) + 
  theme_light() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 4),axis.text.y=element_text(size = 6),
        strip.text.x = element_text(size = 4, margin = margin(0.1,0,0.1,0, "mm")),
        plot.title = element_text(size=8),
        legend.title = element_text(size=8), legend.key.size = unit(3, 'mm'), legend.text = element_text(size=6)) + 
  ggtitle("Large stature sites") 

savePlot("TotalVar_keep 1st bout - large stature.pdf", 'pdf')


windows(18, 25, pointsize = 10)
data_var_comps %>% filter(stature == "small") %>% 
  ggplot(aes(siteID, iccs, color=comps, fill=comps)) +  
   geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + facet_wrap(vars(sizeCategoryDisplay, nlcd_multi), nrow = 4, scales="free_x") +
  geom_col(size = 0.2) +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_manual(values=c("deepskyblue", "#009E73","#E69F00","white")) +
  xlab('Site ID') + 
  ylab(bquote(Total~variance~(R^2~conditional))) + 
  theme_light() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 6),axis.text.y=element_text(size = 6),
        strip.text.x = element_text(size = 4, margin = margin(0.1,0,0.1,0, "mm")),
        plot.title = element_text(size=8),
        legend.title = element_text(size=8), legend.key.size = unit(3, 'mm'), legend.text = element_text(size=6)) + 
  ggtitle("Small stature sites") 

savePlot("TotalVar_keep 1st bout - small stature.pdf", 'pdf')
 




###################################################################################################################
################# variance partitioning for multiyear sites #######################################################
###################################################################################################################

data_in_multiyear <- data_merge  %>% filter(years == 2) %>% filter(siteID != "GUAN" & siteID != "KONZ" & siteID != "MOAB" & siteID != "ORNL" & siteID != "SJER" & siteID != "STEI" & siteID != "TOOL" & siteID != "WREF")
 # manual removal of additional sites is because despite having erroneous years == 2 they really only had one year (perhaps years variable assigned before some filtering)

site_formulas_multiyear <- site_formulas
site_formulas_multiyear$form <- gsub('response_log ~ ','response_log ~ (1|year) + ',site_formulas_multiyear$form) # single year forumula, with year fixed effect tacked on
site_formulas_multiyear$form <- gsub('\\(1\\|clipID:plotID\\) \\+ ','',site_formulas_multiyear$form)
site_formulas_multiyear$form <- gsub('\\(1\\|clipID:plotID:nlcdClass\\) \\+ ','\\(1\\|plotID:nlcdClass\\) \\+ ',site_formulas_multiyear$form)
site_formulas_multiyear$form <- gsub(' \\+ \\(1\\|clipID\\)',' \\+ \\(1\\|plotID\\)',site_formulas_multiyear$form)
site_formulas_multiyear$form <- gsub('\\(1\\|plotID:nlcdClass\\) \\+ \\(1\\|plotID:nlcdClass\\) \\+ ','\\(1\\|plotID:nlcdClass\\) \\+ ',site_formulas_multiyear$form)

site_formulas_multiyear$form <- ifelse((site_formulas_multiyear$siteID == "UKFS" | site_formulas_multiyear$siteID == "TALL" | site_formulas_multiyear$siteID == "LENO"), 
                                       gsub('\\(1\\|year\\) \\+ ','',site_formulas_multiyear$form),site_formulas_multiyear$form)
site_formulas_multiyear$form <- ifelse((site_formulas_multiyear$siteID == "DEJU"), 
                                       gsub('\\(1\\|plotID:nlcdClass\\) \\+ \\(1\\|nlcdClass\\)','\\(1\\|plotID\\)',site_formulas_multiyear$form),site_formulas_multiyear$form)

site_formulas_multiyear$form <- ifelse((site_formulas_multiyear$siteID == "UNDE" | site_formulas_multiyear$siteID == "OAES") & site_formulas_multiyear$sizeCategory == "0-1", 
                                       gsub('\\(1\\|year\\) \\+ ','',site_formulas_multiyear$form),site_formulas_multiyear$form)
site_formulas_multiyear$form <- ifelse((site_formulas_multiyear$siteID == "ONAQ" & site_formulas_multiyear$sizeCategory == "2-10"), 
                                       gsub('\\+ \\(1\\|nlcdClass\\)','',site_formulas_multiyear$form),site_formulas_multiyear$form)
site_formulas_multiyear$form <- ifelse((site_formulas_multiyear$siteID == "JORN" | site_formulas_multiyear$siteID == "DCFS") & site_formulas_multiyear$sizeCategory == "2-10", 
                                       gsub('\\(1\\|year\\) \\+ ','',site_formulas_multiyear$form),site_formulas_multiyear$form)


data_nested_multiyear <- data_in_multiyear %>%
  group_by(sizeCategory, siteID) %>%
  nest() %>%
  left_join(site_formulas_multiyear, by = c("sizeCategory","siteID"))

# fit models by looping through rows in nested data

# copy over nested data into a results nested data object
data_nested_multiyear <- data_nested_multiyear # %>% filter(sizeCategory == "total")

# initialize new column for mod results
data_nested_multiyear$mods <- NA
data_nested_multiyear$r2_cond <- NA
data_nested_multiyear$r2_marg <- NA
data_nested_multiyear$AIC <- NA
data_nested_multiyear$comps <- NA
data_nested_multiyear$iccs <- NA

# loop though and fit models, save results in list
for(i in 1:nrow(data_nested_multiyear)){
  message(i)
  message(data_nested_multiyear$siteID[[i]])

  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_nested_multiyear$data[[i]],
                                form_x = data_nested_multiyear$form[[i]])
  print(mod_result_i)
  
  if(class(mod_result_i)=="lmerMod"){

   data_nested_multiyear$AIC[i] <- extractAIC(mod_result_i)[2] # 
    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_nested_multiyear$mods[i] <- list(mod_result_i)
      
      # get r2 stats
      r2_i <- performance::r2(mod_result_i,tolerance = 1e-7)
      data_nested_multiyear$r2_cond[i] <- r2_i[[1]]["Conditional R2"] #total var explained
      data_nested_multiyear$r2_marg[i] <- r2_i[[1]]["Marginal R2"] #var explained by fixed effects
      
      # get ICCs
      icc_i <- performance::icc(mod_result_i,by_group = TRUE,tolerance = 1e-7)
      data_nested_multiyear$comps[i] <- list(icc_i$Group)
      data_nested_multiyear$iccs[i] <- list(icc_i$ICC)
      
    }
  }
}

# make a table that can write out as a csv
data_var_comps_multiyear <- data_nested_multiyear %>%
  select(-c(data,mods)) %>%
  arrange(siteID) %>% 
  unnest(c(comps,iccs))

data_var_comps_multiyear$comps <- gsub("plotID:nlcdClass", "plotID", data_var_comps_multiyear$comps)
data_var_comps_multiyear$comps <- gsub("nlcdClass", "NLCD", data_var_comps_multiyear$comps)


data_var_comps_multiyear$sizeCategoryDisplay <- gsub("0-1", "0-1 mm", data_var_comps_multiyear$sizeCategory)
data_var_comps_multiyear$sizeCategoryDisplay <- gsub("1-2", "1-2 mm", data_var_comps_multiyear$sizeCategoryDisplay)
data_var_comps_multiyear$sizeCategoryDisplay <- gsub("0-2", "0-2 mm", data_var_comps_multiyear$sizeCategoryDisplay)
data_var_comps_multiyear$sizeCategoryDisplay <- gsub("2-10", "2-10 mm", data_var_comps_multiyear$sizeCategoryDisplay)


data_mods_nofixed_AIC_multiyear <- data_nested_multiyear %>%
  select(c(siteID,form, AIC)) 

################
# data_var_comps has 87 rows, of which only 47 have iccs > 0
nlcd_count <- site_counts %>% distinct(siteID, nlcd_multi)
data_var_comps_multiyear <- data_var_comps_multiyear %>% left_join(stature_year, by = "siteID") %>% left_join(nlcd_count, by = "siteID")

data_var_comps_multiyear <- data_var_comps_multiyear %>% # filter(!is.na(iccs) & iccs > 0) %>%
# d_plot <- data_var_comps %>% filter(!is.na(iccs) & iccs > 0 & n_per_plot >= 3.75) %>%  # only eventIDs with ~ 2 clips per plot x ~ 2 cores per clip (25 of 47)
# d_plot <- data_var_comps %>% filter(!is.na(iccs) & iccs > 0 & n_per_plot < 3.75 & clips_per_plot > 1) %>%  # only eventIDs with less than full design (above) but more than just core reps (below) (5 of 47)
# d_plot <- data_var_comps %>% filter(!is.na(iccs) & iccs > 0 & clips_per_plot == 1 ) %>%  # only eventIDs with 1 clip per plot (17 of 47)
  select(stature_yr, stature, years, nlcd_multi, siteID, y_var, comps, iccs, sizeCategoryDisplay) %>%
  mutate(`Total Variance` = 'r2_cond') 

#d_plot <- d_plot %>% filter(n_per_plot >= 3.75) # only eventIDs with ~ 2 cores per clip x 2 clips per plot
#d_plot <- d_plot %>% filter(n_per_plot >= 1.1 & n_per_plot < 3.5 & clips_per_plot > 1.1) # only eventIDs with ~ 2 cores per clip x 2 clips per plot

# Fig 1 -- boxplot -- variance components by site
#graphics.off()
windows(18, 20, pointsize = 10)

data_var_comps_multiyear %>% 
  ggplot(aes(siteID, iccs, color=comps, fill=comps)) +  
   geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + facet_wrap(vars(sizeCategoryDisplay, stature), nrow = 4, scales="free_x") +
  geom_col() +
  scale_color_manual(values=c("black", "black", "black","black","black")) +
  scale_fill_manual(values=c("#009E73","deepskyblue","black","white")) +
  ylab('Total variance (R^2 conditional)') +
  xlab('Site ID') + 
  theme_light() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1))
 
savePlot("TotalVar_keep both bouts.pdf", 'pdf')









#########################################################################################################
############## calculate true biomass mean of observed data #############################################


#  Calculate mean and variability of mass of cores within each plotID (and clipID)
data_merge$mass_g_m2 <- data_merge$response

bbc_per_m2 <- data_merge %>% filter(plotID != "SJER_053") %>% # was missing one of four cores, which made whole script crash during later resampling step
  group_by(sizeCategory, year, year2, domainID, siteID, plotID, clipID) %>% #, wst10cmDist, wst1cmDist, litterDepth) %>%
  summarise(dryMass = mean(mass_g_m2), sdDM = sd(mass_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass) %>% rename(bbc_g_m2 = dryMass)

bbc_true <- bbc_per_m2 %>% group_by(sizeCategory, year, year2, domainID, siteID) %>% summarise(bbc_g_m2_true = mean(bbc_g_m2, na.rm = TRUE)) %>%
  mutate("bbc-10%" = bbc_g_m2_true * 0.9, "bbc+10%" = bbc_g_m2_true * 1.1)
bbc_true_wide <- bbc_true %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(sizeCategory, domainID, siteID), names_from = year2, names_prefix = "mass_g_m2_true_", values_from = bbc_g_m2_true ) %>%
     mutate(mass_change_true = mass_g_m2_true_2 - mass_g_m2_true_1, mass_change_per_true = 100*mass_change_true/mass_g_m2_true_1,
            "mass_change-10%" = mass_change_true - (abs(mass_change_true) * 0.1), "mass_change+10%" = mass_change_true  + (abs(mass_change_true) * 0.1),
            "mass_change_per-5%" = mass_change_per_true - (abs(mass_change_per_true) * 0.05), "mass_change_per+5%" = mass_change_per_true  + (abs(mass_change_per_true) * 0.05))



############################################################################################################################################
################### Random resampling (very time consuming); later summarized for comparison to observed mean ##############################
repetitions <- 1000


#bbc_fullb <- data_merge  %>% filter(plotID != "SJER_053") # was missing one of four cores, which made whole script crash


## sample two cores in each of 2 clipIDs (both clips always get two cores - this resamples the current full sample design)
# resampling summary: 1) create list of unique clipIDs, 2) slice 2 samples from within each clipID. 
# Enforces that BOTH clipIDs within plot are selected, but could draw the same north or south core twice in addition to randomly reselecting 1 north and 1 south as in the original data 
# resample of full design for large stature sites
sample_set <- 2
calc_2core2clip <- data.frame()
reps <- data.frame()
groups <- unique(data_merge$clipID)
sizeCategories <- unique(data_merge$sizeCategory)


# Create a directory to store intermediate results
dir.create("calc_2core2clip_chunks", showWarnings = FALSE)

for (j in seq_along(sizeCategories)) {
  for (i in seq_along(groups)) {
    tryCatch({
      reps <- data_merge %>%
        filter(clipID == groups[i], sizeCategory == sizeCategories[j]) %>%
        infer::rep_slice_sample(n = sample_set, replace = TRUE, reps = repetitions)  %>% dplyr::ungroup()
      reps$sample <- rep(1:sample_set, times = repetitions)
      reps$set <- "calc_2core2clip"
      
      # Save each chunk to disk
      filename <- paste0("calc_2core2clip_chunks/reps_", sizeCategories[j], "_", groups[i], ".rds")
      saveRDS(reps, filename)
      
    }, error = function(e) {
      message(paste("Error with group", groups[i], "and sizeCategory", sizeCategories[j], ":", e$message))
    })
  }
}

# Combine all saved chunks
calc_2core2clip_list <- list.files("calc_2core2clip_chunks", full.names = TRUE)
calc_2core2clip <- do.call(bind_rows, lapply(calc_2core2clip_list, readRDS))
save(calc_2core2clip, file = "calc_2core2clip.rda")

calc_2core1clip <- calc_2core2clip %>% filter(stature == "small") %>% mutate(set = "calc_2core1clip")
#save(calc_2core1clip, file = "calc_2core1clip.rda")

#####################################
Sys.time()
## sample one randomly selected core in each of 2 clipIDs ()
sample_set <- 1
calc_1core2clip <- data.frame()
reps <- data.frame()
groups <- unique(data_merge$clipID)
sizeCategories <- unique(data_merge$sizeCategory)

# Create a directory to store intermediate results
dir.create("calc_1core2clip_chunks", showWarnings = FALSE)

for (j in seq_along(sizeCategories)) {
  for (i in seq_along(groups)) {
    tryCatch({
      reps <- data_merge %>%
        filter(clipID == groups[i], sizeCategory == sizeCategories[j]) %>%
        infer::rep_slice_sample(n = sample_set, replace = TRUE, reps = repetitions) %>% dplyr::ungroup()
      reps$sample <- rep(1:sample_set, times = repetitions)
      reps$set <- "calc_1core2clip"
      
      # Save each chunk to disk
      filename <- paste0("calc_1core2clip_chunks/reps_", sizeCategories[j], "_", groups[i], ".rds")
      saveRDS(reps, filename)
      
    }, error = function(e) {
      message(paste("Error with group", groups[i], "and sizeCategory", sizeCategories[j], ":", e$message))
    })
  }
}

# Combine all saved chunks
calc_1core2clip_list <- list.files("calc_1core2clip_chunks", full.names = TRUE)
calc_1core2clip <- do.call(bind_rows, lapply(calc_1core2clip_list, readRDS))
#save(calc_1core2clip, file = "calc_1core2clip.rda")

#################################
Sys.time() 
## sample one core (any of the 4 clip by core combos could be randomly picked)
# resampling summary: 1) create list of unique plotIDs, 2) slice just 1 sample from within each plotID
# Enforces that just one clipID and one core per plot are selected. 
# lowest intensity sampling scenario for both large and small statured sites

sample_set <- 1
calc_1core1clip <- data.frame()
reps <- data.frame()
groups <- unique(data_merge$plotID)
sizeCategories <- unique(data_merge$sizeCategory)


# Create a directory to store intermediate results
dir.create("calc_1core1clip_chunks", showWarnings = FALSE)

for (j in seq_along(sizeCategories)) {
  for (i in seq_along(groups)) {
    tryCatch({
      reps <- data_merge %>%
        filter(plotID == groups[i] & sizeCategory == sizeCategories[j]) %>%
        infer::rep_slice_sample(n = sample_set, replace = TRUE, reps = repetitions)    %>% dplyr::ungroup()
      reps$sample<-rep(c(1:sample_set),times=repetitions)
      reps$set <- "calc_1core1clip"
      
      # Save each chunk to disk
      filename <- paste0("calc_1core1clip_chunks/reps_", sizeCategories[j], "_", groups[i], ".rds")
      saveRDS(reps, filename)
      
    }, error = function(e) {
      message(paste("Error with group", groups[i], "and sizeCategory", sizeCategories[j], ":", e$message))
    })
  }
}

# Combine all saved chunks
calc_1core1clip_list <- list.files("calc_1core1clip_chunks", full.names = TRUE)
calc_1core1clip <- do.call(bind_rows, lapply(calc_1core1clip_list, readRDS))
#save(calc_1core1clip, file = "calc_1core1clip.rda")


# load(file='calc_2core2clip.rda')
# load(file='calc_2core1clip.rda')
# load(file='calc_1core2clip.rda')
# load(file='calc_1core1clip.rda')


############################################################################################################
############# Summarize change in mass with different sampling scenarios (for comparison with mean) ########
###################################
##########  2core2clip ###########


rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub)
#   Calculate mean and variability of mass (within a plot)
calc_2core2clip$mass_g_m2 <- calc_2core2clip$response
bbc_2core2clip<- calc_2core2clip%>%
  group_by(sizeCategory, year, year2, domainID, siteID, plotID, replicate) %>%  # sampleID, wst10cmDist, wst1cmDist, litterDepth, 
  summarise(dryMass = mean(mass_g_m2), sdDM = sd(mass_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass) %>% rename(mass_g_m2 = dryMass)
#save(bbc_2core2clip, file = "bbc_2core2clip.rda")
#load(file='bbc_2core2clip.rda')

bbc_per_m2_multiyear_sub <- bbc_2core2clip%>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(sizeCategory,domainID,siteID,replicate) %>% do(tidy(glm(mass_g_m2 ~ year2, data = .) )) # %>% filter(term == "year22")
################# no p-value, not sure why full model results are not visible
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_2core2clip_summary <- bbc_year_sub %>% group_by(sizeCategory,domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_2core2clip = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_2core2clip = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_2core2clip = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )
#save(bbc_year_2core2clip_summary, file = "bbc_year_2core2clip_summary.rda")
#load(file='bbc_year_2core2clip_summary.rda')


check_2core2clip<- bbc_2core2clip%>% group_by(sizeCategory, year, year2, domainID, siteID, replicate) %>% summarise(mass_g_m2 = mean(mass_g_m2, na.rm = TRUE)) # calculate mean mass per site (mean per plot is input)

  check_2core2clip_wide <- check_2core2clip %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(sizeCategory, domainID, siteID, replicate), names_from = year2, names_prefix = "mass_g_m2_", values_from = mass_g_m2 ) %>%
     mutate(mass_change = mass_g_m2_2 - mass_g_m2_1, mass_change_per = 100*mass_change/mass_g_m2_1)
    check_2core2clip_wide <- merge(check_2core2clip_wide, bbc_true_wide, by = c("sizeCategory","domainID", "siteID"))
    check_2core2clip_wide$bbcFlag <- ifelse(check_2core2clip_wide$mass_change > check_2core2clip_wide$`mass_change-10%` & check_2core2clip_wide$mass_change < check_2core2clip_wide$`mass_change+10%`, 1, 0)
    check_2core2clip_wide$bbcFlag_per <- ifelse(check_2core2clip_wide$mass_change_per > check_2core2clip_wide$`mass_change_per-5%` & check_2core2clip_wide$mass_change_per < check_2core2clip_wide$`mass_change_per+5%`, 1, 0)
  check_2core2clip_wide_summary <- check_2core2clip_wide %>% group_by(sizeCategory, domainID, siteID) %>%
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_2core2clip = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      bbc_reps_change_inRange = sum(bbcFlag_per), reps_Number = n(), bbc_per_change_inRange_2core2clip = round((bbc_reps_change_inRange/reps_Number)*100, digits=1)) %>%
      select(c(-bbc_reps_inRange, -reps_Number, -bbc_reps_change_inRange) )
 
check_2core2clip <- merge(check_2core2clip, bbc_true, by = c("sizeCategory","year", "year2", "domainID", "siteID"))
check_2core2clip$bbcFlag <- ifelse(check_2core2clip$mass_g_m2 > check_2core2clip$`bbc-10%` & check_2core2clip$mass_g_m2 < check_2core2clip$`bbc+10%`, 1, 0)

check_2core2clip_summary <- check_2core2clip %>% group_by(sizeCategory, year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_2core2clip = round((bbc_reps_inRange/reps_Number)*100, digits=1)) %>%
  select(c(-bbc_reps_inRange, -reps_Number) )



##########  2core1clip ###########

rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub)
#   Calculate mean and variability of mass
calc_2core1clip$mass_g_m2 <- calc_2core1clip$response
bbc_2core1clip<- calc_2core1clip%>%
  group_by(sizeCategory, year, year2, domainID, siteID, plotID, replicate) %>%  # sampleID, wst10cmDist, wst1cmDist, litterDepth,
  summarise(dryMass = mean(mass_g_m2), sdDM = sd(mass_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass) %>% rename(mass_g_m2 = dryMass)
#save(bbc_2core1clip, file = "bbc_2core1clip.rda")
#load(file='bbc_2core1clip.rda')

bbc_per_m2_multiyear_sub <- bbc_2core1clip%>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(sizeCategory,domainID,siteID,replicate) %>% do(tidy(glm(mass_g_m2 ~ year2, data = .) )) # %>% filter(term == "year22")
################# no p-value, not sure why full model results are not visible
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_2core1clip_summary <- bbc_year_sub %>% group_by(sizeCategory,domainID, siteID) %>%
  summarise(reps_Number = n(), estimate_2core2clip = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_2core2clip = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_2core2clip = round((year_reps_sig/reps_Number)*100, digits=1)) %>%
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )
#save(bbc_year_2core1clip_summary, file = "bbc_year_2core1clip_summary.rda")
#load(file='bbc_year_2core1clip_summary.rda')


check_2core1clip<- bbc_2core1clip%>% group_by(sizeCategory, year, year2, domainID, siteID, replicate) %>% summarise(mass_g_m2 = mean(mass_g_m2, na.rm = TRUE))

  check_2core1clip_wide <- check_2core1clip %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(sizeCategory, domainID, siteID, replicate), names_from = year2, names_prefix = "mass_g_m2_", values_from = mass_g_m2 ) %>%
     mutate(mass_change = mass_g_m2_2 - mass_g_m2_1, mass_change_per = 100*mass_change/mass_g_m2_1)
    check_2core1clip_wide <- merge(check_2core1clip_wide, bbc_true_wide, by = c("sizeCategory","domainID", "siteID"))
    check_2core1clip_wide$bbcFlag <- ifelse(check_2core1clip_wide$mass_change > check_2core1clip_wide$`mass_change-10%` & check_2core1clip_wide$mass_change < check_2core1clip_wide$`mass_change+10%`, 1, 0)
    check_2core1clip_wide$bbcFlag_per <- ifelse(check_2core1clip_wide$mass_change_per > check_2core1clip_wide$`mass_change_per-5%` & check_2core1clip_wide$mass_change_per < check_2core1clip_wide$`mass_change_per+5%`, 1, 0)
  check_2core1clip_wide_summary <- check_2core1clip_wide %>% group_by(sizeCategory, domainID, siteID) %>%
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_2core1clip = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      bbc_reps_change_inRange = sum(bbcFlag_per), reps_Number = n(), bbc_per_change_inRange_2core1clip = round((bbc_reps_change_inRange/reps_Number)*100, digits=1)) %>%
      select(c(-bbc_reps_inRange, -reps_Number, -bbc_reps_change_inRange) )

check_2core1clip <- merge(check_2core1clip, bbc_true, by = c("sizeCategory","year", "year2", "domainID", "siteID"))
check_2core1clip$bbcFlag <- ifelse(check_2core1clip$mass_g_m2 > check_2core1clip$`bbc-10%` & check_2core1clip$mass_g_m2 < check_2core1clip$`bbc+10%`, 1, 0)

check_2core1clip_summary <- check_2core1clip %>% group_by(sizeCategory, year, year2, domainID, siteID) %>%
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_2core1clip = round((bbc_reps_inRange/reps_Number)*100, digits=1)) %>%
  select(c(-bbc_reps_inRange, -reps_Number) )




###################################
##########  1core2clip  ###########


rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub)
#   Calculate mean and variability of mass
calc_1core2clip$mass_g_m2 <- calc_1core2clip$response
bbc_1core2clip <- calc_1core2clip %>%
  group_by(sizeCategory, year, year2, domainID, siteID, plotID, replicate) %>%  # sampleID, wst10cmDist, wst1cmDist, litterDepth, 
  summarise(dryMass = mean(mass_g_m2), sdDM = sd(mass_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass) %>% rename(mass_g_m2 = dryMass)
#save(bbc_1core2clip, file = "bbc_1core2clip.rda")
#load(file='bbc_1core2clip.rda')

bbc_per_m2_multiyear_sub <- bbc_1core2clip %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(sizeCategory, domainID,siteID,replicate) %>% do(tidy(glm(mass_g_m2 ~ year2, data = .) )) # %>% filter(term == "year22")
  ################# no p-value, not sure why full model results are not visible
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_1core2clip_summary <- bbc_year_sub %>% group_by(sizeCategory, domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_1core2clip = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_1core2clip = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_1core2clip = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )
save(bbc_year_1core2clip_summary, file = "bbc_year_1core2clip_summary.rda")
#load(file='bbc_year_1core2clip_summary.rda')


check_1core2clip <- bbc_1core2clip %>% group_by(sizeCategory, year, year2, domainID, siteID, replicate) %>% summarise(mass_g_m2 = mean(mass_g_m2, na.rm = TRUE)) 

  check_1core2clip_wide <- check_1core2clip  %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(sizeCategory, domainID, siteID, replicate), names_from = year2, names_prefix = "mass_g_m2_", values_from = mass_g_m2 ) %>%
     mutate(mass_change = mass_g_m2_2 - mass_g_m2_1, mass_change_per = 100*mass_change/mass_g_m2_1)
    check_1core2clip_wide <- merge(check_1core2clip_wide, bbc_true_wide, by = c("sizeCategory","domainID", "siteID"))
    check_1core2clip_wide$bbcFlag <- ifelse(check_1core2clip_wide$mass_change > check_1core2clip_wide$`mass_change-10%` & check_1core2clip_wide$mass_change < check_1core2clip_wide$`mass_change+10%`, 1, 0)
    check_1core2clip_wide$bbcFlag_per <- ifelse(check_1core2clip_wide$mass_change_per > check_1core2clip_wide$`mass_change_per-5%` & check_1core2clip_wide$mass_change_per < check_1core2clip_wide$`mass_change_per+5%`, 1, 0)
  check_1core2clip_wide_summary <- check_1core2clip_wide %>% group_by(sizeCategory, domainID, siteID) %>%
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_1core2clip = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      bbc_reps_change_inRange = sum(bbcFlag_per), reps_Number = n(), bbc_per_change_inRange_1core2clip = round((bbc_reps_change_inRange/reps_Number)*100, digits=1)) %>%
      select(c(-bbc_reps_inRange, -reps_Number, -bbc_reps_change_inRange) )
 
check_1core2clip <- merge(check_1core2clip, bbc_true, by = c("sizeCategory","year", "year2", "domainID", "siteID"))
check_1core2clip$bbcFlag <- ifelse(check_1core2clip$mass_g_m2 > check_1core2clip$`bbc-10%` & check_1core2clip$mass_g_m2 < check_1core2clip$`bbc+10%`, 1, 0)
#check_1core2clip$bbcFlag <- ifelse(check_1core2clip$mass_g_m2 > check_1core2clip$`bbc-20%` & check_1core2clip$mass_g_m2 < check_1core2clip$`bbc+20%`, 1, 0)

check_1core2clip_summary <- check_1core2clip %>% group_by(sizeCategory, year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_1core2clip = round((bbc_reps_inRange/reps_Number)*100, digits=1)) %>%
  select(c(-bbc_reps_inRange, -reps_Number) )


###################################
##########  1core1clip  ###########


rm(bbc_per_m2_multiyear_sub); rm(bbc_year_sub)
#   Calculate mean and variability of mass
calc_1core1clip$mass_g_m2 <- calc_1core1clip$response
bbc_1core1clip <- calc_1core1clip %>%
  group_by(sizeCategory, year, year2, domainID, siteID, plotID, replicate) %>%  # sampleID, wst10cmDist, wst1cmDist, litterDepth, 
  summarise(dryMass = mean(mass_g_m2), sdDM = sd(mass_g_m2), nDM = n(), seDM = sdDM/sqrt(nDM), relseDM = 100 * seDM / dryMass) %>% rename(mass_g_m2 = dryMass)
#save(bbc_1core1clip, file = "bbc_1core1clip.rda")
#load(file='bbc_1core1clip.rda')

bbc_1core1clip_moist <- bbc_1core1clip # %>% left_join(sls_soilMoisture, by = c("year","plotID"))
#bbc_1core1clip_moist <- bbc_1core1clip_moist %>% left_join(sls_soilMoisture_ave, by = c("year","siteID"))
#bbc_1core1clip_moist$soilMoisture <- ifelse(is.na(bbc_1core1clip_moist$soilMoisture), bbc_1core1clip_moist$soilMoisture_ave, bbc_1core1clip_moist$soilMoisture)

bbc_per_m2_multiyear_sub <- bbc_1core1clip %>% filter(year2 ==1 | year2 ==2)

bbc_year_sub <- bbc_per_m2_multiyear_sub %>% group_by(sizeCategory, domainID,siteID,replicate) %>% do(tidy(glm(mass_g_m2 ~ year2, data = .) )) # %>% filter(term == "year22")
  ################# no p-value, not sure why full model results are not visible
bbc_year_sub$yearFlag <- ifelse(bbc_year_sub$p.value < 0.05, 1, 0)
bbc_year_sub$yearDir <- ifelse(bbc_year_sub$estimate > 0, 1, 0)
bbc_year_1core1clip_summary <- bbc_year_sub %>% group_by(sizeCategory, domainID, siteID) %>% 
  summarise(reps_Number = n(), estimate_1core1clip = mean(estimate), year_reps_pos = sum(yearDir), year_per_pos_1core1clip = round((year_reps_pos/reps_Number)*100, digits=1),
            year_reps_sig = sum(yearFlag), year_per_sig_1core1clip = round((year_reps_sig/reps_Number)*100, digits=1)) %>% 
   select(c(-reps_Number,-year_reps_pos,-year_reps_sig) )
#save(bbc_year_1core1clip_summary, file = "bbc_year_1core1clip_summary.rda")
#load(file='bbc_year_1core1clip_summary.rda')


check_1core1clip <- bbc_1core1clip %>% group_by(sizeCategory, year, year2, domainID, siteID, replicate) %>% summarise(mass_g_m2 = mean(mass_g_m2, na.rm = TRUE)) 

  check_1core1clip_wide <- check_1core1clip  %>% filter(year2 == 1 | year2 ==2) %>% pivot_wider(id_cols = c(sizeCategory, domainID, siteID, replicate), names_from = year2, names_prefix = "mass_g_m2_", values_from = mass_g_m2 ) %>%
     mutate(mass_change = mass_g_m2_2 - mass_g_m2_1, mass_change_per = 100*mass_change/mass_g_m2_1)
    check_1core1clip_wide <- merge(check_1core1clip_wide, bbc_true_wide, by = c("sizeCategory","domainID", "siteID"))
    check_1core1clip_wide$bbcFlag <- ifelse(check_1core1clip_wide$mass_change > check_1core1clip_wide$`mass_change-10%` & check_1core1clip_wide$mass_change < check_1core1clip_wide$`mass_change+10%`, 1, 0)
    check_1core1clip_wide$bbcFlag_per <- ifelse(check_1core1clip_wide$mass_change_per > check_1core1clip_wide$`mass_change_per-5%` & check_1core1clip_wide$mass_change_per < check_1core1clip_wide$`mass_change_per+5%`, 1, 0)
  check_1core1clip_wide_summary <- check_1core1clip_wide %>% group_by(sizeCategory, domainID, siteID) %>%
      summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_1core1clip = round((bbc_reps_inRange/reps_Number)*100, digits=1),
      bbc_reps_change_inRange = sum(bbcFlag_per), reps_Number = n(), bbc_per_change_inRange_1core1clip = round((bbc_reps_change_inRange/reps_Number)*100, digits=1)) %>%
      select(c(-bbc_reps_inRange, -reps_Number, -bbc_reps_change_inRange) )
 
check_1core1clip <- merge(check_1core1clip, bbc_true, by = c("sizeCategory", "year", "year2", "domainID", "siteID"))
check_1core1clip$bbcFlag <- ifelse(check_1core1clip$mass_g_m2 > check_1core1clip$`bbc-10%` & check_1core1clip$mass_g_m2 < check_1core1clip$`bbc+10%`, 1, 0)

check_1core1clip_summary <- check_1core1clip %>% group_by(sizeCategory, year, year2, domainID, siteID) %>% 
  summarise(bbc_reps_inRange = sum(bbcFlag), reps_Number = n(), bbc_per_inRange_1core1clip = round((bbc_reps_inRange/reps_Number)*100, digits=1)) %>%
  select(c(-bbc_reps_inRange, -reps_Number) )


############## Combine the summaries of mean biomass for the different sampling scenarios #####################################

bbc_site_list <- list(check_2core2clip_summary, check_2core1clip_summary, check_1core2clip_summary, check_1core1clip_summary)  # check_4cores_summary, check_2core1clip_summary,     
bbc_site_summary <- bbc_site_list %>% reduce(full_join, by=c("year", "year2", "domainID", "siteID","sizeCategory") )
 View(bbc_site_summary)
 write.table(bbc_site_summary, paste0("year_site_bbc_",Sys.Date(),".txt"), sep = "\t", row.names=F)

bbc_site_summary_long <- bbc_site_summary %>% 
  pivot_longer(
    cols = !c(sizeCategory,year, year2, domainID, siteID), 
    names_to = "subsamples", 
    values_to = "percent"
  )


bbc_site_summary_long$subsamples <- gsub("bbc_per_inRange_", "", bbc_site_summary_long$subsamples)
 bbc_site_summary_long$subsamples <- gsub("2core1clip", 2.1, bbc_site_summary_long$subsamples)
bbc_site_summary_long$subsamples <- gsub("2core2clip", 4, bbc_site_summary_long$subsamples)
bbc_site_summary_long$subsamples <- gsub("1core2clip", 2, bbc_site_summary_long$subsamples)
bbc_site_summary_long$subsamples <- gsub("1core1clip", 1, bbc_site_summary_long$subsamples)
bbc_site_summary_long$subsamples <- as.numeric(bbc_site_summary_long$subsamples)
bbc_site_summary_long$year2 <- as.character(bbc_site_summary_long$year2)
bbc_site_summary_long <- bbc_site_summary_long %>% relocate(subsamples, .before=year2)
bbc_site_summary_long <- bbc_site_summary_long %>% left_join(stature_class, by = "siteID")
bbc_site_summary_long <- bbc_site_summary_long %>% filter(!(stature == "small" & subsamples == 4))


windows(9,6)
bbc_site_summary_large <- bbc_site_summary_long %>% filter(stature == "large") %>% filter(subsamples != 2.1) 
ggplot(data = bbc_site_summary_large, aes(x=subsamples, y=percent)) +  
  scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
  geom_point(aes(colour=sizeCategory)) + 
   geom_line(aes(colour=sizeCategory, linetype = year2)) + 
#  scale_linetype_manual(values=c("twodash","solid")) + 
  geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) +  # ylim(75, 100) + 
  scale_x_reverse() + 
  facet_wrap(~ siteID , ncol = 8) + 
  labs(x = "Sample design (4 = resampled current (2 clip x 2 core),  2 = 2 clip x 1 core,  1 = 1 clip x 1 core)") + 
  labs(y = "Percent of resample iterations within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
    theme_minimal() + ggtitle("Large stature sites") # + theme(legend.position="none") 

 savePlot('site_means_largeStature','pdf')  # currently Figure 3


windows(9,6)
bbc_site_summary_small <- bbc_site_summary_long %>% filter(stature == "small") %>% filter(subsamples != 2)
bbc_site_summary_small$subsamples <- as.numeric(gsub(2.1, 2, bbc_site_summary_small$subsamples))
ggplot(data = bbc_site_summary_small, aes(x=subsamples, y=percent)) + 
  scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
  geom_point(aes(colour=sizeCategory)) + 
   geom_line(aes(colour=sizeCategory, linetype = year2)) + 
#  scale_linetype_manual(values=c("twodash","solid")) + 
  geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) +  # ylim(75, 100) + 
  scale_x_reverse(breaks=c(2, 1)) + 
  facet_wrap(~ siteID , ncol = 8) + 
  labs(x = "Sample design (2 = current (1 clip x 2 core),  1 = 1 clip x 1 core)") + 
  labs(y = "Percent of resample iterations within +/- 10% of observed mean") + theme(axis.title = element_text(size = 40)) + theme(axis.text.y = element_text(size = 20))  + 
    theme_minimal()  + ggtitle("Small stature sites") # + theme(legend.position="right")

 savePlot('site_means_smallStature','pdf')  # currently Figure 4

# write.table(bbc_year_summary_subset, paste0("year_summary_bbc_subset_",Sys.Date(),".txt"), sep = "\t", row.names=F)



