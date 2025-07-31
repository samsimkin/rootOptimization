library(tidyverse)
library(lme4)

# The prerequisite for running this script is running the "BBC_field_02 variance partitioning and bootstrapping.R" script.


# full data, look at year, sizeCategory and interaction in multiyear data

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

# test fn_fit_models
fn_fit_models_fixed <- function(data_x, form_x) {
  cat(as.data.frame(data_x)$siteID[1], "\n")
  tryCatch(
    {
      mod <- glm(
        formula = as.formula(form_x),
        data = as.data.frame(data_x))
      
      return(mod)
    },
    error = function(cond) NA
  )
}

rm(bbc_fixed)
rm(data_mods_fixed_nested)

site_formulas <- read.delim("site_formulas.txt", header = TRUE, sep = "\t")


bbc_fixed <- data_merge %>% mutate(response_log = log(response)) %>% select(-years) # years var was calculated before dropping sites without full design
     year_class <- bbc_fixed %>% distinct(eventID, siteID) %>% group_by(siteID) %>% summarise(years = n() )
bbc_fixed <- bbc_fixed %>% left_join(year_class, by = "siteID") %>% filter(years == 2)


site_formulas_w_fixed_effects <- site_formulas
site_formulas_w_fixed_effects_noSize <- site_formulas_w_fixed_effects %>% select(-sizeCategory) %>% distinct()
site_formulas_w_fixed_effects_noSize$form <- 'response_log ~ (1|nlcdClass)' # toggle on or off
site_formulas_w_fixed_effects_noSize$form <- ifelse(site_formulas_w_fixed_effects_noSize$siteID == "MLBS" | site_formulas_w_fixed_effects_noSize$siteID == "DSNY" | 
                                                      site_formulas_w_fixed_effects_noSize$siteID == "WOOD", 'response_log ~ (1|plotID)', site_formulas_w_fixed_effects_noSize$form)
site_formulas_w_fixed_effects_noSize$form_fixed <- 'response_log ~ year' # toggle on or off

data_mods_fixed_nested <- bbc_fixed %>%
  group_by(siteID, sizeCategory) %>%
  nest() %>%
  left_join(site_formulas_w_fixed_effects_noSize, by = "siteID") 

# initialize new columns for mod results
data_mods_fixed_nested$mods <- NA
data_mods_fixed_nested$r2_cond <- NA
data_mods_fixed_nested$AIC <- NA
data_mods_fixed_nested$mods_fixed <- NA
data_mods_fixed_nested$r2_fixed_cond <- NA
data_mods_fixed_nested$r2_fixed_marg <- NA
data_mods_fixed_nested$AIC_w_fixed <- NA
data_mods_fixed_nested$year_coeff_p <- NA

data_mods_fixed_nested$modComp_p <- NA

data_mods_fixed_nested <- data_mods_fixed_nested %>% filter(!(siteID == "BART" & sizeCategory == "2-10"))
# loop though and fit models, save results in list
for(i in 1:nrow(data_mods_fixed_nested)){
  message(i)
  message(data_mods_fixed_nested$siteID[[i]])
  tryCatch({
  mod_fixed_result_i <- NA
  mod_fixed_result_i <- fn_fit_models_fixed(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form_fixed[[i]])
   coeff_fixed_result_i<- as.data.frame(coef(summary(mod_fixed_result_i)))
    year_coeff_i <- coeff_fixed_result_i %>% filter(grepl("year", row.names(coeff_fixed_result_i)))
 
     data_mods_fixed_nested$year_coeff_p[i] <- year_coeff_i$'Pr(>|t|'
  print(mod_fixed_result_i)
      data_mods_fixed_nested$AIC_w_fixed[i] <- extractAIC(mod_fixed_result_i)[2] # 
    
    if(!performance::check_singularity(mod_fixed_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods_fixed[i] <- list(mod_fixed_result_i)
      
      # get r2 stats
      r2_fixed_i <- performance::r2(mod_fixed_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_fixed_cond[i] <- r2_fixed_i[[1]]["Conditional R2"] #total var explained
      data_mods_fixed_nested$r2_fixed_marg[i] <- r2_fixed_i[[1]]["Marginal R2"] #var explained by fixed effects
    }

  
  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form[[i]])
  print(mod_result_i)
  
      data_mods_fixed_nested$AIC[i] <- extractAIC(mod_result_i)[2] # 
    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods[i] <- list(mod_result_i)
      
      # get r2 stats
      r2_i <- performance::r2(mod_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_cond[i] <- r2_i[[1]]["Conditional R2"] #total var explained
 
    }
  modComp_ANOVA_i <- anova(mod_fixed_result_i, mod_result_i)
 #   data_mods_fixed_nested$modComp_p[i] <-  modComp_ANOVA_i$"Pr(>Chisq)"[2]

  }, error = function(e) {
    message(paste("Error in row", i, ":", e$message))
  })
}

data_mods_fixed_nested$AIC_improve <- data_mods_fixed_nested$AIC - data_mods_fixed_nested$AIC_w_fixed # compare model with moisture variable to model with just random nlcdClass
mod_comp_full <- data_mods_fixed_nested



###########################
# data subset 2core2clip, look at year, sizeCategory and interaction in multiyear data for large stature sites

rm(bbc_fixed)
rm(data_mods_fixed_nested)

bbc_fixed <- calc_2core2clip %>% mutate(response_log = log(response))  %>% filter(replicate <= 50) %>% filter(siteID == "BART" | siteID == "DEJU" | siteID == "SCBI" | siteID == "TALL") %>%
  select(-years) # years var was calculated before dropping sites without full design
     year_class <- bbc_fixed %>% distinct(eventID, siteID) %>% group_by(siteID) %>% summarise(years = n() )
# site filtering is to large sites with a year effect in raw data
bbc_fixed <- bbc_fixed %>% left_join(year_class, by = "siteID") %>% filter(years == 2) # %>% filter(nlcd_count > 1) 


data_mods_fixed_nested <- bbc_fixed %>%
  group_by(siteID, sizeCategory,replicate) %>%
  nest() %>%
  left_join(site_formulas_w_fixed_effects_noSize, by = "siteID") 

# initialize new columns for mod results
data_mods_fixed_nested$mods <- NA
data_mods_fixed_nested$r2_cond <- NA
data_mods_fixed_nested$AIC <- NA
data_mods_fixed_nested$mods_fixed <- NA
data_mods_fixed_nested$r2_fixed_cond <- NA
data_mods_fixed_nested$r2_fixed_marg <- NA
data_mods_fixed_nested$AIC_w_fixed <- NA
data_mods_fixed_nested$year_coeff_p <- NA

data_mods_fixed_nested$modComp_p <- NA

# loop though and fit models, save results in list
for(i in 1:nrow(data_mods_fixed_nested)){
  message(i)
  message(data_mods_fixed_nested$siteID[[i]])

  tryCatch({
    mod_fixed_result_i <- NA
  mod_fixed_result_i <- fn_fit_models_fixed(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form_fixed[[i]])
   coeff_fixed_result_i<- as.data.frame(coef(summary(mod_fixed_result_i)))
    year_coeff_i <- coeff_fixed_result_i %>% filter(grepl("year", row.names(coeff_fixed_result_i))) 
     data_mods_fixed_nested$year_coeff_p[i] <- year_coeff_i$'Pr(>|t|'
  print(mod_fixed_result_i)
      data_mods_fixed_nested$AIC_w_fixed[i] <- extractAIC(mod_fixed_result_i)[2] # 
    
    if(!performance::check_singularity(mod_fixed_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods_fixed[i] <- list(mod_fixed_result_i)
      
      # get r2 stats
      r2_fixed_i <- performance::r2(mod_fixed_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_fixed_cond[i] <- r2_fixed_i[[1]]["Conditional R2"] #total var explained
      data_mods_fixed_nested$r2_fixed_marg[i] <- r2_fixed_i[[1]]["Marginal R2"] #var explained by fixed effects
    }

  
  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form[[i]])
  print(mod_result_i)
  
      data_mods_fixed_nested$AIC[i] <- extractAIC(mod_result_i)[2] # 
    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods[i] <- list(mod_result_i)
      
      # get r2 stats
      r2_i <- performance::r2(mod_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_cond[i] <- r2_i[[1]]["Conditional R2"] #total var explained
 
    }
  modComp_ANOVA_i <- anova(mod_fixed_result_i, mod_result_i)

  }, error = function(e) {
    message(paste("Error in row", i, ":", e$message))
  })
}

data_mods_fixed_nested$AIC_improve <- data_mods_fixed_nested$AIC - data_mods_fixed_nested$AIC_w_fixed # compare model with moisture variable to model with just random nlcdClass
mod_comp_2core2clip <- data_mods_fixed_nested


mod_comp_2core2clip$yr_eff_lt_05_2core2clip <- if_else(mod_comp_2core2clip$year_coeff_p <0.05, 1, 0, 0)
mod_comp_per_2core2clip <- mod_comp_2core2clip %>% group_by(siteID, sizeCategory) %>% summarise(
                                                                          yr_eff_lt_05_per_2core2clip = 100 * sum(yr_eff_lt_05_2core2clip) / n())


###########################
# data subset 1core2clip, look at year, sizeCategory and interaction in multiyear data for large stature sites

rm(bbc_fixed)
rm(data_mods_fixed_nested)

bbc_fixed <- calc_1core2clip %>% mutate(response_log = log(response))  %>% filter(replicate <= 50) %>% filter(siteID == "BART" | siteID == "DEJU" | siteID == "SCBI" | siteID == "TALL") %>%
  select(-years) # years var was calculated before dropping sites without full design
# site filtering is to large sites with a year effect in raw data
     year_class <- bbc_fixed %>% distinct(eventID, siteID) %>% group_by(siteID) %>% summarise(years = n() )
bbc_fixed <- bbc_fixed %>% left_join(year_class, by = "siteID") %>% filter(years == 2) # %>% filter(nlcd_count > 1) 


data_mods_fixed_nested <- bbc_fixed %>%
  group_by(siteID, sizeCategory,replicate) %>%
  nest() %>%
  left_join(site_formulas_w_fixed_effects_noSize, by = "siteID") 

# initialize new columns for mod results
data_mods_fixed_nested$mods <- NA
data_mods_fixed_nested$r2_cond <- NA
data_mods_fixed_nested$AIC <- NA
data_mods_fixed_nested$mods_fixed <- NA
data_mods_fixed_nested$r2_fixed_cond <- NA
data_mods_fixed_nested$r2_fixed_marg <- NA
data_mods_fixed_nested$AIC_w_fixed <- NA
data_mods_fixed_nested$year_coeff_p <- NA

data_mods_fixed_nested$modComp_p <- NA

# loop though and fit models, save results in list
for(i in 1:nrow(data_mods_fixed_nested)){
  message(i)
  message(data_mods_fixed_nested$siteID[[i]])

  tryCatch({
  mod_fixed_result_i <- NA
  mod_fixed_result_i <- fn_fit_models_fixed(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form_fixed[[i]])
   coeff_fixed_result_i<- as.data.frame(coef(summary(mod_fixed_result_i)))
    year_coeff_i <- coeff_fixed_result_i %>% filter(grepl("year", row.names(coeff_fixed_result_i))) 
     data_mods_fixed_nested$year_coeff_p[i] <- year_coeff_i$'Pr(>|t|'
  print(mod_fixed_result_i)
      data_mods_fixed_nested$AIC_w_fixed[i] <- extractAIC(mod_fixed_result_i)[2] # 
    
    if(!performance::check_singularity(mod_fixed_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods_fixed[i] <- list(mod_fixed_result_i)
      
      # get r2 stats
      r2_fixed_i <- performance::r2(mod_fixed_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_fixed_cond[i] <- r2_fixed_i[[1]]["Conditional R2"] #total var explained
      data_mods_fixed_nested$r2_fixed_marg[i] <- r2_fixed_i[[1]]["Marginal R2"] #var explained by fixed effects
    }

  
  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form[[i]])
  print(mod_result_i)
  
      data_mods_fixed_nested$AIC[i] <- extractAIC(mod_result_i)[2] # 
    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods[i] <- list(mod_result_i)
      
      # get r2 stats
      r2_i <- performance::r2(mod_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_cond[i] <- r2_i[[1]]["Conditional R2"] #total var explained
 
    }
  modComp_ANOVA_i <- anova(mod_fixed_result_i, mod_result_i)

  }, error = function(e) {
    message(paste("Error in row", i, ":", e$message))
  })
}

data_mods_fixed_nested$AIC_improve <- data_mods_fixed_nested$AIC - data_mods_fixed_nested$AIC_w_fixed # compare model with moisture variable to model with just random nlcdClass
mod_comp_1core2clip <- data_mods_fixed_nested


mod_comp_1core2clip$yr_eff_lt_05_1core2clip <- if_else(mod_comp_1core2clip$year_coeff_p <0.05, 1, 0, 0)
mod_comp_per_1core2clip <- mod_comp_1core2clip %>% group_by(siteID, sizeCategory) %>% summarise(
                                                                          yr_eff_lt_05_per_1core2clip = 100 * sum(yr_eff_lt_05_1core2clip) / n())



###########################
# data subset 2core1clip, look at year, sizeCategory and interaction in multiyear data for small stature sites

rm(bbc_fixed)
rm(data_mods_fixed_nested)

bbc_fixed <- calc_2core1clip %>% mutate(response_log = log(response))  %>% filter(replicate <= 50) %>% filter(siteID == "DCFS" | siteID == "DSNY" | siteID == "JORN" | siteID == "OAES" | siteID == "ONAQ" | siteID == "WOOD") %>%
  select(-years) # years var was calculated before dropping sites without full design
# site filtering is to small sites with a year effect in raw data
     year_class <- bbc_fixed %>% distinct(eventID, siteID) %>% group_by(siteID) %>% summarise(years = n() )
bbc_fixed <- bbc_fixed %>% left_join(year_class, by = "siteID") %>% filter(years == 2) # %>% filter(nlcd_count > 1) 


data_mods_fixed_nested <- bbc_fixed %>%
  group_by(siteID, sizeCategory,replicate) %>%
  nest() %>%
  left_join(site_formulas_w_fixed_effects_noSize, by = "siteID") 

# initialize new columns for mod results
data_mods_fixed_nested$mods <- NA
data_mods_fixed_nested$r2_cond <- NA
data_mods_fixed_nested$AIC <- NA
data_mods_fixed_nested$mods_fixed <- NA
data_mods_fixed_nested$r2_fixed_cond <- NA
data_mods_fixed_nested$r2_fixed_marg <- NA
data_mods_fixed_nested$AIC_w_fixed <- NA
data_mods_fixed_nested$year_coeff_p <- NA

data_mods_fixed_nested$modComp_p <- NA

# loop though and fit models, save results in list
for(i in 1:nrow(data_mods_fixed_nested)){
  message(i)
  message(data_mods_fixed_nested$siteID[[i]])

  tryCatch({
  mod_fixed_result_i <- NA
  mod_fixed_result_i <- fn_fit_models_fixed(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form_fixed[[i]])
   coeff_fixed_result_i<- as.data.frame(coef(summary(mod_fixed_result_i)))
    year_coeff_i <- coeff_fixed_result_i %>% filter(grepl("year", row.names(coeff_fixed_result_i))) 
     data_mods_fixed_nested$year_coeff_p[i] <- year_coeff_i$'Pr(>|t|'
  print(mod_fixed_result_i)
      data_mods_fixed_nested$AIC_w_fixed[i] <- extractAIC(mod_fixed_result_i)[2] # 
    
    if(!performance::check_singularity(mod_fixed_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods_fixed[i] <- list(mod_fixed_result_i)
      
      # get r2 stats
      r2_fixed_i <- performance::r2(mod_fixed_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_fixed_cond[i] <- r2_fixed_i[[1]]["Conditional R2"] #total var explained
      data_mods_fixed_nested$r2_fixed_marg[i] <- r2_fixed_i[[1]]["Marginal R2"] #var explained by fixed effects
    }

  
  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form[[i]])
  print(mod_result_i)
  
      data_mods_fixed_nested$AIC[i] <- extractAIC(mod_result_i)[2] # 
    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods[i] <- list(mod_result_i)
      
      # get r2 stats
      r2_i <- performance::r2(mod_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_cond[i] <- r2_i[[1]]["Conditional R2"] #total var explained
 
    }
  modComp_ANOVA_i <- anova(mod_fixed_result_i, mod_result_i)

  }, error = function(e) {
    message(paste("Error in row", i, ":", e$message))
  })
}

data_mods_fixed_nested$AIC_improve <- data_mods_fixed_nested$AIC - data_mods_fixed_nested$AIC_w_fixed # compare model with moisture variable to model with just random nlcdClass
mod_comp_2core1clip <- data_mods_fixed_nested


mod_comp_2core1clip$yr_eff_lt_05_2core1clip <- if_else(mod_comp_2core1clip$year_coeff_p <0.05, 1, 0, 0)
mod_comp_per_2core1clip <- mod_comp_2core1clip %>% group_by(siteID, sizeCategory) %>% summarise(
                                                                          yr_eff_lt_05_per_2core1clip = 100 * sum(yr_eff_lt_05_2core1clip) / n())



###########################
# data subset 1core1clip, look at year, sizeCategory and interaction in multiyear data for small AND large stature sites

rm(bbc_fixed)
rm(data_mods_fixed_nested)

bbc_fixed <- calc_1core1clip %>% mutate(response_log = log(response))  %>% filter(replicate <= 50) %>% 
  filter(siteID == "JORN" | siteID == "ONAQ" | siteID == "DCFS" | siteID == "OAES" | siteID == "DEJU" | siteID == "SCBI" | siteID == "TALL" | siteID == "BART" | siteID == "DSNY" | siteID == "WOOD") %>%
  select(-years) # years var was calculated before dropping sites without full design
# site filtering is to small and large sites with a year effect in raw data
     year_class <- bbc_fixed %>% distinct(eventID, siteID) %>% group_by(siteID) %>% summarise(years = n() )
bbc_fixed <- bbc_fixed %>% left_join(year_class, by = "siteID") %>% filter(years == 2) # %>% filter(nlcd_count > 1) 


data_mods_fixed_nested <- bbc_fixed %>%
  group_by(siteID, sizeCategory,replicate) %>%
  nest() %>%
  left_join(site_formulas_w_fixed_effects_noSize, by = "siteID") 

# initialize new columns for mod results
data_mods_fixed_nested$mods <- NA
data_mods_fixed_nested$r2_cond <- NA
data_mods_fixed_nested$AIC <- NA
data_mods_fixed_nested$mods_fixed <- NA
data_mods_fixed_nested$r2_fixed_cond <- NA
data_mods_fixed_nested$r2_fixed_marg <- NA
data_mods_fixed_nested$AIC_w_fixed <- NA
data_mods_fixed_nested$year_coeff_p <- NA

data_mods_fixed_nested$modComp_p <- NA

# loop though and fit models, save results in list
for(i in 1:nrow(data_mods_fixed_nested)){
  message(i)
  message(data_mods_fixed_nested$siteID[[i]])

  tryCatch({
  mod_fixed_result_i <- NA
  mod_fixed_result_i <- fn_fit_models_fixed(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form_fixed[[i]])
   coeff_fixed_result_i<- as.data.frame(coef(summary(mod_fixed_result_i)))
    year_coeff_i <- coeff_fixed_result_i %>% filter(grepl("year", row.names(coeff_fixed_result_i))) 
     data_mods_fixed_nested$year_coeff_p[i] <- year_coeff_i$'Pr(>|t|'
  print(mod_fixed_result_i)
      data_mods_fixed_nested$AIC_w_fixed[i] <- extractAIC(mod_fixed_result_i)[2] # 
    
    if(!performance::check_singularity(mod_fixed_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods_fixed[i] <- list(mod_fixed_result_i)
      
      # get r2 stats
      r2_fixed_i <- performance::r2(mod_fixed_result_i,tolerance = 1e-7)
      data_mods_fixed_nested$r2_fixed_cond[i] <- r2_fixed_i[[1]]["Conditional R2"] #total var explained
      data_mods_fixed_nested$r2_fixed_marg[i] <- r2_fixed_i[[1]]["Marginal R2"] #var explained by fixed effects
    }

  
  mod_result_i <- NA
  mod_result_i <- fn_fit_models(data_x = data_mods_fixed_nested$data[[i]],
                                form_x = data_mods_fixed_nested$form[[i]])
  print(mod_result_i)

    
    if(!performance::check_singularity(mod_result_i, tolerance = 1e-7)){
      data_mods_fixed_nested$mods[i] <- list(mod_result_i)
 
    }
  modComp_ANOVA_i <- anova(mod_fixed_result_i, mod_result_i)

  }, error = function(e) {
    message(paste("Error in row", i, ":", e$message))
  })
}

data_mods_fixed_nested$AIC_improve <- data_mods_fixed_nested$AIC - data_mods_fixed_nested$AIC_w_fixed # compare model with moisture variable to model with just random nlcdClass
mod_comp_1core1clip <- data_mods_fixed_nested


mod_comp_1core1clip$yr_eff_lt_05_1core1clip <- if_else(mod_comp_1core1clip$year_coeff_p <0.05, 1, 0, 0)
mod_comp_per_1core1clip <- mod_comp_1core1clip %>% group_by(siteID, sizeCategory) %>% summarise(
                                                                          yr_eff_lt_05_per_1core1clip = 100 * sum(yr_eff_lt_05_1core1clip) / n())

remove.packages(lmerTest)

#################### combine and summarize difference in AIC between models with N percent (or other fixed variable)
#write_delim(mod_comp_full, paste0("mod_comp_yr_full_",Sys.Date(),".txt"), delim = "\t")

mod_comp_yr_raw <- rbind(mod_comp_2core2clip, mod_comp_1core2clip, mod_comp_per_2core1clip, mod_comp_1core1clip)
#write_delim(mod_comp_yr_raw, paste0("mod_comp_yr_raw.txt_",Sys.Date(),".txt"), delim = "\t")

mod_comp_yr <- mod_comp_per_1core1clip %>% left_join(mod_comp_per_2core1clip, by = c("siteID","sizeCategory")) %>% left_join(mod_comp_per_1core2clip, by = c("siteID","sizeCategory")) %>% left_join(mod_comp_per_2core2clip, by = c("siteID","sizeCategory")) # %>% 
write_delim(mod_comp_yr, paste0("mod_comp_yr_",Sys.Date(),".txt"), delim = "\t")

mod_comp_yr_plot <- mod_comp_yr %>% 
  select(siteID, sizeCategory, yr_eff_lt_05_per_2core2clip, yr_eff_lt_05_per_1core2clip, yr_eff_lt_05_per_2core1clip, yr_eff_lt_05_per_1core1clip) %>% ungroup() %>% 
   rename("4" = "yr_eff_lt_05_per_2core2clip", "2" = "yr_eff_lt_05_per_1core2clip", "2.1" = "yr_eff_lt_05_per_2core1clip", "1" = "yr_eff_lt_05_per_1core1clip") %>%
  left_join(stature_year, by = "siteID")
mod_comp_yr_plot <- mod_comp_yr_plot %>% select(-stature_yr, -years) %>%
  pivot_longer(
    cols = !c(siteID,sizeCategory,stature), 
    names_to = "subsamples", 
    values_to = "percent"
  ) %>% mutate(year2 = 1)
mod_comp_yr_plot <- mod_comp_yr_plot %>% filter(!(stature == "large" & subsamples == "2.1")) %>% filter(!(stature == "small" & (subsamples == "2" | subsamples == "4")))
mod_comp_yr_plot$subsamples <- (ifelse(mod_comp_yr_plot$subsamples == "2.1", "2", mod_comp_yr_plot$subsamples))
mod_comp_yr_plot$year2 <- as.character(mod_comp_yr_plot$year2)

mod_comp_yr_plot$subsamples <- as.integer(mod_comp_yr_plot$subsamples)


mod_comp_full_p <- mod_comp_full %>% select(siteID, sizeCategory, year_coeff_p)
mod_comp_yr_plot <- mod_comp_yr_plot %>% full_join(mod_comp_full_p, by = c("siteID","sizeCategory")) %>% filter(year_coeff_p < 0.05)
#mod_comp_yr_plot$siteID <- factor(mod_comp_yr_plot$siteID, levels=c("BART", "DCFS", "DEJU", "DSNY", "SCBI", "JORN", "TALL", "OAES", "ONAQ", "WOOD"))
mod_comp_yr_plot$siteID <- factor(mod_comp_yr_plot$siteID, levels=c("DCFS", "DSNY", "JORN", "OAES", "ONAQ", "WOOD",   "BART", "DEJU", "SCBI", "TALL"))  

# plot of year effect for sites with multiple years of data
windows(9,6)
mod_comp_yr_plot %>% 
  ggplot(aes(subsamples, percent, color=sizeCategory)) +  
  scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
  geom_point() + 
  geom_line(aes(colour=sizeCategory)) +
  geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) +  # ylim(75, 100) + 
  facet_wrap(vars(siteID, stature), ncol = 4, scales = "free_x") + 
  scale_x_reverse(breaks = c(4,2,1)) + 
  labs(x = "Subsamples per plot") + labs(y = "Percent of resamples in which \n year effect term has p <0.05") + theme(axis.title = element_text(size = 40)) + 
  theme_light() 

savePlot("year effect for sites with multiple years data.pdf", 'pdf')  
