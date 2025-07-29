#' #######################################################################################
#' Power Analysis Script
#' Originally authored by Eric R Sokol (esokol@battelleecology.org) with input from Andrew L. and 
#' for use with HBP data (last updated for that purpose 2018 Dec 6)
#' Sam Simkin modified for BBC data, starting 2024 Aug 3
#' #######################################################################################

# User defined variables -- change these to modify script

y_var_name = 'response_log'
NUM_SIM = 50
YEARLY_PERCENT_CHANGE = -0.2

###################################
# load packages
library(tidyverse)
library(lme4)

############################
# load linear mixed models fit to empirical data, which used year as a fixed effect
# these will be used as a starting point for the power analysis


########################################################################################

# try all models that converged only for "All herbaceous plants"
if(! exists('data_nested')) load(file='data_nested.rda')
data_mods_nested_to_analyze <- data_nested %>% filter(!is.na(mods)) %>%   #  & sizeCategory == "total"
   mutate(n_boutNumber = 1) # start with previously determined model structure and then impose year effect in this script
data_mods_nested_to_analyze$form <- gsub('response_log ~ ','response_log ~ year + ',data_mods_nested_to_analyze$form)

stature_year <- read.delim("stature_year.txt", header = TRUE, sep = "\t")
data_mods_nested_to_analyze <- data_mods_nested_to_analyze %>% left_join (stature_year, by = "siteID")


#######################################
######################################
# - FUNCTIONS
######################################
##################################
#' Simulate percent cover based on a fit lmer object
#'
#' @param plots_per_site 
#' @param clipIDs per plot 
#' @param years must be two
#' @param bouts_per_year 
#' @param fit_model is a lme4 object 
#' @param year_percent_change percent change between two years, e.g. 0.05 for 5%
#'
#' @return A data.frame with plots, clipID, bout, year, and root biomass
#' @details 
#'
#' @examples

sim_change_data <- function(
  #############################
  # design factors that we want to simulate
  # plots_per_site = 40, 
  # subplots_per_plot = 4, 
  # years = 2, 
  # bouts_per_year = 1, 
  plots_per_site = 30, 
  clipIDs_per_plot = 2,
  cores_per_clip = 2,
  years = 2, 
  bouts_per_year = 1, 
  fit_model = mod,  # the model fit with empirical data -- use this to assign fixed and random effects to factors above
  year_percent_change = -0.2,   # parameter of interest for the power analysis
  y_var_name = y_var_name) {
  
  ###############
  # function code starts here
  
  # make sure there are only two years. The power analysis will be nonsensical with more than 2 years
  stopifnot(years >= 2, years <= 2)
  
  # extract input data from model object
  fit_model_df <- fit_model@frame
  
  # assumes fit_model is a lme4 model -- extract variance and correlation components
  var_components = as.data.frame(lme4::VarCorr(fit_model))
  
  # fixed effects
  coefficients = fixef(fit_model)
  
  if(grepl('(?i)log', y_var_name)){
    y0 <- expm1(coefficients[["(Intercept)"]]) # if log1p (natural log of 1 +x) transformation had been used (I changed input data to be this on 2024-08-02)
    year_effect <-  log1p(y0 * year_percent_change + y0) - log1p(y0)
  }else{
    y0 <- coefficients[["(Intercept)"]]
    year_effect <- y0 * year_percent_change
  }
  
  # get variable names used in model
  rand_var_names <- var_components$grp
  clipID_var_name <- rand_var_names[grepl('^clipID', rand_var_names)]
  plot_var_name <- rand_var_names[grepl('^plot', rand_var_names)]
  

  # extract bout fixed effect -- set to 0 if it's not in the input model
  # the actual bout effect size is unimportant
  if( any(grepl('(?i)bout',names(coefficients))) ){
    bout_effect = coefficients[[which(grepl('(?i)bout',names(coefficients)))]]
  } else {
    bout_effect = 0
  }
  
  # nlcdClass rand effect --  assign to plots in the correct proportion at each site.
  if( 'nlcdClass' %in% rand_var_names ){
    nlcdClass_type_relative_distribution <- fit_model_df %>% group_by(nlcdClass) %>% 
      summarize(
        nlcdClass_num_plots = nrow(.)) %>% 
      ungroup() %>%
      mutate(
        nlcdClass_proportion = nlcdClass_num_plots / (sum(nlcdClass_num_plots))
      )
    
    nlcdClass_sd <- var_components %>% filter(grp == 'nlcdClass') %>% select(sdcor) %>% unlist()
    
    nlcdClass_df <- tibble(
      nlcdClass = nlcdClass_type_relative_distribution$nlcdClass,
      
      # assign random effects to nlcdClasses
      nlcdClass_rand_eff = rnorm(nrow(nlcdClass_type_relative_distribution), 0, nlcdClass_sd)
    )
    
    nlcdClass_plots_df <- tibble(
      plot = seq_len(plots_per_site),
      nlcdClass = sample(
        nlcdClass_type_relative_distribution$nlcdClass,
        size = plots_per_site,
        replace = TRUE,
        prob = nlcdClass_type_relative_distribution$nlcdClass_proportion)) %>% 
      left_join(nlcdClass_df, by = 'nlcdClass') %>%
      select(nlcdClass, nlcdClass_rand_eff, plot)
  }else{
    nlcdClass_plots_df <- tibble(
      nlcdClass = 'generic',
      nlcdClass_rand_eff = 0,
      plot = seq_len(plots_per_site))
  } 
  
  # calculate plot rand effect
if(length(plot_var_name) > 0){
   plot_sd <- var_components %>% filter(grp == plot_var_name) %>% select(sdcor) %>% unlist()
   } else {
    plot_sd <- 0
   } 
  
  plot_df <- tibble(
    plot = seq_len(plots_per_site),
    plot_rand_eff = rnorm(plots_per_site, 0, plot_sd)
  )
  
  # calculate clipID rand effect
  if(length(clipID_var_name) > 0){
    clipID_sd <- var_components %>% filter(grp == clipID_var_name) %>% select(sdcor) %>% unlist()
  } else {
    clipID_sd <- 0
  }
  
  clipID_df <- crossing(
    plot = seq_len(plots_per_site),
    clipID = seq_len(clipIDs_per_plot)
  ) %>% 
    mutate(
      clipID_rand_eff = rnorm(
        plots_per_site * clipIDs_per_plot, 0, clipID_sd)
    )
  
  # calculate residual error
  residual_sd <- var_components %>% filter(grp == 'Residual') %>% select(sdcor) %>% unlist()
  
  # create a simulation scenario for year 0
  sim_df <- crossing(
    plot = seq_len(plots_per_site),
    clipID = seq_len(clipIDs_per_plot),
    year = as.integer(c(0,1)),
    bout = seq_len(bouts_per_year) - 1
  ) %>% 
    mutate(
      residual_error = rnorm(length(plot), 0, residual_sd),
      fixed_effect = coefficients[["(Intercept)"]] + 
        (year_effect * year) +
        (bout_effect * bout)
    ) 
  
  sim_df <- nlcdClass_plots_df %>% 
    left_join(sim_df, by = 'plot')
  
  sim_df <- sim_df %>% 
    left_join(plot_df, by = "plot") %>% 
    left_join(clipID_df, by = c("plot", "clipID")) %>% 
    mutate(
      plot = as.factor(plot),
      clipID = as.factor(clipID),
      bout = as.factor(bout),
      year = as.factor(year))
  
  sim_df <- sim_df %>% mutate(
    y = fixed_effect + nlcdClass_rand_eff + plot_rand_eff + clipID_rand_eff + residual_error
  ) %>% 
    data.table::setnames(c('y','plot','clipID','bout'), c(y_var_name, 'plotID', 'clipID', 'boutNumber'))
  
  return(sim_df)
}

#' Get tval (or zval) for a fitted glmmTMB or lme model for a simulated data set
#'
#' @param x data set for a given model run
#'
#' @return a z value or t value
#' @details Also requires form to be in global env
#' @examples
get_tval_lmer <- function(x, form_updated) {
  tryCatch({
    x <- as.data.frame(x)
    # form updates -- remove terms with only 1 level
    n_plots <- x$plotID %>% unique() %>% length()
    n_clipIDs <- x$clipID %>% unique() %>% length()
    n_bouts <- x$boutNumber %>% unique() %>% length()
    
    # form_updated <- form
    if(! n_bouts > 1) form_updated <- gsub('\\+ boutNumber ', '', form_updated)
    if(! n_clipIDs > 1){
      term_list <- form_updated %>% stringr::str_split('\\+') %>% unlist() %>% trimws()
      form_updated <- term_list[!grepl('clipID', term_list)] %>% paste0(., collapse = ' + ')
    }
    
 #   form_updated <- as.formula(form_updated)
    form_updated <- as.formula(paste(form_updated, collapse = " "))
    sim_fit_mod <- lme4::lmer(formula = form_updated, data = x)
    sim_fit_summary <- summary(sim_fit_mod)$coefficients %>% as.data.frame()
    sim_fit_summary$`t value`[grepl('(?i)year',row.names(sim_fit_summary))]
    
  },
  error=function(err) NA)
}
###############################
##############################
# -- END FUNCTIONS
######################
#####################

# initialize results data frame that script will write to csv
data_power_analysis_results <- data.frame()

# loop for each site/species combination
for(i_row in 1:nrow(data_mods_nested_to_analyze)){
    
  try({
    # get data for i_row_th data set
    data_mods_nested_i <- data_mods_nested_to_analyze[i_row,]
    df <- data_mods_nested_i$data[[1]]
    form <- data_mods_nested_i$form[[1]]
    mod <- data_mods_nested_i$mods[[1]]
    n_bout <- data_mods_nested_i$n_boutNumber[[1]]
    stature <- data_mods_nested_i$stature[[1]]
    # form <- data_mods_nested_i$form[[1]] %>% gsub('\\(1\\|year\\)', 'year', .) %>% as.formula()
    # mod <- glmmTMB(formula = form, data = df, family = tweedie)
    
if(stature == "small"){
design_settings = 
      crossing(
        PLOTS_PER_SITE = c(30,20,15,10),
        clipIDs_per_plot = c(1),
        cores_per_clip = c(2,1),
        NUM_YEARS = 2,
        BOUTS_PER_YEAR = n_bout,
        YEAR_PERCENT = -0.2) %>% 
      filter(!(PLOTS_PER_SITE < 30 & cores_per_clip == 1))
}
else {
design_settings = 
      crossing(
        PLOTS_PER_SITE = c(20,15,10),
        clipIDs_per_plot = c(2,1),
        cores_per_clip = c(2,1),
        NUM_YEARS = 2,
        BOUTS_PER_YEAR = n_bout,
        YEAR_PERCENT = -0.2) %>% 
      filter(!(PLOTS_PER_SITE < 20 & (cores_per_clip == 1 | clipIDs_per_plot ==1))) %>% filter(!(clipIDs_per_plot == 1 & cores_per_clip ==2))
}

design_settings <- design_settings %>% 
      mutate(
        Power = NA_real_,
        n_simulations_tried = NA_integer_,
        n_simulations_used = NA_integer_
      ) %>%
      distinct()
    
    # loop through each design setting, make NUM_SIM replicates of that simulation scenario
    for (i_design in 1:nrow(design_settings)) {
      cat(i_design, "of", nrow(design_settings))
        
      # for design i, make a list of NUM_SIM simulated data sets -- i.e., replication for each design
      sim_dfs_list = lapply(
        seq_len(NUM_SIM),
        function(x)
          sim_change_data(
            plots_per_site = design_settings$PLOTS_PER_SITE[i_design],
            clipIDs_per_plot = design_settings$clipIDs_per_plot[i_design],
            cores_per_clip = design_settings$cores_per_clip[i_design],
            years = design_settings$NUM_YEARS[i_design],
            bouts_per_year = design_settings$BOUTS_PER_YEAR[i_design],
            fit_model = mod,
            year_percent_change = design_settings$YEAR_PERCENT[i_design],
            y_var_name = y_var_name
          )
      )
      
      
        system.time({
          t_vals = vapply(sim_dfs_list, 
                          get_tval_lmer, 
                          pi,
                          form_updated = form)
        })
      
      # Count the number of models fitted to simulated data sets that detect a significant difference among years
      # critical value of 2 is conservative for alpha 0.05 as long as df > 60
      design_settings$Power[i_design] <- mean(abs(t_vals) > 2, na.rm = TRUE) # get absolute value of t-values, see which ones >2 (~ sig) (1 if true) and take mean (fraction of simulations with sig year effect)
      design_settings$n_simulations_tried[i_design] <- length(t_vals)
      design_settings$n_simulations_used[i_design] <- length(na.omit(t_vals))
      cat(" power =", design_settings$Power[i_design], "\n")
      
      data_power_analysis_results_i <- data.frame(
        siteID = data_mods_nested_i$siteID,
        sizeCategory = data_mods_nested_i$sizeCategory,
#        herbGroup = data_mods_nested_i$herbGroup,
        design_settings
      )
    }
    
    data_power_analysis_results <- bind_rows(data_power_analysis_results, data_power_analysis_results_i)  
  }) # END TRY
  print(paste(i_row, 'completed out of', nrow(data_mods_nested_to_analyze)))
} #END LOOP for site/species combination
#data_power_analysis_results <- data_power_analysis_results %>% filter(!(stature))

######################

# write out results to a csv to the local working directory. Add date time to end of filename so you don't accidentally overwrite previous results.


write_file_name <- paste0('RESULTS_power_analysis__',
                          'ALL_convergent__',
                          'NREPS_',NUM_SIM,
                          Sys.time() %>% format('%Y-%m-%d_%H%M%S'),'.csv')
write.csv(data_power_analysis_results, file = write_file_name, row.names = FALSE)











########################################################
## Visualize power analysis

library(tidyverse)
options(stringsAsFactors = FALSE)

# load mods testing for year fixed effect
# if(! exists('data_mods_nested')) load('MODS_for_power_analysis.RDA')
# 
# site_ref_data <- read_csv('all_sites_ref_sheet.csv')

stature_year <- read.delim("stature_year.txt", header = TRUE, sep = "\t")

if(!exists('data_power_analysis_results')){
  # get most recent version of RESULTS_power_analysis
   list_of_files <- list.files()
   list_of_results_files <- list_of_files[grepl('(?i)results_power_analysis', list_of_files)]
 file_to_read <- (sort(list_of_results_files, decreasing = TRUE))[1]
  data_power_analysis_results <- readr::read_csv(file_to_read)
}

windows(9,6) 

### large stature sites

df <- data_power_analysis_results %>% left_join (stature_year, by = "siteID") %>% filter(!is.na(Power)) %>% as.data.frame()

df <- df %>% filter(!(clipIDs_per_plot == 1 & cores_per_clip ==2))
 df$subsamples <- df$clipIDs_per_plot * df$cores_per_clip

df %>% filter(stature == "large" & PLOTS_PER_SITE == 20) %>%
#   ggplot(aes(x=subsamples, y=Power, group=1)) +    # group=1 basically says that there is no group or that everything is in one group, not sure why it doesn't default to this if you don't specify actual grouping variable
   ggplot(aes(x=subsamples, y=Power, color=sizeCategory)) +  
#   geom_line(linetype = "solid") + geom_point() +
   scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
   geom_line(aes(colour=sizeCategory, linetype = sizeCategory)) + geom_point() +
#   geom_line(aes(colour=sizeCategory)) + geom_point() +
#  scale_linetype_manual(values=c("twodash","twodash","twodash","solid")) + 
  scale_linetype_manual(values=c("twodash","twodash","solid","twodash")) + 
  geom_hline(yintercept = 0.80, col="gray", lty=2, lwd=1) +  
  scale_x_reverse(breaks = c(4,2,1)) + 
  facet_wrap(~ siteID, ncol = 8) + 
  labs(x = "Sample design (4 = resampled current (2 clip x 2 core),  2 = 2 clip x 1 core,  1 = 1 clip x 1 core)") + labs(y = "Power") + 
    theme_minimal() + ggtitle("Large stature sites") # + theme(legend.position="none")

 savePlot('power - large stature 2x2 2x1 and 1x1 designs','pdf')

df %>% filter(stature == "large" & subsamples == 4) %>%
   ggplot(aes(x=PLOTS_PER_SITE, y=Power, color=sizeCategory)) +  
   scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
   geom_line(aes(colour=sizeCategory, linetype = sizeCategory)) + geom_point() +
#  scale_linetype_manual(values=c("twodash","twodash","twodash","solid")) + 
  scale_linetype_manual(values=c("twodash","twodash","solid","twodash")) + 
  geom_hline(yintercept = 0.80, col="gray", lty=2, lwd=1) +  
  scale_x_reverse(breaks = c(30,20,15,10)) + 
  facet_wrap(~ siteID, ncol = 8) + 
  labs(x = "Number of plots") + labs(y = "Power") + 
    theme_minimal() + ggtitle("Large stature sites") # + theme(legend.position="none")

 savePlot('power - large stature plot number','pdf')


### small stature sites

df <- data_power_analysis_results %>% left_join (stature_year, by = "siteID") %>% filter(!is.na(Power) & clipIDs_per_plot == 1) %>% as.data.frame()
 df$subsamples <- df$clipIDs_per_plot * df$cores_per_clip

df %>% filter(stature == "small" & PLOTS_PER_SITE == 30) %>% 
   ggplot(aes(x=subsamples, y=Power, color=sizeCategory)) +  
   scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
   geom_line(aes(colour=sizeCategory, linetype = sizeCategory)) +  geom_point(aes(colour=sizeCategory)) +
#  scale_linetype_manual(values=c("twodash","twodash","twodash","solid")) + 
  scale_linetype_manual(values=c("twodash","twodash","solid","twodash")) + 
  geom_hline(yintercept = 0.80, col="gray", lty=2, lwd=1) +  
  scale_x_reverse(breaks = c(2,1)) + 
  facet_wrap(~ siteID, ncol = 8) + 
  labs(x = "Sample design (2 = current (1 clip x 2 core),  1 = 1 clip x 1 core)") + labs(y = "Power") + 
    theme_minimal() + ggtitle("Small stature sites") # + theme(legend.position="none")
 savePlot('power - small stature 1clipx2core and 1x1 designs','pdf')

df %>% filter(stature == "small" & cores_per_clip == 2) %>%
   ggplot(aes(x=PLOTS_PER_SITE, y=Power, color=sizeCategory)) +  
   scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
   geom_line(aes(colour=sizeCategory, linetype = sizeCategory)) + geom_point() +
#  scale_linetype_manual(values=c("twodash","twodash","twodash","solid")) + 
  scale_linetype_manual(values=c("twodash","twodash","solid","twodash")) + 
  geom_hline(yintercept = 0.80, col="gray", lty=2, lwd=1) +  
  scale_x_reverse(breaks = c(30,20,15,10)) + 
  facet_wrap(~ siteID, ncol = 8) + 
  labs(x = "Number of Plots") + labs(y = "Power") + 
    theme_minimal() + ggtitle("Small stature sites") # + theme(legend.position="none")

 savePlot('power - small stature plot number 2 cores per clip','pdf')

graphics.off()


