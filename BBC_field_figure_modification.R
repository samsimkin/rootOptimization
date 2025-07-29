## Field optimization


### Variance partitioning (from field script 2)

desired_order <- c(
  "0-2 mm,    >1 nlcd",
  "0-2 mm,    1 nlcd",
  "2-10 mm,    >1 nlcd",
  "2-10 mm,    1 nlcd",
  "total,    >1 nlcd",
  "total,    1 nlcd"
)

data_var_comps_mod <- data_var_comps %>%
  mutate(
    facet_label = interaction(sizeCategoryDisplay, nlcd_multi, sep = ",    "),
    facet_label = factor(facet_label, levels = desired_order)
  )
data_var_comps_mod$comps <- ifelse(data_var_comps_mod$comps == "clipID", "cell", data_var_comps_mod$comps)

# tall stature # Fig 1 -- boxplot -- variance components by site
windows(18, 18, pointsize = 10)
data_var_comps_mod %>% filter(stature == "large") %>% 
ggplot(aes(siteID, iccs, color = comps, fill = comps)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5) +
  facet_wrap(vars(facet_label), ncol = 2, scales = "free_x") +
  geom_col(size = 0.2) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_fill_manual(values = c("deepskyblue", "#009E73", "#E69F00", "white")) +
  xlab('Site ID') +
  ylab(bquote(Total ~ variance ~ (R^2 ~ conditional))) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 8, margin = margin(0.1, 0, 0.1, 0, "mm")),
    plot.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(size = 6)
  )  +   ggtitle("Tall stature sites")

savePlot("TotalVar_keep 1st bout - tall stature 2025.pdf", 'pdf')   # currently Supp Figure 1

# short stature # Fig 2 -- boxplot -- variance components by site
windows(18, 18, pointsize = 10)
data_var_comps_mod %>% filter(stature == "small") %>% 
ggplot(aes(siteID, iccs, color = comps, fill = comps)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5) +
  facet_wrap(vars(facet_label), ncol = 2, scales = "free_x") +
  geom_col(size = 0.2) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_fill_manual(values = c("deepskyblue", "#009E73", "#E69F00", "white")) +
  xlab('Site ID') +
  ylab(bquote(Total ~ variance ~ (R^2 ~ conditional))) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 8, margin = margin(0.1, 0, 0.1, 0, "mm")),
    plot.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(size = 6)
  )  +   ggtitle("Short stature sites")

savePlot("TotalVar_keep 1st bout - short stature 2025.pdf", 'pdf')   # currently Supp Figure 2


#short and tall stature for sites with multiple bouts
desired_order <- c(
  "0-2 mm,    tall",
  "0-2 mm,    short",
  "2-10 mm,    tall",
  "2-10 mm,    short",
  "total,    tall",
  "total,    short"
)

data_var_comps_multiyear_mod <- data_var_comps_multiyear
data_var_comps_multiyear_mod$stature <- ifelse(data_var_comps_multiyear_mod$stature == "large", "tall", data_var_comps_multiyear_mod$stature)
data_var_comps_multiyear_mod$stature <- ifelse(data_var_comps_multiyear_mod$stature == "small", "short", data_var_comps_multiyear_mod$stature)
  
data_var_comps_multiyear_mod <- data_var_comps_multiyear_mod %>%
  mutate(
    facet_label = interaction(sizeCategoryDisplay, stature, sep = ",    "),
    facet_label = factor(facet_label, levels = desired_order)
  )

windows(18, 20, pointsize = 10)
data_var_comps_multiyear_mod %>% 
  ggplot(aes(siteID, iccs, color=comps, fill=comps)) +  
   geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + facet_wrap(vars(facet_label), nrow = 4, scales="free_x") +
  geom_col() +
  scale_color_manual(values=c("black", "black", "black","black","black")) +
  scale_fill_manual(values=c("#009E73","deepskyblue","black","white")) +
  ylab(bquote(Total ~ variance ~ (R^2 ~ conditional))) +
  xlab('Site ID') + 
  theme_light() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1))
 
savePlot("TotalVar_keep both bouts 2025.pdf", 'pdf')




### means replication (from field script 2)

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
  labs(x = "Sample design (4 = resampled current (2 cell x 2 core),  \n 2 = 2 cell x 1 core,  1 = 1 cell x 1 core)") + 
  labs(y = "Percent of resample iterations \n within +/- 10% of observed mean") + 
    theme_minimal() + 
    theme(axis.title = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 12))  + 
    ggtitle("Tall stature sites") # + theme(legend.position="none") 

 savePlot('site_means_tallStature 2025','pdf')  # currently Figure 3


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
  labs(x = "Sample design (2 = current (1 cell x 2 core),  1 = 1 cell x 1 core)") + 
  labs(y = "Percent of resample iterations \nwithin +/- 10% of observed mean") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))  + 
    ggtitle("Short stature sites") # + theme(legend.position="right")

 savePlot('site_means_shortStature 2025','pdf')  # currently Figure 4
 
 
### year comparison (from script 3)

mod_comp_yr_plot_mod <- mod_comp_yr_plot
mod_comp_yr_plot_mod$stature <- ifelse(mod_comp_yr_plot_mod$stature == "large", "tall", mod_comp_yr_plot_mod$stature)
mod_comp_yr_plot_mod$stature <- ifelse(mod_comp_yr_plot_mod$stature == "small", "short", mod_comp_yr_plot_mod$stature)
mod_comp_yr_plot_mod <- mod_comp_yr_plot_mod %>%
  mutate(
    facet_label = interaction(siteID, stature, sep = ",    ")
  )
  
 windows(9,6)
 mod_comp_yr_plot_mod %>% 
  ggplot(aes(subsamples, percent, color=sizeCategory)) +  
#  scale_color_manual(values=c("#009E73","#E69F00","black", "white")) +
  scale_color_manual(values=c("#009E73","black", "#E69F00", "white")) + # to handle fact that 2-10 class is now missing
  geom_point() + 
  geom_line(aes(colour=sizeCategory)) +
  geom_hline(yintercept = 90, col="gray", lty=2, lwd=1) +  # ylim(75, 100) + 
  facet_wrap(vars(facet_label), ncol = 4, scales = "free_x") + 
  scale_x_reverse(breaks = c(4,2,1)) + 
  labs(x = "Subsamples per plot") + 
  labs(y = "Percent of resamples in which \n year effect term has p < 0.05") + 
#   theme_minimal() +
   theme_light() +
  theme(strip.text.x = element_text(size = 12, margin = margin(0.1, 0, 0.1, 0, "mm"))) +
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 14)) 


savePlot("year effect for sites with multiple years data 2025.pdf", 'pdf')  


### power analysis (from script 4)


### large stature sites

df <- data_power_analysis_results %>% left_join (stature_year, by = "siteID") %>% filter(!is.na(Power)) %>% as.data.frame()

df <- df %>% filter(!(clipIDs_per_plot == 1 & cores_per_clip ==2))
 df$subsamples <- df$clipIDs_per_plot * df$cores_per_clip

windows(9,6) 
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
  labs(x = "Sample design (4 = resampled current (2 cell x 2 core), \n  2 = 2 cell x 1 core,  1 = 1 cell x 1 core)") + labs(y = "Power") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))  + 
  ggtitle("Tall stature sites") # + theme(legend.position="none")

 savePlot('power - tall stature 2x2 2x1 and 1x1 designs 2025','pdf')

windows(9,6) 
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
    theme_minimal() + 
  theme(axis.title = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))  + 
  ggtitle("Tall stature sites") # + theme(legend.position="none")

 savePlot('power - tall stature plot number 2025','pdf')


### small stature sites

windows(9,6) 
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
  labs(x = "Sample design (2 = current (1 cell x 2 core),  1 = 1 cell x 1 core)") + labs(y = "Power") + 
    theme_minimal() + 
  theme(axis.title = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))  + 
  ggtitle("Short stature sites") # + theme(legend.position="none")
 savePlot('power - short stature 1clipx2core and 1x1 designs 2025','pdf')

windows(9,6) 
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
    theme_minimal() + 
  theme(axis.title = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))  + 
  ggtitle("Short stature sites") # + theme(legend.position="none")

 savePlot('power - short stature plot number 2 cores per clip 2025','pdf')

graphics.off()