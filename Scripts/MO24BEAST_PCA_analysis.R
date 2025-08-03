# load libraries
library(tidyverse)
library(readr)
library(factoextra)
library(devtools)
library(ggbiplot)
library(here)
library(mvnormtest)
library(ggpubr)
library(corrplot)
library(plotly)
library(ggfortify)
#source("http://www.sthda.com/upload/rquery_cormat.r")
library(vegan)

# import metadata 
#metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
chem_per_tank <- read_csv(here("Data", "Chemistry", "chem_summary_data.csv"))

chem_full <- read_csv(here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))
chem_full <- chem_full %>%
  select(c(TIME, TREATMENT, TA, pH, NPOC_uM)) %>%
  filter(TIME %in% c("12:00:00", "21:00:00")) %>%
  drop_na()

physio_full <- read_csv(here("Data", "MO24BEAST_physio_metadata.csv")) 
physio_full <- physio_full %>%
  drop_na()
  
# select only needed chemistry parameters for first PCA plot
# chem per tank summary data 
metadata_chem <- chem_per_tank %>%
  select(c(DOC_rangemean, DOC_mean, pH_rangemean, pH_mean, TA_rangemean, TA_mean)) %>%
  dplyr::rename(Dissolved.Organic.Carbon.Range = DOC_rangemean,
                Dissolved.Organic.Carbon.Mean = DOC_mean,
                pH.Range = pH_rangemean, 
                pH.Mean = pH_mean, 
                Total.Alkalinity.Range = TA_rangemean, 
                Total.Alkalinity.Mean = TA_mean)

metadata_chem_deltas <- chem_per_tank %>%
  select(c(deltaTA_rangemean, deltaTA_mean, deltapH_rangemean, deltapH_mean, deltaDOC_rangemean, deltaDOC_mean)) %>%
  dplyr::rename(Delta.Total.Alkalinty.Range = deltaTA_rangemean,
                Delta.Total.Alkalinity.Mean = deltaTA_mean,
                Delta.pH.Range = deltapH_rangemean, 
                Delta.pH.Mean = deltapH_mean,
                Delta.Dissolved.Organic.Carbon.Range = deltaDOC_rangemean,
                Delta.Dissolved.Organic.Carbon.Mean = deltaDOC_mean)

metadata_chem_raw <- chem_full %>%
  select(c(TREATMENT, TIME, TA, pH, NPOC_uM)) %>%
  dplyr::rename(Total.Alkalinity = TA, 
                Dissolved.Organic.Carbon = NPOC_uM) %>%
  mutate(log_TA = log(Total.Alkalinity), 
         log_pH = log(pH), 
         log_DOC = log(Dissolved.Organic.Carbon))

# select only physiology parameters for second PCA plot
metadata_physio <- physio_full %>%
  select(c(TREATMENT, endo_per_cm2, chla.ug.cm2, mean_tissue_biomass, R, GP)) %>%
  dplyr::rename(Endosymbionts = endo_per_cm2, #rename columns for "prettier" names when plotting
         Chlorophyll.a = chla.ug.cm2,
         Tissue.Biomass = mean_tissue_biomass) %>%
  mutate(Endosymbionts = log(Endosymbionts), 
         Chlorophyll.a = log(Chlorophyll.a), 
         Tissue.Biomass = log(Tissue.Biomass), 
         R = log(R), 
         GP = log(GP)) %>%
  drop_na()

metadata_physio2 <- metadata_physio %>%
  select(c(Endosymbionts, Chlorophyll.a, Tissue.Biomass, R, GP))

############## PCA OF SUMMARY CHEM DATA ################
# convert to z-scores since remaining variables are not on the same scale 
metadata_chem_scaled <- scale(metadata_chem, scale = TRUE, center = TRUE)

metadata_chem_PCA <- princomp(metadata_chem_scaled, cor=FALSE) # create PCA model using princomp (principal component)
summary(metadata_chem_PCA)

fviz_eig(metadata_chem_PCA) # just a more aesthetic way to look at the scree plot 
metadata_chem_PCA$loadings # shows loadings of each variable on each principal component 
metadata_chem_PCA$scores # gives principal components for each component on each other


# rough biplot 
biplot(metadata_chem_PCA, xlab="PC1", ylab="PC2")

# PCA plot of summary chem values # 
metadata_chem_plot <- ggbiplot(metadata_chem_PCA, obs.scale=1, var.scale=1, groups=chem_per_tank$TREATMENT, ellipse=TRUE, ellipse.fill = FALSE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  geom_point(aes(colour=factor(chem_per_tank$TREATMENT)), size = 1) + 
  labs(color = "Dominant \nBenthic \nCommunity") +
  coord_fixed(xlim = c(-5,5), ylim = c(-5,5)) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen",
                                "Coral_Dom" = "coral","Rubble_Dom" = "tan")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = 0.0, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0.0, color = "black", alpha = 0.5)
metadata_chem_plot

#ggsave(here("Output","PCA", "summary_chem_PCA_plot.png"), plot = metadata_chem_plot)

# PCA plot of delta chem values (delta pH, delta TA, delta DOC) #
chem_deltas_scaled <- scale(metadata_chem_deltas, scale = TRUE, center = TRUE)

chem_deltas_PCA <- princomp(chem_deltas_scaled, cor=FALSE) 
summary(chem_deltas_PCA)
fviz_eig(chem_deltas_PCA)

chem_deltas_plot <- ggbiplot(chem_deltas_PCA, obs.scale=1, var.scale=1, groups=chem_per_tank$TREATMENT, ellipse=TRUE, ellipse.fill = FALSE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  geom_point(aes(colour=factor(chem_per_tank$TREATMENT)), size = 1) + 
  labs(color = "Dominant \nBenthic \nCommunity") +
  coord_fixed(xlim = c(-5,5), ylim = c(-5,5)) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen",
                                "Coral_Dom" = "coral","Rubble_Dom" = "tan")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = 0.0, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0.0, color = "black", alpha = 0.5)
chem_deltas_plot

#ggsave(here("Output", "PCA", "deltas_chem_PCA_plot.png"), plot = chem_deltas_plot)

# PCA plot of raw chem data # 
chem_raw_scaled <- scale(metadata_chem_raw[,-c(1:2)], scale = TRUE, center = TRUE)

chem_raw_PCA <- princomp(chem_raw_scaled[,4:6], cor=FALSE) 
summary(chem_raw_PCA)
fviz_eig(chem_raw_PCA)

chem_data_full <- bind_cols(as_tibble(chem_raw_PCA$scores), metadata_chem_raw)

chem_raw_plot <- ggbiplot(chem_raw_PCA, obs.scale=1, var.scale=1, groups=chem_full$TREATMENT, ellipse=TRUE, ellipse.fill = FALSE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  geom_point(aes(colour=factor(chem_full$TREATMENT)), size = 1) + 
  labs(color = "Dominant \nBenthic \nCommunity") +
  #coord_fixed(xlim = c(-5,5), ylim = c(-5,5)) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen",
                                "Coral_Dom" = "coral","Rubble_Dom" = "tan")) +
  theme_bw() +
  #acet_wrap(~TIME)
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = 0.0, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0.0, color = "black", alpha = 0.5)
chem_raw_plot

#ggsave(here("Output", "PCA", "raw_chem_PCA_plot.png"), plot = chem_raw_plot)


chem_data_full %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = TREATMENT)) + 
  geom_point() + 
  stat_ellipse() + 
  facet_wrap(~TIME)

chem_PCA_matrix <- chem_data_full %>%
  select(Comp.1, Comp.2) %>%
  as.matrix()

chem_per_tank <- chem_per_tank %>%
  bind_cols(chem_PCA_matrix[,1:2])



disper_chem_MANOVA <- betadisper(dist(chem_data_full[,c("log_pH", "log_TA", "log_DOC")]), chem_data_full$TREATMENT, type = "centroid")
disper_chem_MANOVA
anova(disper_chem_MANOVA)


#chem_PCA_scores <- metadata_chem_PCA$scores
summary(metadata_chem_MANOVA, test = "Wilks", tol=0)

###################################################

############ PCA OF PHYSIOLOGICAL DATA ##################
# convert to z-scores since remaining variables are not on the same scale 
metadata_physio_scaled <- scale(metadata_physio2, scale = TRUE, center = TRUE) # physio data is logged

metadata_physio_PCA <- princomp(metadata_physio_scaled, cor=FALSE) # create PCA model using princomp (principal component)
summary(metadata_physio_PCA)

fviz_eig(metadata_physio_PCA) # just a more aesthetic way to look at the scree plot 
metadata_physio_PCA$loadings # shows loadings of each variable on each principal component 
metadata_physio_PCA$scores # gives principal components for each component on each other

# rough biplot 
biplot(metadata_physio_PCA, xlab="PC1", ylab="PC2")

# code for a prettier biplot! 
metadata_physio$TREATMENT <- factor(metadata_physio$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

metadata_physio_plot <- ggbiplot::ggbiplot(metadata_physio_PCA, obs.scale=1, var.scale=1, groups=metadata_physio$TREATMENT, ellipse=TRUE, ellipse.fill = FALSE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  geom_point(aes(color=factor(metadata_physio$TREATMENT)), size = 1) + 
  labs(color = "Dominant \nBenthic \nCommunity") +
  coord_fixed(xlim = c(-5,5), ylim = c(-5,5)) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = 0.0, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0.0, color = "black", alpha = 0.5)
  
metadata_physio_plot

#ggsave(here("Output", "PCA", "metadata_PCA_physio_plot.png"), plot = metadata_physio_plot)

# Test if ellipses are statistically different from each other # 
## CHEM MANOVA ## 
chem_pca_data <- as.data.frame(metadata_chem_PCA$scores)
chem_pca_data$TREATMENT <- chem_per_tank$TREATMENT
chem_PCA_matrix <- as_tibble(chem_pca_data)[,c(1:2,7)]
chem_PCA_matrix1 <- as.matrix(chem_PCA_matrix[,1:2])

chem_per_tank <- chem_per_tank %>%
  bind_cols(chem_PCA_matrix[,1:2])
metadata_chem_MANOVA <- manova(chem_PCA_matrix1 ~ TREATMENT, data = chem_PCA_matrix)


#chem_PCA_scores <- metadata_chem_PCA$scores
summary(metadata_chem_MANOVA, test = "Wilks", tol=0) #setting tol=0 here means that all PC's are retained for analysis
                                                      #tol=0 means the tolerance for permissible change in eigenvalues is zero
                                                    #useful here bc of the low number of observations. there are 8 variables with only 4 observations for each

# each of the ellipses for the chemistry parameters are statistically different from each other, 
# which can be seen from the plot itself 

#################################
#doing qq plots for each variable to check for normality#
ggqqplot(chem_full, "DOC_rangemean", facet.by="TREATMENT")
ggqqplot(chem_full, "DOC_mean", facet.by="TREATMENT")
ggqqplot(chem_full, "TON_rangemean", facet.by="TREATMENT")
ggqqplot(chem_full, "TON_mean", facet.by="TREATMENT")
ggqqplot(chem_full, "pH_rangemean", facet.by="TREATMENT") #pH and TA data especially not normal
ggqqplot(chem_full, "pH_mean", facet.by="TREATMENT") #but most likely due to sampling frequency for TA? 
ggqqplot(chem_full, "TA_rangemean", facet.by="TREATMENT")
ggqqplot(chem_full, "TA_mean", facet.by="TREATMENT")
##############################


## PHYSIO MANOVA ## 
physio_PCA_scores <- metadata_physio_PCA$scores
metadata_physio_MANOVA <- manova(physio_PCA_scores ~ TREATMENT, data = metadata_physio)
summary(metadata_physio_MANOVA, test = "Wilks") 
# ellipses in physio PCA are NOT statistically different from each other

# qq plots for each physio parameter to check for normality # 
ggqqplot(physio_full, "endo_per_cm2", facet.by="TREATMENT") # normal except for 1 outlier on algae dom
ggqqplot(physio_full, "chla.ug.cm2", facet.by="TREATMENT") # normal with algae and rubble outliers
ggqqplot(physio_full, "chlc2.ug.cm2", facet.by="TREATMENT") # normal with a few rubble and coral outliers
ggqqplot(physio_full, "mean_tissue_biomass", facet.by="TREATMENT") # normal with a rubble outlier 
ggqqplot(physio_full, "R", facet.by="TREATMENT") #normal with 1 control, algae, and coral outlier
ggqqplot(physio_full, "NP", facet.by="TREATMENT") # normal with a couple outliers in each 
ggqqplot(physio_full, "GP", facet.by="TREATMENT") # normal with a few outliers in control and algae dom


rquery.cormat(metadata_physio_scaled)



###############################
# join chem and physio full to plot every single chem and physio parameter with PC1 and PC2

# take chem and physio data frames and plot each variable with PC1 # 
# CHEM VARIABLES # all chem variables should have 4 points per treatment 
# NEC range mean and PC1
NEC_PC1 <- chem_per_tank %>% 
  ggplot(aes(x = NEC_rangemean, y = Comp.1)) +
  geom_point() + 
  geom_smooth(aes(x = NEC_rangemean, y = Comp.1)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                #"Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEC_PC1
#ggsave(here("Output", "PCA", "NEC_rangemean_PC1.png"), plot = NEC_PC1)


# NEC mean and PC1
NECmean_PC1 <- chem_per_tank %>% 
  ggplot(aes(x = NEC_mean, y = Comp.1)) +
  geom_point() + 
  geom_smooth(aes(x = NEC_mean, y = Comp.1)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                #"Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NECmean_PC1
#ggsave(here("Output", "PCA", "NEC_mean_PC1.png"), plot = NECmean_PC1)

# NEC range mean and PC2
NEC_PC2 <- chem_per_tank %>% 
  ggplot(aes(x = NEC_rangemean, y = Comp.2, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEC_PC2
#ggsave(here("Output", "PCA", "NEC_rangemean_PC2.png"), plot = NEC_PC2)

# NEC mean and PC2
NECmean_PC2 <- chem_per_tank %>% 
  ggplot(aes(x = NEC_mean, y = Comp.2, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NECmean_PC2
#ggsave(here("Output", "PCA", "NEC_mean_PC2.png"), plot = NECmean_PC2)

# NEP range mean and PC1
NEP_PC1 <- chem_per_tank %>% 
  ggplot(aes(x = NEP_rangemean, y = Comp.1)) +
  geom_point() + 
  geom_smooth(aes(x = NEP_rangemean, y = Comp.1)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                               # "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_PC1
#ggsave(here("Output", "PCA", "NEP_rangemean_PC1.png"), plot = NEP_PC1)

# NEP mean and PC1
NEPmean_PC1 <- chem_per_tank %>% 
  ggplot(aes(x = NEP_mean, y = Comp.1)) +
  geom_point() + 
  geom_smooth(aes(x = NEP_mean, y = Comp.1)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                               # "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEPmean_PC1
#ggsave(here("Output", "PCA", "NEP_mean_PC1.png"), plot = NEPmean_PC1)

NEP_pH_ranges <- chem_per_tank %>% 
  ggplot(aes(x = NEP_rangemean, y = pH_rangemean)) +
  geom_point() + 
  geom_smooth(aes(x = NEP_rangemean, y = pH_rangemean)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
  # "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_pH_ranges
#ggsave(here("Output", "NEP_Plots", "NEP_pH_ranges.png"), plot = NEP_pH_ranges)

NEP_TA_ranges <- chem_per_tank %>% 
  ggplot(aes(x = NEP_rangemean, y = TA_rangemean)) +
  geom_point() + 
  geom_smooth(aes(x = NEP_rangemean, y = TA_rangemean)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
  # "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_TA_ranges
#ggsave(here("Output", "NEP_Plots", "NEP_TA_ranges.png"), plot = NEP_TA_ranges)

NEP_pH_means <- chem_per_tank %>% 
  ggplot(aes(x = NEP_mean, y = pH_mean)) +
  geom_point() + 
  geom_smooth(aes(x = NEP_mean, y = pH_mean)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
  # "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_pH_means
#ggsave(here("Output", "NEP_Plots", "NEP_pH_means.png"), plot = NEP_pH_means)

NEP_TA_means <- chem_per_tank %>% 
  ggplot(aes(x = NEP_mean, y = TA_mean)) +
  geom_point() + 
  geom_smooth(aes(x = NEP_mean, y = TA_mean)) +
  #scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
  # "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_TA_means
#ggsave(here("Output", "NEP_Plots", "NEP_TA_means.png"), plot = NEP_TA_means)


# NEP range mean and PC2 
NEP_PC2 <- chem_per_tank %>% 
  ggplot(aes(x = NEP_rangemean, y = Comp.2, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_PC2 
#ggsave(here("Output", "PCA", "NEP_rangemean_PC2.png"), plot = NEP_PC2)

# NEP mean and PC2  
NEPmean_PC2 <- chem_per_tank %>% 
  ggplot(aes(x = NEP_mean, y = Comp.2, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEPmean_PC2
#ggsave(here("Output", "PCA", "NEP_mean_PC2.png"), plot = NEPmean_PC2)


