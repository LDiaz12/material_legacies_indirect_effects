# load libraries
library(tidyverse)
library(readr)
library(factoextra)
#library(devtools)
library(ggbiplot)
library(here)
library(mvnormtest)
library(ggpubr)
library(corrplot)
library(plotly)
library(ggfortify)
#source("http://www.sthda.com/upload/rquery_cormat.r")


# import metadata 
#metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
chem_full <- read_csv(here("Data", "Chemistry", "chem_full.csv"))
physio_full <- read_csv(here("Data", "MO24BEAST_physio_metadata.csv")) 
physio_full <- physio_full %>%
  drop_na()
#physio_full <- physio_full %>%
  
# select only chemistry parameters for first PCA plot
metadata_chem <- chem_full %>%
  select(c(DOC_rangemean, DOC_mean, pH_rangemean, pH_mean, TA_rangemean, TA_mean,
           NEC_rangemean, NEC_mean, NEP_rangemean, NEP_mean)) %>%
  dplyr::rename(Dissolved.Organic.Carbon.Range = DOC_rangemean,
                Dissolved.Organic.Carbon.Mean = DOC_mean,
                pH.Range = pH_rangemean, 
                pH.Mean = pH_mean, 
                Total.Alkalinity.Range = TA_rangemean, 
                Total.Alkalinity.Mean = TA_mean, 
                NEC.Range = NEC_rangemean,
                NEC.Mean = NEC_mean, 
                NEP.Range = NEP_rangemean, 
                NEP.Mean = NEP_mean)
                

# select only physiology parameters for second PCA plot
metadata_physio <- physio_full %>%
  select(c(endo_per_cm2, chla.ug.cm2, chlc2.ug.cm2, mean_tissue_biomass, R, NP, GP)) %>%
  dplyr::rename(Endosymbionts = endo_per_cm2, #rename columns for "prettier" names when plotting
         Chlorophyll.a = chla.ug.cm2,
         Chlorophyll.c = chlc2.ug.cm2,
         Tissue.Biomass = mean_tissue_biomass) 

############## PCA OF CHEM DATA ################
# convert to z-scores since remaining variables are not on the same scale 
metadata_chem_scaled <- scale(metadata_chem, scale = TRUE, center = TRUE)
#hist(metadata_chem)

metadata_chem_PCA <- princomp(metadata_chem_scaled, cor=FALSE) # create PCA model using princomp (principal component)
summary(metadata_chem_PCA)

fviz_eig(metadata_chem_PCA) # just a more aesthetic way to look at the scree plot 
metadata_chem_PCA$loadings # shows loadings of each variable on each principal component 
metadata_chem_PCA$scores # gives principal components for each component on each other
chem_pca_data <- as.data.frame(metadata_chem_PCA$scores)
chem_pca_data$TREATMENT <- chem_full$TREATMENT

# rough biplot 
biplot(metadata_chem_PCA, xlab="PC1", ylab="PC2")

# code for a prettier biplot! 
#chem_full$TREATMENT <- factor(chem_full$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

############################
# normality tests
ggplot(chem_full, aes(x=pH_mean)) +
  geom_histogram(bins = 20) # not normally distributed 

ggplot(chem_full, aes(x=pH_rangemean)) +
  geom_histogram(bins = 20) # not normal

ggplot(chem_full, aes(x=TON_mean)) +
  geom_histogram(bins = 125)

ggplot(chem_full, aes(x=TON_rangemean)) +
  geom_histogram(bins = 100)

ggplot(chem_full, aes(x=DOC_mean)) +
  geom_histogram(bins = 125)

ggplot(chem_full, aes(x=DOC_rangemean)) +
  geom_histogram(bins = 125)

ggplot(chem_full, aes(x=TA_mean)) +
  geom_histogram(bins = 40)

ggplot(chem_full, aes(x=TA_rangemean)) +
  geom_histogram(bins = 40)
################################

#plot chem PCA 
p_blank <- autoplot(metadata_chem_PCA, data = chem_full, color = "TREATMENT") +
  coord_fixed(xlim = c(-0.6,0.6), ylim= c(-0.6,0.6)) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated",
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen",
                                                                    "coral","tan")) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = 0.0, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0.0, color = "black", alpha = 0.5)
p_blank

#ggsave(here("Output","PCA", "metadata_PCA_chem_BLANK.png"), plot = p_blank)


metadata_chem_plot2 <- ggbiplot(metadata_chem_PCA, obs.scale=1, var.scale=1, groups=chem_full$TREATMENT, ellipse=TRUE, ellipse.fill = FALSE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  geom_point(aes(colour=factor(chem_full$TREATMENT)), size = 1) + 
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
metadata_chem_plot2

#ggsave(here("Output","PCA", "metadata_PCA_chem_plot.png"), plot = metadata_chem_plot2)


###################################################

############ PCA OF PHYSIOLOGICAL DATA ##################
# convert to z-scores since remaining variables are not on the same scale 
metadata_physio_scaled <- scale(metadata_physio, scale = TRUE, center = TRUE)

metadata_physio_PCA <- princomp(metadata_physio_scaled, cor=FALSE) # create PCA model using princomp (principal component)
summary(metadata_physio_PCA)

fviz_eig(metadata_physio_PCA) # just a more aesthetic way to look at the scree plot 
metadata_physio_PCA$loadings # shows loadings of each variable on each principal component 
metadata_physio_PCA$scores # gives principal components for each component on each other

# rough biplot 
biplot(metadata_physio_PCA, xlab="PC1", ylab="PC2")

# code for a prettier biplot! 
physio_full$TREATMENT <- factor(physio_full$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

metadata_physio_plot <- ggbiplot(metadata_physio_PCA, obs.scale=1, var.scale=1, groups=physio_full$TREATMENT, ellipse=TRUE, ellipse.fill = FALSE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  geom_point(aes(colour=factor(physio_full$TREATMENT)), size = 1) + 
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

ggsave(here("Output", "PCA", "metadata_PCA_physio_plot.png"), plot = metadata_physio_plot)

# Test if ellipses are statistically different from each other # 
## CHEM MANOVA ## 
chem_PCA_matrix <- as_tibble(chem_pca_data)[,c(1:2,11)]
chem_PCA_matrix1 <- as.matrix(chem_PCA_matrix[,1:2])

chem_full <- chem_full %>%
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
metadata_physio_MANOVA <- manova(physio_PCA_scores ~ TREATMENT, data = physio_full)
summary(metadata_physio_MANOVA, test = "Wilks") 
# ellipses in physio PCA are NOT statistically different from each other, which can be seen in the plot and the overlapping 

# qq plots for each physio parameter to check for normality # 
ggqqplot(physio_full, "endo_per_cm2", facet.by="TREATMENT") # normal except for 1 outlier on algae dom
ggqqplot(physio_full, "chla.ug.cm2", facet.by="TREATMENT") # normal with algae and rubble outliers
ggqqplot(physio_full, "chlc2.ug.cm2", facet.by="TREATMENT") # normal with a few rubble and coral outliers
ggqqplot(physio_full, "mean_tissue_biomass", facet.by="TREATMENT") # normal with a rubble outlier 
ggqqplot(physio_full, "R", facet.by="TREATMENT") #normal with 1 control, algae, and coral outlier
ggqqplot(physio_full, "NP", facet.by="TREATMENT") # normal with a couple outliers in each 
ggqqplot(physio_full, "GP", facet.by="TREATMENT") # normal with a few outliers in control and algae dom


# join chem and physio full to plot every single chem and physio parameter with PC1 and PC2
chem_physio <- chem_full %>%
  full_join(physio_full) 

# take chem and physio data frames and plot each variable with PC1 # 
# CHEM VARIABLES # all chem variables should have 4 points per treatment 
# DOC range mean #
DOC_RM_PC1 <- chem_physio %>% # 'DOC rangemean and comp 1
  ggplot(aes(x = Comp.1, y = DOC_rangemean, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
DOC_RM_PC1
#ggsave(here("Output", "PCA", "DOC_rangemean_PC1.png"), plot = DOC_RM_PC1)

# DOC mean #
DOC_mean_PC1 <- chem_physio %>% 
  ggplot(aes(x = Comp.1, y = DOC_mean, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
DOC_mean_PC1 # add on tank numbers to these plots 
#ggsave(here("Output", "PCA", "DOC_mean_PC1.png"), plot = DOC_mean_PC1)

# pH range mean #
pH_RM_PC1 <- chem_physio %>% 
  ggplot(aes(x = Comp.1, y = pH_rangemean, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
pH_RM_PC1
#ggsave(here("Output", "PCA", "pH_rangemean_PC1.png"), plot = pH_RM_PC1)

# pH mean # 
pH_mean_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = pH_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
pH_mean_PC1
#ggsave(here("Output", "PCA", "pH_mean_PC1.png"), plot = pH_mean_PC1)

# TA range mean # 
TA_RM_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = TA_rangemean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
TA_RM_PC1
#ggsave(here("Output", "PCA", "TA_rangemean_PC1.png"), plot = TA_RM_PC1)

# TA mean # 
TA_mean_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = TA_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
TA_mean_PC1
#ggsave(here("Output", "PCA", "TA_mean_PC1.png"), plot = TA_mean_PC1)

# NEC range mean # 
NEC_RM_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = NEC_rangemean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEC_RM_PC1
#ggsave(here("Output", "PCA", "NEC_rangemean_PC1.png"), plot = NEC_RM_PC1)

# NEC Mean # 
NEC_mean_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = NEC_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEC_mean_PC1
#ggsave(here("Output", "PCA", "NEC_mean_PC1.png"), plot = NEC_mean_PC1)

# NEP range mean #
NEP_RM_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = NEP_rangemean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_RM_PC1
#ggsave(here("Output", "PCA", "NEP_rangemean_PC1.png"), plot = NEP_RM_PC1)

# NEP mean # 
NEP_mean_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = NEP_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_mean_PC1
#ggsave(here("Output", "PCA", "NEP_mean_PC1.png"), plot = NEP_mean_PC1)


# PHYSIO VARIABLES # 
# Endosymbionts #
endos_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = endo_per_cm2, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
endos_PC1
#ggsave(here("Output", "PCA", "endos_PC1.png"), plot = endos_PC1)

# Chl a #
chla_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = Chla_norm, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
chla_PC1
#ggsave(here("Output", "PCA", "chla_PC1.png"), plot = chla_PC1)

# Mean Tissue Biomass # 
mean_biomass_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = mean_tissue_biomass, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
mean_biomass_PC1
#ggsave(here("Output", "PCA", "mean_tissue_biomass_PC1.png"), plot = mean_biomass_PC1)

# Respiration Rate #
respiration_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = R, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
respiration_PC1
#ggsave(here("Output", "PCA", "respiration_PC1.png"), plot = respiration_PC1)

# Net Photosynthesis # 
netphoto_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = NP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
netphoto_PC1
#ggsave(here("Output", "PCA", "netphoto_PC1.png"), plot = netphoto_PC1)

# Gross Photo # 
grossphoto_PC1 <- chem_physio %>%
  ggplot(aes(x = Comp.1, y = GP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
grossphoto_PC1
#ggsave(here("Output", "PCA", "grossphoto_PC1.png"), plot = grossphoto_PC1)

## Now we're going to do the same thing with each of the variables but with PC2 ## 

# DOC range mean #
DOC_RM_PC2 <- chem_physio %>% # 'DOC rangemean and comp 2
  ggplot(aes(x = Comp.2, y = DOC_rangemean, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
DOC_RM_PC2
#ggsave(here("Output", "PCA", "DOC_rangemean_PC2.png"), plot = DOC_RM_PC2)

# DOC mean #
DOC_mean_PC2 <- chem_physio %>% 
  ggplot(aes(x = Comp.2, y = DOC_mean, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
DOC_mean_PC2 # add on tank numbers to these plots 
#ggsave(here("Output", "PCA", "DOC_mean_PC2.png"), plot = DOC_mean_PC2)

# pH range mean #
pH_RM_PC2 <- chem_physio %>% 
  ggplot(aes(x = Comp.2, y = pH_rangemean, color = TREATMENT)) +
  geom_point() + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
pH_RM_PC2
#ggsave(here("Output", "PCA", "pH_rangemean_PC2.png"), plot = pH_RM_PC2)

# pH mean # 
pH_mean_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = pH_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
pH_mean_PC2
#ggsave(here("Output", "PCA", "pH_mean_PC2.png"), plot = pH_mean_PC2)

# TA range mean # 
TA_RM_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = TA_rangemean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
TA_RM_PC2
#ggsave(here("Output", "PCA", "TA_rangemean_PC2.png"), plot = TA_RM_PC2)

# TA mean # 
TA_mean_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = TA_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
TA_mean_PC2
#ggsave(here("Output", "PCA", "TA_mean_PC2.png"), plot = TA_mean_PC2)

# NEC range mean # 
NEC_RM_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = NEC_rangemean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEC_RM_PC2
#ggsave(here("Output", "PCA", "NEC_rangemean_PC2.png"), plot = NEC_RM_PC2)

# NEC Mean # 
NEC_mean_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = NEC_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEC_mean_PC2
#ggsave(here("Output", "PCA", "NEC_mean_PC2.png"), plot = NEC_mean_PC2)

# NEP range mean #
NEP_RM_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = NEP_rangemean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_RM_PC2
#ggsave(here("Output", "PCA", "NEP_rangemean_PC2.png"), plot = NEP_RM_PC2)

# NEP mean # 
NEP_mean_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = NEP_mean, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
NEP_mean_PC2
#ggsave(here("Output", "PCA", "NEP_mean_PC2.png"), plot = NEP_mean_PC2)


# PHYSIO VARIABLES # 
# Endosymbionts #
endos_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = endo_per_cm2, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
endos_PC2
#ggsave(here("Output", "PCA", "endos_PC2.png"), plot = endos_PC2)

# Chl a #
chla_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = Chla_norm, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
chla_PC2
#ggsave(here("Output", "PCA", "chla_PC2.png"), plot = chla_PC2)

# Mean Tissue Biomass # 
mean_biomass_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = mean_tissue_biomass, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
mean_biomass_PC2
#ggsave(here("Output", "PCA", "mean_tissue_biomass_PC2.png"), plot = mean_biomass_PC2)

# Respiration Rate #
respiration_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = R, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
respiration_PC2
#ggsave(here("Output", "PCA", "respiration_PC2.png"), plot = respiration_PC2)

# Net Photosynthesis # 
netphoto_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = NP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
netphoto_PC2
#ggsave(here("Output", "PCA", "netphoto_PC2.png"), plot = netphoto_PC2)

# Gross Photo # 
grossphoto_PC2 <- chem_physio %>%
  ggplot(aes(x = Comp.2, y = GP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  theme_bw()
grossphoto_PC2
#ggsave(here("Output", "PCA", "grossphoto_PC2.png"), plot = grossphoto_PC2)











