# load libraries
library(tidyverse)
library(readr)
library(factoextra)
#library(devtools)
library(ggbiplot)
library(here)
#install.packages("here")
#install.packages("mvnormtest")
library(mvnormtest)
library(ggpubr)
library(corrplot)
library(plotly)
library(ggfortify)
#source("http://www.sthda.com/upload/rquery_cormat.r")


# import metadata 
#metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
chem_full <- read_csv(here("Data", "chem_full.csv"))
physio_full <- read_csv(here("Data", "MO24BEAST_physio_metadata.csv"))
physio_full <- physio_full %>%
  
  drop_na()
# select only chemistry parameters for first PCA plot
metadata_chem <- chem_full %>%
  select(c(DOC_rangemean, DOC_mean, TON_rangemean, TON_mean, pH_rangemean, pH_mean, TA_rangemean, TA_mean)) %>%
  dplyr::rename(Dissolved.Organic.Carbon.Range = DOC_rangemean,
                Dissolved.Organic.Carbon.Mean = DOC_mean,
                Total.Organic.Nitrogen.Range = TON_rangemean,
                Total.Organic.Nitrogen.Mean = TON_mean,
                pH.Range = pH_rangemean, 
                pH.Mean = pH_mean, 
                Total.Alkalinity.Range = TA_rangemean, 
                Total.Alkalinity.Mean = TA_mean)
                

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
chem_full$TREATMENT <- factor(chem_full$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

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
ggsave(here("Output","PCA", "metadata_PCA_chem_BLANK.png"), plot = p_blank)

chem_loadings <- metadata_chem_PCA$loadings
scale_factor <- 0.8
chem_loadings_df <- data.frame(x = chem_loadings[,1] * scale_factor - 0.05,
                               y = chem_loadings[,2] * scale_factor + 0.05,
                               label = rownames(chem_loadings))

metadata_chem_plot <- autoplot(metadata_chem_PCA, data = chem_full, color = "TREATMENT", loadings = TRUE, loadings.color = "black") +
  coord_fixed(xlim = c(-0.6,0.6), ylim= c(-0.6,0.6)) +
  labs(color = "Dominant \nBenthic \nCommunity") +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated",
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen",
                                                                    "coral","tan")) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = 0.0, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0.0, color = "black", alpha = 0.5) +
  geom_text(data = chem_loadings_df, aes(x=x, y=y, label = label), size = 2.5, color = "red", 
            angle = -30, fontface = "bold")
metadata_chem_plot

ggsave(here("Output","PCA", "metadata_PCA_chem_plot.png"), plot = metadata_chem_plot)
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
metadata_physio_plot <- ggbiplot(metadata_physio_PCA, obs.scale=1, var.scale=1, groups=physio_full$TREATMENT, ellipse=TRUE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
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
chem_PCA_matrix <- as_tibble(chem_pca_data)[,c(1:2,9)]
chem_PCA_matrix1 <- as.matrix(chem_PCA_matrix[,1:2])

chem_full <- chem_full %>%
  bind_cols(chem_PCA_matrix[,1:2])
metadata_chem_MANOVA <- manova(chem_PCA_matrix1 ~ TREATMENT, data = chem_PCA_matrix)


#chem_PCA_scores <- metadata_chem_PCA$scores
summary(metadata_chem_MANOVA, test = "Wilks", tol=0) #setting tol=0 here means that all PC's are retained for analysis
                                                      #tol=0 means the tolerance for permissible change in eigenvalues is zero
                                                    #useful here bc of the low number of observations. there are 8 variables with only 4 observations for each

#doing qq plots for each variable to check for normality#
ggqqplot(chem_full, "DOC_rangemean", facet.by="TREATMENT")
ggqqplot(chem_full, "DOC_mean", facet.by="TREATMENT")
ggqqplot(chem_full, "TON_rangemean", facet.by="TREATMENT")
ggqqplot(chem_full, "TON_mean", facet.by="TREATMENT")
ggqqplot(chem_full, "pH_rangemean", facet.by="TREATMENT") #pH and TA data especially not normal
ggqqplot(chem_full, "pH_mean", facet.by="TREATMENT") #but most likely due to sampling frequency for TA? 
ggqqplot(chem_full, "TA_rangemean", facet.by="TREATMENT")
ggqqplot(chem_full, "TA_mean", facet.by="TREATMENT")

rquery.cormat(chem_PCA_matrix) #covariance plot
#from the above, TA_rangemean and TA_mean have an r value of 0.89; pH_rangemean and pH_mean have an r value of 0.79. 
#it's possible that if I log transform those variables, they would no longer be collinear, but
#I don't think I can because of negative values from the ranges? 

chem_univariates <- aov(chem_PCA_matrix~TREATMENT, data=chem_full)
summary(chem_univariates)
# is the significance being driven by pH and TA range mean/mean? 

chem_physio <- chem_full %>%
  full_join(physio_full)

chem_physio %>%
  ggplot(aes(x = Comp.1, y = endo_per_cm2, color = TREATMENT)) +
  geom_point()


## PHYSIO MANOVA ## 
physio_PCA_scores <- metadata_physio_PCA$scores
metadata_physio_MANOVA <- manova(physio_PCA_scores ~ TREATMENT, data = physio_full)
summary(metadata_physio_MANOVA, test = "Wilks") 

# qq plots for each physio parameter to check for normality # 
ggqqplot(physio_full, "endo_per_cm2", facet.by="TREATMENT") # normal execpt for 1 outlier on algae dom
ggqqplot(physio_full, "chla.ug.cm2", facet.by="TREATMENT") # normal with algae and rubble outliers
ggqqplot(physio_full, "chlc2.ug.cm2", facet.by="TREATMENT") # normal with a few rubble and coral outliers
ggqqplot(physio_full, "mean_tissue_biomass", facet.by="TREATMENT") # normal with a rubble outlier 
ggqqplot(physio_full, "R", facet.by="TREATMENT") #normal with 1 control, algae, and coral outlier
ggqqplot(physio_full, "NP", facet.by="TREATMENT") # normal with a couple outliers in each 
ggqqplot(physio_full, "GP", facet.by="TREATMENT") # normal with a few outliers in control and algae dom

physio_PCA_matrix <- as.matrix(metadata_physio)
rquery.cormat(physio_PCA_matrix)





