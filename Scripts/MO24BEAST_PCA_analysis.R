# load libraries
library(tidyverse)
library(readr)
library(factoextra)
library(devtools)
library(ggbiplot)
library(here)
here()
#install.packages("here")

# import metadata 
metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
metadata_chem <- metadata %>%
  select(c(mean_NPOC_mgL, mean_NPOC_uM, mean_TN_mgL, mean_TN_uM, pH_rangemean, pH_mean, TA_rangemean, TA_mean))


# PCA pf chem data 
# convert to z-scores since remaining variables are not on the same scale 
metadata_chem_scaled <- scale(metadata_chem, scale = TRUE, center = TRUE)

metadata_chem_PCA <- princomp(metadata_chem_scaled, cor=FALSE) # create PCA model using princomp (principal component)
summary(metadata_chem_PCA)

fviz_eig(metadata_chem_PCA) # just a more aesthetic way to look at the scree plot 
metadata_chem_PCA$loadings # shows loadings of each variable on each principal component 
metadata_chem_PCA$scores # gives principal components for each component on each other

# rough biplot 
biplot(metadata_chem_PCA, xlab="PC1", ylab="PC2")

# code for a prettier biplot! I think group by treatment
metadata_chem_plot <- ggbiplot(metadata_chem_PCA, obs.scale=1, var.scale=1, groups=metadata$TREATMENT, ellipse=TRUE, varname.size=3, varname.adjust=1.2, circle=FALSE) +
  scale_color_discrete(name='') +
  geom_point(aes(colour=factor(metadata$TREATMENT)), size = 1) + #Color codes by treatment
  theme(legend.direction = 'horizontal', legend.position='bottom', legend.text=element_text(size=8))
metadata_chem_plot

ggsave(here("Output", "metadata_PCA_chem_plot.png"), plot = metadata_chem_plot)


