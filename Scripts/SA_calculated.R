library(here)
library(tidyverse)

here()
std_curve <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "Laurel_Standard_Curve.csv"))
surface_area_data <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA.csv"))

# calculating surface area of sphere from known diameters
# area of sphere is 4pi(r)^2
std_curve<- std_curve %>%
  mutate(sphere_radius = (diameter.cm)/2, 
         area = 4 * pi * (sphere_radius)^2)

## plot for linear regression of wax weight and dowel diameter
plot <- std_curve %>%
  ggplot(aes(x = wax_weight.g, y = area)) + 
  geom_point() +
  geom_smooth(method = "lm")
plot

curve_model <- lm(area ~ wax_weight.g, data = std_curve)
curve_model
intercept <- curve_model$coefficients[1] #extracts intercept from coefficients 
slope <- curve_model$coefficients[2] # extracts slope from coefficients 

## calculate surface area of corals 

surface_area_data <- surface_area_data %>%
  mutate(SA_cm_2 = slope*weight_of_wax_g + intercept)

write_csv(surface_area_data, here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))

# cm2 = slope(waxweight)+intercept

