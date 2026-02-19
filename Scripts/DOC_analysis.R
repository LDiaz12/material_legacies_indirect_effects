library(tidyverse)
library(here)
library(moments)
library(emmeans)
library(agricolae)
library(ggrepel)
library(lubridate)
library(ggplot2)

# read in DOC data files # 
DOC_1 <- read_csv(here("Data", "DOC", "DOC_JUN2_JUN6.csv"))
DOC_2 <- read_csv(here("Data", "DOC", "DOC_JUN6_JUN26.csv"))

# combine both datasets # 
DOC_full <- bind_rows(DOC_1, DOC_2) %>%
  mutate(DATE = mdy(DATE))

# create a datetime column # 
DOC_full <- DOC_full %>%
  mutate(DATETIME = ymd_hms(paste(DATE, TIME)))
#write_csv(DOC_full, here("Data", "DOC", "DOC_full_data.csv")) # all DOC data joined together

DOC_full <- read_csv(here("Data", "DOC", "DOC_full_data.csv"))

## filtering for 12:00 and 21:00 sampling times and reframing to add TA daily mean and daily range between 12 and 9 
DOC_full_summary <- DOC_full %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(DOC_range = NPOC_uM[TIME == hms("12:00:00")] - NPOC_uM[TIME == hms("21:00:00")],
          DOC_dailymean = mean(NPOC_uM, na.rm = TRUE),
          TON_range = TN_uM[TIME == hms("12:00:00")] - TN_uM[TIME == hms("21:00:00")],
          TON_dailymean = mean(TN_uM, na.rm = TRUE))

## create TA plotdata ##
DOC_tankmeans <- DOC_full_summary %>%
  group_by(TANK_NUM) %>%
  summarize(DOC_rangemean = mean(DOC_range, na.rm = TRUE),
            DOC_mean = mean(DOC_dailymean, na.rm = TRUE),
            TON_rangemean = mean(TON_range, na.rm = TRUE), 
            TON_mean = mean(TON_dailymean, na.rm = TRUE))
#write_csv(DOC_tankmeans, here("Data", "DOC", "DOC_tank_means_data.csv")) 

