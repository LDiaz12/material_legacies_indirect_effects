library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(ggridges)
library(agricolae)
library(lme4)
library(lmerTest)
library(moments)
library(performance)
#install.packages("dplyr")

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Chemistry", "TrisCalSummer2024.csv"))
pHData<-read_csv(here("Data", "Chemistry", "CarbonateChemistry.csv"))
TableID<-read_csv(here("Data", "TableID.csv"))


pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = pHcalib))) %>% # linear regression of mV and temp of the tris
  reframe(broom::tidy(fitpH)) %>% 
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TEMPINLAB*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TEMPINSITU)%>%
  drop_na(mV) %>%
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=SALINITY,T=TEMPINLAB))  # calculate pH of the samples using the pH seacarb function

#Now calculate pH
pHSlope <-pHSlope%>%
  drop_na(TEMPINSITU, TEMPINLAB, SALINITY, pH) %>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = 2200, Tinsi = TEMPINSITU, Tlab = TEMPINLAB, 
                            S = SALINITY, Pt = 0.1, k1k2 = "m10", kf = "dg")) %>%
  select(!pH)%>%
  rename(pH = pH_insitu) %>% 
  ungroup() %>%
  select(-c(mV, TrisCalDate, TTris, `(Intercept)`, mVTris)) ## warnings are fine, ignore them

# remove the inflow data and join it with the tanks that had that specific inflow water
InflowData <- pHSlope %>%
  filter(TANK_NUM %in% c("Inflow1","Inflow2")) %>%
  select(-c(FLOW_LEFT, FLOW_RIGHT, Notes, DO_MG_L, SALINITY, TEMPINSITU))  %>% ### remove the values that I don't need -- You will eventually need to keep TA which is why I dropped these instead of coding for the ones that I need
  rename(pH_inflow = pH,
         TA_inflow = TA) %>%# rename the pH to show that it is inflow pH
  mutate(INFLOW_TABLE = ifelse(TANK_NUM == "Inflow1",1,2)) %>% # give them inflow numbers to pair easily with the TankID 
  ungroup()%>%
  select(DATE,TIME, INFLOW_TABLE, pH_inflow, TA_inflow) # drop the Tank ID column to be able to join the data correctly by inflow #


SurfaceArea <- 22.5*22.5

Data<-pHSlope %>%
  ungroup()%>%
  filter(!TANK_NUM %in% c("Inflow1","Inflow2"))%>% # filter out the inflow data
  mutate(TANK_NUM = as.numeric(TANK_NUM))%>% # convert to numeric since the inflow data is now dropped
  left_join(TableID) %>%
  left_join(InflowData) %>% # join with the inflow data for easier calculations of rates
  mutate(DATETIME = ymd_hms(paste(DATE,TIME)), # make a datetime
         pHDiff = pH - pH_inflow, # calculate the difference between the inflow and the pH in each tank 
         totalflow = FLOW_RIGHT+FLOW_LEFT,
         residence_time = (1/totalflow)*(10000/60),# convert ml/min to hours by multiplying by the volume of water in ml and divide by 60
         deltaTA = TA_inflow - TA, # calculate the difference between in and outflow
         NEC = (deltaTA/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea) ### for a real rate should probably normalize the delta TA to the delta control just like in respo
         )
#write_csv(Data, here("Data", "Chemistry", "Full_Chem_Data.csv"))
### Now Make a plot showing how the Tank pH differed from the inflow pH over time

tank_pH_diffs <- Data %>% ## something weird going on with this plot... 
  ggplot(aes(x = DATETIME, y = pHDiff, color = TREATMENT, group = TANK_NUM, na.rm = TRUE))+
  geom_point()+
  geom_line()
tank_pH_diffs +
  scale_color_hue(labels = c("Algae-Dominated", "Control", "Rubble-Dominated", "Coral-Dominated"))

#ggsave(plot = tank_pH_diffs, filename = here("Output", "tank_pH_diffs.png"), width = 11, height = 9)


# pull out the delta TA from the controls and take average by inflow and normalize the NEC rates to it-- this accounts for changes in TA due to background water
control_deltaTA<-Data %>%
  filter(TREATMENT == "Control") %>%
  select(DATETIME, deltaTA, INFLOW_TABLE ) %>%
  group_by(DATETIME, INFLOW_TABLE)%>%
  summarise(deltaTA_blank = mean(deltaTA, na.rm = TRUE)) 

Data<-Data %>%
  left_join(control_deltaTA) %>%
  mutate(NEC = ((deltaTA-deltaTA_blank)/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea) ### for a real rate should probably normalize the delta TA to the delta control just like in respo
  )

### Select what I want for TA data ### 
TA_Data <- Data %>%
  select(DATETIME, DATE, TIME, TA, TREATMENT, TANK_NUM) %>%
  group_by(DATETIME, DATE, TIME, TREATMENT, TANK_NUM) %>%
  drop_na()
TA_Data

## create rect intervals for the NEC plot below for light vs dark times ##
rect_intervals <- tibble::tibble(
  xmin = as.POSIXct(c("2024-06-02 06:00:00", "2024-06-02 18:00:00", 
                      "2024-06-03 06:00:00", "2024-06-03 18:00:00", 
                      "2024-06-04 06:00:00", "2024-06-04 18:00:00",
                      "2024-06-05 06:00:00", "2024-06-05 18:00:00",
                      "2024-06-06 06:00:00", "2024-06-06 18:00:00",
                      "2024-06-07 06:00:00", "2024-06-07 18:00:00",
                      "2024-06-08 06:00:00", "2024-06-08 18:00:00",
                      "2024-06-09 06:00:00", "2024-06-09 18:00:00",
                      "2024-06-10 06:00:00", "2024-06-10 18:00:00",
                      "2024-06-11 06:00:00", "2024-06-11 18:00:00",
                      "2024-06-12 06:00:00", "2024-06-12 18:00:00",
                      "2024-06-13 06:00:00", "2024-06-13 18:00:00",
                      "2024-06-14 06:00:00", "2024-06-14 18:00:00",
                      "2024-06-15 06:00:00", "2024-06-15 18:00:00",
                      "2024-06-16 06:00:00", "2024-06-16 18:00:00",
                      "2024-06-17 06:00:00", "2024-06-17 18:00:00",
                      "2024-06-18 06:00:00", "2024-06-18 18:00:00",
                      "2024-06-19 06:00:00", "2024-06-19 18:00:00",
                      "2024-06-20 06:00:00", "2024-06-20 18:00:00")),
  xmax = as.POSIXct(c("2024-06-02 18:00:00", "2024-06-03 06:00:00", 
                      "2024-06-03 18:00:00", "2024-06-04 06:00:00", 
                      "2024-06-04 18:00:00", "2024-06-05 06:00:00",
                      "2024-06-05 18:00:00", "2024-06-06 06:00:00",
                      "2024-06-06 18:00:00", "2024-06-07 06:00:00",
                      "2024-06-07 18:00:00", "2024-06-08 06:00:00",
                      "2024-06-08 18:00:00", "2024-06-09 06:00:00",
                      "2024-06-09 18:00:00", "2024-06-10 06:00:00",
                      "2024-06-10 18:00:00", "2024-06-11 06:00:00",
                      "2024-06-11 18:00:00", "2024-06-12 06:00:00",
                      "2024-06-12 18:00:00", "2024-06-13 06:00:00",
                      "2024-06-13 18:00:00", "2024-06-14 06:00:00",
                      "2024-06-14 18:00:00", "2024-06-15 06:00:00",
                      "2024-06-15 18:00:00", "2024-06-16 06:00:00",
                      "2024-06-16 18:00:00", "2024-06-17 06:00:00",
                      "2024-06-17 18:00:00", "2024-06-18 06:00:00",
                      "2024-06-18 18:00:00", "2024-06-19 06:00:00",
                      "2024-06-19 18:00:00", "2024-06-20 06:00:00",
                      "2024-06-20 18:00:00", "2024-06-21 06:00:00")),
  fill = c("lightyellow", "lightgrey", "lightyellow", "lightgrey", "lightyellow",
           "lightgrey", "lightyellow", "lightgrey", "lightyellow", "lightgrey", 
           "lightyellow", "lightgrey", "lightyellow","lightgrey", "lightyellow", 
           "lightgrey", "lightyellow","lightgrey", "lightyellow", "lightgrey", 
           "lightyellow","lightgrey", "lightyellow", "lightgrey", "lightyellow",
           "lightgrey", "lightyellow", "lightgrey", "lightyellow","lightgrey", 
           "lightyellow", "lightgrey", "lightyellow","lightgrey", "lightyellow", 
           "lightgrey", "lightyellow","lightgrey"))

### Plot TA ###
TA_plot <- TA_Data %>%
  group_by(TREATMENT, DATETIME)%>%
  summarise(mean_TA = mean(TA, na.rm = TRUE),
            se_TA = sd(TA, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = DATETIME, y = mean_TA, color = TREATMENT, na.rm = TRUE))+
  geom_rect(data = rect_intervals, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 1/5, color = NA, inherit.aes = FALSE) + 
  scale_fill_identity() +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_TA - se_TA,
                    ymax = mean_TA + se_TA), color = "black", width = 0.1) +
  theme_classic() +
  labs(x="Date & Time",
       y = "Mean Total Alkalinity (umol/kg)") +
  geom_line() +
  theme(plot.title = element_text(size = 14))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
TA_plot

## filtering for 12:00 and 21:00 sampling times and reframing to add TA daily mean and daily range between 12 and 9 
TA_data2 <- TA_Data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(TA_range = TA[TIME == hms("12:00:00")] - TA[TIME == hms("21:00:00")],
          TA_dailymean = mean(TA, na.rm = TRUE))

## create TA plotdata ##
TA_plotdata <- TA_data2 %>%
  group_by(TREATMENT) %>%
  summarize(TA_rangemean = mean(TA_range, na.rm = TRUE), # means and ranges throughout entire 25-day experimental period
            TA_rangese = sd(TA_range, na.rm = TRUE)/sqrt(n()),
            TA_mean = mean(TA_dailymean, na.rm = TRUE),
            TA_se = sd(TA_dailymean, na.rm = TRUE)/sqrt(n()))
TA_plotdata

#write_csv(TA_plotdata, here("Data", "Chemistry", "Cleaned_TA_Data_per_Treatment.csv"))


# filtering out outlier
TA_data2_filtered <- TA_data2 %>%
  filter(TA_range < 400)

## plot daily TA range ## 
TA_plotdata2 <- TA_data2_filtered %>%
  group_by(TREATMENT) %>%
  summarize(TA_rangemean = mean(TA_range, na.rm = TRUE), # means and ranges throughout entire 25-day experimental period
            TA_rangese = sd(TA_range, na.rm = TRUE)/sqrt(n()),
            TA_mean = mean(TA_dailymean, na.rm = TRUE),
            TA_se = sd(TA_dailymean, na.rm = TRUE)/sqrt(n()))
TA_plotdata2

TA_plotdata2$TREATMENT <- factor(TA_plotdata2$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

TA_range_plot2 <- TA_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = TA_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean Total Alkalinity Range" ~ (µmol ~ kg^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = TA_data2_filtered, aes(x = TREATMENT, y = TA_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = TA_rangemean - TA_rangese,
                    ymax = TA_rangemean + TA_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
TA_range_plot2
#ggsave(plot = TA_range_plot2, filename = here("Output", "TA_range_plot2.png"), width = 9, height = 7)

# mean range TA stats #
TA_range_model <- lmer(TA_range ~ TREATMENT +(1|TANK_NUM), data=TA_data2)
check_model(TA_range_model) #not great homogeneity of variance; potential outlier point 33 
summary(TA_range_model)
anova(TA_range_model)

# TA per treatment
TA_treatmentdata <- TA_data2 %>%
  group_by(TREATMENT) %>%
  summarize(TA_rangemean = mean(TA_range, na.rm = TRUE),
            TA_rangese = sd(TA_range, na.rm = TRUE)/sqrt(n()),
            TA_mean = mean(TA_dailymean, na.rm = TRUE),
            TA_se = sd(TA_dailymean, na.rm = TRUE)/sqrt(n()))


# create new plot with the outlier removed - ask Nyssa if there's a cleaner way to do this
TA_range_plot <- TA_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = TA_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean Total Alkalinity Range" ~ (µmol ~ kg^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = TA_data2_filtered, aes(x = TREATMENT, y = TA_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = TA_rangemean - TA_rangese,
                    ymax = TA_rangemean + TA_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
TA_range_plot
#ggsave(plot = TA_range_plot, filename = here("Output", "TA_range_plot.png"), width = 9, height = 7)


# stats with outlier filtered out
TA_range_model2 <- lmer(TA_range ~ TREATMENT + (1|TANK_NUM), data = TA_data2_filtered)
check_model(TA_range_model2) # model fits much better and meets assumptions with outlier removed
summary(TA_range_model2) # now seeing significant effect of coral dom communities on TA daily range 
anova(TA_range_model2) # significant effect of community type on daily range in TA

# mean TA plot #
TA_mean_plot <- TA_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = TA_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean TA" ~ (µmol ~ kg^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = TA_data2_filtered, aes(x = TREATMENT, y = TA_dailymean), alpha = 0.7) +
  geom_errorbar(aes(ymin = TA_mean - TA_se,
                    ymax = TA_mean + TA_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
TA_mean_plot
ggsave(plot = TA_mean_plot, filename = here("Output", "TA_mean_plot.png"), width = 9, height = 6)

## TA daily mean stats ##
TA_mean_model <- lmer(TA_dailymean ~ TREATMENT + (1|TANK_NUM), data=TA_data2_filtered)
check_model(TA_mean_model)
summary(TA_mean_model) # rubble dom and coral dom community significant at p < 0.05 for both
anova(TA_mean_model) # significant effect of benthic community on mean TA throughout experimental period 
# this makes sense given the high calcification rates of corals and CCA

# model without random effects for post hoc groupings
TA_mean_model2 <- lm(TA_dailymean ~ TREATMENT, data=TA_data2_filtered)
HSD.test(TA_mean_model2, "TREATMENT", console = TRUE)
# coral and rubble dom communities are most similar, algae is most similar to controls which makes 
# sense bc turbs are non calcifying 

###################################
### NEC DATA ### 
NEC_data <- Data %>%
  select(DATETIME, DATE, TIME, NEC, TREATMENT, TANK_NUM) %>%
  group_by(DATETIME, DATE, TIME, TREATMENT, TANK_NUM) %>%
  drop_na()

## plot mean NEC over course of experiment and color by light vs dark time intervals ## 
NEC_plot <- NEC_data %>%
  group_by(TREATMENT, DATETIME)%>%
  summarise(mean_NEC = mean(NEC, na.rm = TRUE),
            se_NEC = sd(NEC, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = DATETIME, y = mean_NEC, color = TREATMENT, na.rm = TRUE))+
  geom_rect(data = rect_intervals, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 1/5, color = NA, inherit.aes = FALSE) + 
  scale_fill_identity() +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_NEC - se_NEC,
                    ymax = mean_NEC + se_NEC), color = "black", width = 0.1) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  labs(x="Date & Time",
       y = "Mean Net Ecosystem Calcification (NEC) (mmol/m2h)") +
  geom_line() +
  theme(plot.title = element_text(size = 14))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_plot

## filter NEC data for 12:00 and 21:00 sampling and reframe for range and daily mean of these times ## 
NEC_data2 <- NEC_data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(NEC_range = NEC[TIME == hms("12:00:00")] - NEC[TIME == hms("21:00:00")],
          NEC_dailymean = mean(NEC, na.rm = TRUE))

NEC_plotdata <- NEC_data2 %>%
  group_by(TREATMENT, TANK_NUM) %>%
  dplyr::summarize(NEC_rangemean = mean(NEC_range, na.rm = TRUE),
            NEC_rangese = sd(NEC_range, na.rm = TRUE)/sqrt(n()),
            NEC_mean = mean(NEC_dailymean, na.rm = TRUE),
            NEC_se = sd(NEC_dailymean, na.rm = TRUE)/sqrt(n()))
NEC_plotdata

## plot daily NEC range ## 
NEC_range_plot <- NEC_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEC_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean NEC Range" ~ (mmol ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEC_data2, aes(x = TREATMENT, y = NEC_range), alpha = 0.7) +
  #geom_errorbar(aes(ymin = NEC_rangemean - NEC_rangese,
                    #ymax = NEC_rangemean + NEC_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_range_plot # crazy outlier, let's filter out 
#ggsave(plot = NEC_range_plot, filename = here("Output", "NEC_range_plot.png"), width = 9, height = 6)

# filter out outlier from coral dom treatment 
NEC_data_2_filtered <- NEC_data2 %>%
  filter(NEC_range > -8, 
         NEC_dailymean < 0.8)

# new plotdata calculations 
NEC_plotdata2 <- NEC_data_2_filtered %>%
  group_by(TREATMENT) %>%
  summarize(NEC_rangemean = mean(NEC_range, na.rm = TRUE),
            NEC_rangese = sd(NEC_range, na.rm = TRUE)/sqrt(n()),
            NEC_mean = mean(NEC_dailymean, na.rm = TRUE),
            NEC_se = sd(NEC_dailymean, na.rm = TRUE)/sqrt(n()))
NEC_plotdata2

# new plot using filtered data with outlier removed 
NEC_range_plot2 <- NEC_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = NEC_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean NEC Range" ~ (mmol ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEC_data_2_filtered, aes(x = TREATMENT, y = NEC_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEC_rangemean - NEC_rangese,
                    ymax = NEC_rangemean + NEC_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_range_plot2
ggsave(plot = NEC_range_plot2, filename = here("Output", "NEC_range_plot2.png"), width = 9, height = 6)


## NEC mean daily range stats with outlier removed ## 
NEC_range_model <- lmer(NEC_range ~ TREATMENT +(1|TANK_NUM), data=NEC_data_2_filtered)
check_model(NEC_range_model) # not great homogeneity of variance 
summary(NEC_range_model)
anova(NEC_range_model) # no significant effect of community type on daily range in NEC 

## plot NEC daily mean ##
NEC_plotdata2$TREATMENT <- factor(NEC_plotdata2$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_mean_plot <- NEC_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = NEC_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean NEC" ~ (mmol ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEC_data_2_filtered, aes(x = TREATMENT, y = NEC_dailymean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEC_mean - NEC_se,
                    ymax = NEC_mean + NEC_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_mean_plot # looks like an outlier... 
ggsave(plot = NEC_mean_plot, filename = here("Output", "NEC_mean_plot.png"), width = 9, height = 6)

## NEC daily mean stats ##
NEC_mean_model <- lmer(NEC_dailymean ~ TREATMENT + (1|TANK_NUM), data=NEC_data_2_filtered)
check_model(NEC_mean_model)
summary(NEC_mean_model) # significant effect of coral and rubble dom communities on the daily mean NEC
anova(NEC_mean_model)# significant effect of community type of daily mean NEC 


######################################
### PH DATA ###
pH_plot <- Data %>%
  group_by(TREATMENT, DATETIME)%>%
  summarise(mean_diff = mean(pHDiff, na.rm = TRUE),
            se_diff = sd(pHDiff, na.rm = TRUE)/sqrt(n()))%>%
  filter(!is.na(mean_diff), !is.na(se_diff), !is.na(DATETIME)) %>%
  ggplot(aes(x = DATETIME, y = mean_diff, color = TREATMENT, na.rm = TRUE))+
  geom_rect(data = rect_intervals, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 1/5, color = NA, inherit.aes = FALSE) + 
  scale_fill_identity() +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff+se_diff), color = "black", width = 0.1) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  labs(x="Date & Time",
       y = "Change in pH") +
  geom_line() +
  annotate("text", x = ymd_hms("2024-06-03 08:00:00"), y = 0.1, label = "Overcast", size = 5) +
  annotate("text", x = ymd_hms("2024-06-06 13:00:00"), y = 0.1, label = "Overcast \n Rain", size = 5) +
  theme(plot.title = element_text(size = 14))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_plot

# create columns for mean DAY pH and mean NIGHT pH # 
pH_day_night <- Data %>% 
  filter(TIME %in% c("12:00:00", "21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(pH_day = pH[TIME == hms("12:00:00")], 
          pH_night = pH[TIME == hms("21:00:00")])
#write_csv(pH_day_night, here("Data", "Chemistry", "pH_day_night.csv"))

# now calculate mean pH for both day and night sampling times # 
pH_day_night_means <- pH_day_night %>% 
  group_by(TREATMENT) %>%
  summarize(pH_day_mean = mean(pH_day, na.rm = TRUE), 
            pH_day_se = sd(pH_day, na.rm = TRUE)/sqrt(n()), 
            pH_night_mean = mean(pH_night, na.rm = TRUE), 
            pH_night_se = sd(pH_night, na.rm = TRUE)/sqrt(n()))
#write_csv(pH_day_night_means, here("Data", "Chemistry", "pH_day_night_means.csv"))


pHday_means_plot <- ggplot(data = pH_day_night_means, aes(x = TREATMENT, y = pH_day_mean)) +
  geom_jitter()
pHday_means_plot

pHnight_means_plot <- ggplot(data = pH_day_night_means, aes(x = TREATMENT, y = pH_night_mean)) +
  geom_jitter()
pHnight_means_plot

## filter only 12:00 and 21:00 sampling times ## 
pH_clean <- Data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(pH_range = pH[TIME == hms("12:00:00")] - pH[TIME == hms("21:00:00")],
            pH_dailymean = mean(pH, na.rm = TRUE))
pH_clean

# add pH day and pH night means to pH_clean 
pH_clean <- pH_clean %>% 
  right_join(pH_day_night_means)

#write_csv(pH_clean, here("Data", "Chemistry", "Cleaned_pH_Data_FULL.csv"))

pH_plotdata<- pH_clean %>%
  group_by(TREATMENT) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))
pH_plotdata
#write_csv(pH_plotdata, here("Data", "Chemistry", "Cleaned_pH_Data_per_Treatment.csv"))

pH_plotdata_tank<- pH_clean %>%
  group_by(TANK_NUM) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))
pH_plotdata_tank

write_csv(pH_plotdata_tank_full, here("Data", "Chemistry", "Cleaned_pH_Data_per_Tank.csv"))


## plot pH range from 12:00 and 21:00 sampling throughout experiment ## 
pH_plotdata$TREATMENT <- factor(pH_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

pH_range_plot <- pH_plotdata %>%
  ggplot(aes(x = TREATMENT, y = pH_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = "Daily Mean pH Range") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = pH_clean, aes(x = TREATMENT, y = pH_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = pH_rangemean - pH_rangese,
                    ymax = pH_rangemean + pH_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_range_plot
ggsave(plot = pH_range_plot, filename = here("Output", "pH_range_plot.png"), width = 9, height = 6)

## pH range stats ##
pH_range_model <- lmer(pH_range ~ TREATMENT +(1|TANK_NUM), data=pH_clean)
check_model(pH_range_model)
summary(pH_range_model) # significant effect of control treatment on daily mean pH but not what we were 
# looking for. expected to see differences between the communities... 
anova(pH_range_model)

# take outlier from rubble dom community out and replot/redo stats 
pH_clean_filtered <- pH_clean %>%
  filter(pH_range > -0.2)

pH_plotdata2 <- pH_clean_filtered %>%
  group_by(TREATMENT) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))
pH_plotdata2

pH_plotdata2$TREATMENT <- factor(pH_plotdata2$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

pH_range_plot2 <- pH_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = pH_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = "Daily Mean pH Range") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = pH_clean_filtered, aes(x = TREATMENT, y = pH_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = pH_rangemean - pH_rangese,
                    ymax = pH_rangemean + pH_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_range_plot2

# stats with filtered data set
pH_range_model2 <- lmer(pH_range ~ TREATMENT +(1|TANK_NUM), data=pH_clean_filtered)
check_model(pH_range_model2)
summary(pH_range_model2)
anova(pH_range_model2) # control still only the significant treatment on mean pH range 

## mean pH plot of 12:00 and 21:00 sampling throughout experiment ## 
pH_mean_plot <- pH_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = pH_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = "Daily Mean pH") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = pH_clean_filtered, aes(x = TREATMENT, y = pH_dailymean), alpha = 0.7) +
  geom_errorbar(aes(ymin = pH_mean - pH_se,
                    ymax = pH_mean + pH_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_mean_plot
ggsave(plot = pH_mean_plot, filename = here("Output", "pH_mean_plot.png"), width = 9, height = 6)

## mean pH treatment stats ## 
mean_pH_model <- lmer(pH_dailymean~TREATMENT +(1|TANK_NUM), data=pH_clean_filtered)
check_model(mean_pH_model)
summary(mean_pH_model)
anova(mean_pH_model)
# no significant effect of community type on daily mean pH 

### Combine biomass data and carb chem data ### 
# read in biomass (AFDW) data 
afdw <- read_csv(here("Data", "Data_Raw", "Growth", "MO24BEAST_AFDW.csv"))

#filter out 'pre' experiment treatments 
afdw_nopre <- afdw %>%
  filter(!TREATMENT == "Pre")

chem_biomass_data <- afdw_nopre %>% # join afdw data frame with chem data 
  right_join(Data) %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, TANK_NUM, 
         TA, deltaTA, NEC, mean_AFDW, mean_tissue_biomass, DATE, TIME, pH)

chem_biomass_data_clean <- chem_biomass_data %>%
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TIME, TANK_NUM) 

chem_biomass_plotdata <- chem_biomass_data_clean %>%
  group_by(TANK_NUM, TREATMENT, mean_tissue_biomass, CORAL_NUM, GENOTYPE) %>%
  summarize(pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()))

# create model for mean pH influence on mean coral tissue biomass
pH_biomass_model_rando <- lmer(mean_tissue_biomass ~ pH_mean + (1|GENOTYPE), data = chem_biomass_plotdata)
check_model(pH_biomass_model_rando)
summary(pH_biomass_model_rando)

pH_biomass_model <- lm(mean_tissue_biomass ~ pH_mean, data = chem_biomass_plotdata)
check_model(pH_biomass_model)
summary(pH_biomass_model) 

# remove influential plot 
chem_biomass_plotdata2 <- chem_biomass_plotdata %>% 
  filter(mean_tissue_biomass < 0.00075)

# new stats with outlier removed 
pH_biomass_model_rando2 <- lmer(mean_tissue_biomass ~ pH_mean + (1|GENOTYPE), data = chem_biomass_plotdata2)
check_model(pH_biomass_model_rando2)
summary(pH_biomass_model_rando2)

biomass_meanpH_plot <- chem_biomass_plotdata2 %>%
  ggplot(aes(x = pH_mean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) +
  labs(y = expression(bold("Mean Tissue Biomass" ~ (g ~ mL^-1 ~ cm^-2))), x= "Mean pH") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
biomass_meanpH_plot
ggsave(plot = biomass_meanpH_plot, filename = here("Output", "biomass_meanpH_plot.png"), width = 9, height = 7)


# mean tissue biomass and mean pH with random effect of tank number 
biomass_meanpH_tank <- lmer(mean_tissue_biomass ~ pH_mean + (1|TANK_NUM), data = chem_biomass_plotdata2)
check_model(biomass_meanpH_tank)
summary(biomass_meanpH_tank) # no sig effect of mean pH on mean tissue biomass

## TA and biomass data ## 
chem_biomass_plotdata <- chem_biomass_data_clean %>%
  group_by(TANK_NUM, TREATMENT, mean_tissue_biomass, CORAL_NUM, GENOTYPE) %>%
  summarize(pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()), 
            TA_mean = mean(TA, na.rm = TRUE), 
            TA_se = sd(TA, na.rm = TRUE)/sqrt(n()))
chem_biomass_plotdata

chem_biomass_plotdata <- chem_biomass_plotdata %>%
  filter(mean_tissue_biomass < 0.00075)

chem_biomass_plotdata$TREATMENT <- factor(chem_biomass_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# plot mean tissue biomass as a function of mean TA
biomass_TA_plot <- chem_biomass_plotdata %>%
  ggplot(aes(x = TA_mean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) +
  labs(y = expression(bold("Mean Tissue Biomass" ~ (g ~ mL^-1 ~ cm^-2))), x= expression(bold("Mean Total Alkalinity" ~ (µmol ~ kg^-1)))) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "#E31A1C", "Control" = "green4", "Coral_Dom" = "dodgerblue2",
                                "Rubble_Dom" = "#6A3D9A"))
biomass_TA_plot

#ggsave(plot = biomass_TA_plot, filename = here("Output", "biomass_TA_plot.png"), width = 9, height = 7)

# create model of mean tissue biomass as a function of mean TA 
TA_biomass_model <- lmer(mean_tissue_biomass ~ TA_mean + (1|TANK_NUM), data = chem_biomass_plotdata)
check_model(TA_biomass_model)
summary(TA_biomass_model)



# light vs diff pH plot
light_pH <- Data %>%
  ggplot(aes(x = Light_nm, y = pHDiff, color = Treatment))+
  geom_point()+
  geom_line() 
  #facet_wrap(~Date)
light_pH +
  scale_color_hue(labels = c("Algae-Dominated", "Control", "Rubble-Dominated", "Coral-Dominated"))
ggsave(plot = light_pH, filename = here("Output", "light_pH.png"), width = 15, height = 9)
