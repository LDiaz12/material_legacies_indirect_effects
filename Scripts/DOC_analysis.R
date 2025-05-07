library(tidyverse)
library(here)
library(moments)
library(emmeans)
library(agricolae)
library(ggrepel)
library(lubridate)

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





## DOC1 plot NPOC mg/l ##
# NPOC = non-purgeable organic carbon 
# indicator of organic matter 
DOC_NPOC_mgL <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = NPOC_mg_L, color = as.factor(TIME))) + 
  geom_point() +
  geom_text_repel(aes(label = TANK_NUM)) + 
  facet_wrap(~DATE)
DOC_NPOC_mgL
#ggsave(plot = DOC_NPOC_mgL, filename = here("Output", "DOCOutput", "DOC_NPOC_mgL.png"), width = 12, height = 9)

# create a "plot data" sheet to use for plotting means and se # 
DOC_plotdata <- DOC_full %>%
  group_by(TREATMENT) %>%
  summarise(NPOC_mg_L_mean = mean(NPOC_mg_L, na.rm = TRUE), 
            NPOC_mg_L_se = sd(NPOC_mg_L, na.rm = TRUE)/sqrt(n()),
            TN_mg_L_mean = mean(TN_mg_L, na.rm = TRUE), 
            TN_mg_L_se = sd(TN_mg_L, na.rm = TRUE)/sqrt(n()), 
            NPOC_uM_mean = mean(NPOC_uM, na.rm = TRUE), 
            NPOC_uM_se = sd(NPOC_uM, na.rm = TRUE)/sqrt(n()), 
            TN_uM_mean = mean(TN_uM, na.rm = TRUE), 
            TN_uM_se = sd(TN_uM, na.rm = TRUE)/sqrt(n()))
DOC_plotdata

# reorder factor levels for plotting to compare control to communities 
DOC_plotdata$TREATMENT <- factor(DOC_plotdata$TREATMENT, levels = c("Control", "Inflow", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

## DOC - mean NPOC mg/L by treatment plot ## 
DOC_mean_NPOC1 <- DOC_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NPOC_mg_L_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Mean NPOC" ~ (mg ~ L^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated",
                            "Inflow" = "Inflow")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = DOC_full, aes(x = TREATMENT, y = NPOC_mg_L), alpha = 0.7) +
  geom_errorbar(aes(ymin = NPOC_mg_L_mean - NPOC_mg_L_se,
                    ymax = NPOC_mg_L_mean + NPOC_mg_L_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan", "Inflow" = "orchid"))
DOC_mean_NPOC1
#ggsave(plot = DOC_mean_NPOC1, filename = here("Output", "DOCOutput", "DOC_mean_NPOC1.png"), width = 9, height = 7)
#####################
# TN_mgL and communities # 
DOC_TN_mgL <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = TN_mg_L, color = as.factor(TIME))) + 
  geom_point() +
  geom_text_repel(aes(label = TANK_NUM)) +
  facet_wrap(~DATE) 
DOC_TN_mgL

#ggsave(plot = DOC_TN_mgL, filename = here("Output", "DOCOutput", "DOC_TN_mgL.png"), width = 12, height = 10)
###########################
########################
# DOC - plot mean TN mg/L by treatment # 
DOC_mean_TN1 <- DOC_plotdata %>%
  ggplot(aes(x = TREATMENT, y = TN_mg_L_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Mean TN" ~ (mg ~ L^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated",
                            "Inflow" = "Inflow")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = DOC_full, aes(x = TREATMENT, y = TN_mg_L), alpha = 0.7) +
  geom_errorbar(aes(ymin = TN_mg_L_mean - TN_mg_L_se,
                    ymax = TN_mg_L_mean + TN_mg_L_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan", "Inflow" = "orchid"))
DOC_mean_TN1
#ggsave(plot = DOC_mean_TN1, filename = here("Output", "DOCOutput", "DOC_mean_TN1.png"), width = 9, height = 7)
###############################
# NPOC uM and communities # 

DOC_NPOC_uM <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = NPOC_uM, color = as.factor(TIME))) + 
  geom_point() +
  scale_y_continuous(trans = 'log10') +
  geom_text_repel(aes(label = TANK_NUM)) +
  facet_wrap(~DATE) 
DOC_NPOC_uM
ggsave(plot = DOC_NPOC_uM, filename = here("Output", "DOCOutput", "DOC_NPOC_uM.png"), width = 12, height = 10)

# DOC - plot mean NPOC uM by community # 

DOC_mean_NPOC2 <- DOC_plotdata %>%
  filter(TREATMENT != "Inflow") %>%
  ggplot(aes(x = TREATMENT, y = NPOC_uM_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Mean NPOC" ~ (uM)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated",
                            "Inflow" = "Inflow")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = DOC_full %>%
                filter(TIME %in% c("12:00:00","21:00:00"),
                       TREATMENT !="Inflow"), aes(x = TREATMENT, y = NPOC_uM), alpha = 0.7) +
  geom_errorbar(aes(ymin = NPOC_uM_mean - NPOC_uM_se,
                    ymax = NPOC_uM_mean + NPOC_uM_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan", "Inflow" = "orchid")) +
  facet_wrap(~TIME)
DOC_mean_NPOC2
#ggsave(plot = DOC_mean_NPOC2, filename = here("Output", "DOCOutput", "DOC_mean_NPOC2.png"), width = 9, height = 7)

################################
# TN uM and communities # 
DOC_TN_uM <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = TN_uM, color = as.factor(TIME))) + 
  geom_point() +
  scale_y_continuous(trans = 'log10') +
  geom_text_repel(aes(label = TANK_NUM)) +
  facet_wrap(~DATE) 
DOC_TN_uM
#ggsave(plot = DOC_TN_uM, filename = here("Output", "DOCOutput", "DOC_TN_uM.png"), width = 12, height = 10)

# DOC - plot mean TN uM by community # 

DOC_mean_TN2 <- DOC_plotdata %>%
  ggplot(aes(x = TREATMENT, y = TN_uM_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Mean TN" ~ (uM)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated",
                            "Inflow" = "Inflow")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = DOC_full, aes(x = TREATMENT, y = TN_uM), alpha = 0.7) +
  geom_errorbar(aes(ymin = TN_uM_mean - TN_uM_se,
                    ymax = TN_uM_mean + TN_uM_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan", "Inflow" = "orchid"))
DOC_mean_TN2
#ggsave(plot = DOC_mean_TN2, filename = here("Output", "DOCOutput", "DOC_mean_TN2.png"), width = 9, height = 7)
############################


## STATS ## 

## NPOC mg/L and treatment ## 
NPOCmgl_treatment_model <- lmer(NPOC_mg_L ~ TREATMENT + (1|TANK_NUM), data = DOC_full)
check_model(NPOCmgl_treatment_model)
summary(NPOCmgl_treatment_model) # coral dom communities slightly significant (p = 0.056)
anova(NPOCmgl_treatment_model)

## NPOC uM and treatment ## 
NPOCuM_treatment_model <- lmer(NPOC_uM ~ TREATMENT + (1|TANK_NUM), data = DOC_full)
check_model(NPOCuM_treatment_model)
summary(NPOCuM_treatment_model) 
anova(NPOCuM_treatment_model)
## same as mg/L stats 

## TN mg/L and treatment ## 
TNmgl_treatment_model <- lmer(TN_mg_L ~ TREATMENT + (1|TANK_NUM), data = DOC_full)
check_model(TNmgl_treatment_model)
summary(TNmgl_treatment_model) 
anova(TNmgl_treatment_model)
## results non significant 

## TN uM and treatment ## 
TNuM_treatment_model <- lmer(TN_uM ~ TREATMENT + (1|TANK_NUM), data = DOC_full)
check_model(TNuM_treatment_model)
summary(TNuM_treatment_model) 
anova(TNuM_treatment_model)
## results non significant 

