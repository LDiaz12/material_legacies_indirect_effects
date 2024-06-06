library(tidyverse)
library(ggtext)
library(car)
jandata <- read_csv("~/Desktop/pH_Slope_Moorea_Jan_24.csv")
## Mean TA values by organism and dawn v dusk ##
mean_TAs <- jandata %>%
  group_by(Organism, Dusk_Dawn) %>%
  summarise(avg = mean(TA)) %>%
  arrange(avg)
mean_TAs

dawn_dusk_TA <- ggplot(data=mean_TAs, aes(fill=Organism, y=avg, x=Dusk_Dawn)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_color_fermenter() + 
  labs(x="Dusk and Dawn Samples", y="Total Alkalinity (mg/L CaCO3)", title="Diel Fluctuations in TA in Different Communities") +
  theme_bw()+
  theme(legend.position="bottom")
dawn_dusk_TA

## TA and mV modelling ##
par(mfrow=c(1,1))
TA_pH_model <- lm(TA~pH, data=jandata)
resid(TA_pH_model)
TA_pH_res <- resid(TA_pH_model)

qqp(TA_pH_res, "norm")
plot(TA_pH_res~fitted(TA_pH_model))
plot(TA_pH_model)
plot(TA~pH, data=jandata)
abline(TA_pH_model)
summary(TA_pH_model)

TA_pH_plot <- ggplot(jandata, aes(x=pH, y=TA))+
  geom_point(aes(colour = Organism, shape = Dusk_Dawn)) +
  geom_smooth(method = "lm", formula = y~x) + 
  labs(y = "Total Alkalinity (umol/kg)", x="pH") +
  theme_minimal()
TA_pH_plot

## Two-Way ANOVA TA Dawn Dusk Orgs ##
TA_TOD_orgs <- lm(TA~Dusk_Dawn*Organism, data=jandata)
anova(TA_TOD_orgs)
plot(TA_TOD_orgs)

#TA_TOD_orgs_plot <- ggplot(jandata, aes(x=Dusk_Dawn, y=TA))+
  #geom_point(aes(colour = Organism)) +
  #geom_smooth(method = "lm", formula = y~x, se=FALSE, fullrange=TRUE) + 
  #labs(y = "TA", x="Dusk or Dawn") +
  #theme_minimal()
#TA_TOD_orgs_plot

TA_TOD_orgs_plot <- ggplot(jandata, aes(y=TA, x=Dusk_Dawn, fill=Organism)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("Surface"="lightblue","Dead"="tan", "Pocillopora" = "coral", "Algae" = "lightgreen")) + 
  theme_bw() + 
  labs(y="Total Alkalinity (umol/kg)", title="Diel Fluctuations in Total Alkalinity in Different Communities") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="right") + 
  theme(text = element_text(size = 14),
        axis.text = element_text(size=12))
TA_TOD_orgs_plot

## Two-way ANOVA DO Dawn Dusk Orgs ##
DO_TOD_orgs <- lm(DO~Dusk_Dawn*Organism, data=jandata)
anova(DO_TOD_orgs)
plot(DO_TOD_orgs)


DO_TOD_orgs_plot <- ggplot(jandata, aes(y=DO, x=Dusk_Dawn, fill=Organism)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("Surface"="lightblue","Dead"="tan", "Pocillopora" = "coral", "Algae" = "lightgreen")) + 
  theme_bw() + 
  labs(y="% Dissolved Oxygen", title="Diel Fluctuations in Dissolved Oxygen in Different Communities") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="right") + 
  theme(text = element_text(size = 14), 
        axis.text = element_text(size = 12)) 
DO_TOD_orgs_plot

## Two-way ANOVA pH Dawn Dusk Orgs ##
pH_TOD_orgs <- lm(pH~Dusk_Dawn*Organism, data=jandata)
anova(pH_TOD_orgs)
plot(pH_TOD_orgs)

pH_TOD_orgs_plot <- ggplot(jandata, aes(y=pH, x=Dusk_Dawn, fill=Organism)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("Surface"="lightblue","Dead"="tan", "Pocillopora" = "coral", "Algae" = "lightgreen")) + 
  theme_bw() + 
  labs(y="pH", title="Diel Fluctuations in pH in Different Communities") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="right") +
  theme(text = element_text(size = 14), 
        axis.text = element_text(size = 12)) 
pH_TOD_orgs_plot

## Test Meso ##
test_meso <- read.csv('~/Desktop/TestMesoJan.csv')
## Mean DOs ##
test_mean_DOs <- test_meso %>%
  group_by(Tank, Dusk_or_Dawn) %>%
  summarise(avg = mean(DO)) %>%
  arrange(avg)
test_mean_DOs
## Plot DOs ##
attach(test_meso)
par(mfrow=c(3,1))
test_DO_mod <- lm(DO~Dusk_or_Dawn*Tank, data=test_meso)
anova(test_DO_mod)
plot(test_DO_mod)
library(emmeans)
graphdata3 <- as.data.frame(emmeans(test_DO_mod, ~ Dusk_or_Dawn*Tank))
graphdata3
test_DO_plot <- ggplot(graphdata3, aes(x=Dusk_or_Dawn, y=emmean, fill=factor(Tank))) + 
  geom_bar(stat="identity", position="dodge", linewidth=1.0) + 
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE), stat="identity", position=position_dodge(width=0.9), width=0.4) +
  labs(x="Dusk or Dawn", y="% Dissolved Oxygen", fill="Tank") + 
  scale_fill_manual(values=c("Control"="lightblue","Treatment"="darkgreen", "Response" = "coral")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
test_DO_plot
## Mean mVs ##
test_mean_mV <- test_meso %>%
  group_by(Tank, Dusk_or_Dawn) %>%
  summarise(avg = mean(mV)) %>%
  arrange(avg)
test_mean_mV
## Plot mV ##
test_mV_mod <- lm(mV~Dusk_or_Dawn*Tank, data=test_meso)
anova(test_mV_mod)
plot(test_mV_mod)


library(emmeans)
graphdata4 <- as.data.frame(emmeans(test_mV_mod, ~ Dusk_or_Dawn*Tank))
graphdata4
test_mV_plot <- ggplot(graphdata4, aes(x=Dusk_or_Dawn, y=emmean, fill=factor(Tank))) + 
  geom_bar(stat="identity", position="dodge", linewidth=1.0) + 
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE), stat="identity", position=position_dodge(width=0.9), width=0.4) +
  labs(x="Dusk or Dawn", y="mV", fill="Tank") + 
  scale_fill_manual(values=c("Control"="lightblue","Treatment"="darkgreen", "Response" = "coral")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
test_mV_plot
## Plot Salinities ##
test_sal_mod <- lm(Salinity~Dusk_or_Dawn*Tank, data=test_meso)
anova(test_sal_mod)
plot(test_sal_mod)
library(emmeans)
graphdata5 <- as.data.frame(emmeans(test_sal_mod, ~ Dusk_or_Dawn*Tank))
graphdata5
test_sal_plot <- ggplot(graphdata5, aes(x=Dusk_or_Dawn, y=emmean, fill=factor(Tank))) + 
  geom_bar(stat="identity", position="dodge", linewidth=1.0) + 
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE), stat="identity", position=position_dodge(width=0.9), width=0.4) +
  labs(x="Dusk or Dawn", y="Salinity (psu)", fill="Tank") + 
  scale_fill_manual(values=c("Control"="lightblue","Treatment"="darkgreen", "Response" = "coral")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
test_sal_plot

library(patchwork)
test_DO_plot + test_mV_plot + test_sal_plot
