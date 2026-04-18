# =============================================================================
# Carbonate Chemistry Analysis — Revised
# Material Legacies / Indirect Effects Mesocosm Experiment
# Mo'orea, French Polynesia | June 2024
#
# Calculates pH (from mV), in-situ pH adjustment, DIC, NEC, and NEP for
# each community tank at day (12:00) and night (21:00) sampling times.
# Produces Fig. 2 (NEC/NEP by treatment) and Fig. 3 (NCP–pH–DOC relationships).
# =============================================================================

library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(car)
library(lubridate)
library(lme4)
library(lmerTest)
library(performance)
library(ggpubr)
library(emmeans)
library(agricolae)
library(patchwork)
library(nlme)

# ── 1. DATA LOADING ───────────────────────────────────────────────────────────

pHcalib  <- read_csv(here("Data", "Chemistry", "TrisCalSummer2024.csv"))
pHData   <- read_csv(here("Data", "Chemistry", "CarbonateChemistry.csv"))
TableID  <- read_csv(here("Data", "TableID.csv"))
DOC_data <- read_csv(here("Data", "DOC", "DOC_full_data.csv"))

# ── 2. pH CALCULATION ─────────────────────────────────────────────────────────
# Fit Tris calibration (mV ~ temperature) per calibration date, then
# calculate sample pH at lab temp, then adjust to in-situ temperature.

pHSlope <- pHcalib |>
  nest_by(TrisCalDate) |>
  mutate(fitpH = list(lm(mVTris ~ TTris, data = pHcalib))) |>
  reframe(broom::tidy(fitpH)) |>
  select(TrisCalDate, term, estimate) |>
  pivot_wider(names_from = term, values_from = estimate) |>
  right_join(pHData, by = "TrisCalDate") |>
  mutate(mVTris = TEMPINLAB * TTris + `(Intercept)`) |>
  drop_na(TEMPINSITU, mV) |>
  mutate(pH = pH(Ex = mV, Etris = mVTris, S = SALINITY, T = TEMPINLAB)) |>
  drop_na(TEMPINSITU, TEMPINLAB, SALINITY, pH) |>
  mutate(pH_insitu = pHinsi(
    pH = pH, ALK = 2200,
    Tinsi = TEMPINSITU, Tlab = TEMPINLAB,
    S = SALINITY, Pt = 0.1,
    k1k2 = "m10", kf = "dg"
  )) |>
  select(-pH) |>
  rename(pH = pH_insitu) |>
  ungroup() |>
  select(-c(mV, TrisCalDate, TTris, `(Intercept)`, mVTris))

# ── 3. DIC CALCULATION ────────────────────────────────────────────────────────
# flag = 8: use pH + TA to derive remaining carbonate system variables.

DIC_calc <- pHSlope |>
  drop_na(pH, TA) |>
  mutate(TA_mol_kg = TA / 1e6)

carb_table <- carb(
  flag   = 8,
  var1   = DIC_calc$pH,
  var2   = DIC_calc$TA_mol_kg,
  S      = DIC_calc$SALINITY,
  T      = DIC_calc$TEMPINSITU,
  P      = 0, Patm = 0, Pt = 0, Sit = 0,
  pHscale = "T", kf = "dg", k1k2 = "m10", ks = "d"
)

# Append DIC (µmol kg⁻¹) back to working data frame
DIC_calc <- DIC_calc |>
  bind_cols(carb_table |> select(DIC) |> mutate(DIC_umol_kg = DIC * 1e6)) |>
  select(-TA_mol_kg, -DIC)

pHSlope2 <- pHSlope |>
  left_join(DIC_calc |> select(any_of(names(pHSlope)), DIC_umol_kg),
            by = intersect(names(pHSlope), names(DIC_calc))) |>
  mutate(TREATMENT = ifelse(is.na(TREATMENT), "Inflow", TREATMENT)) |>
  full_join(DOC_data |> select(-DATETIME), by = intersect(names(pHSlope), names(DOC_data))) |>
  select(-c(NPOC_mg_L, TN_mg_L))

# ── 4. INFLOW PAIRING ─────────────────────────────────────────────────────────
# Separate inflow measurements and join them to their paired community tanks.

InflowData <- pHSlope2 |>
  filter(TANK_NUM %in% c("Inflow1", "Inflow2")) |>
  select(-c(FLOW_LEFT, FLOW_RIGHT, Notes, DO_MG_L, SALINITY, TEMPINSITU)) |>
  rename(
    pH_inflow  = pH,
    TA_inflow  = TA,
    DIC_inflow = DIC_umol_kg,
    DOC_inflow = NPOC_uM
  ) |>
  mutate(INFLOW_TABLE = ifelse(TANK_NUM == "Inflow1", 1, 2)) |>
  ungroup() |>
  select(DATE, TIME, INFLOW_TABLE, pH_inflow, TA_inflow, DIC_inflow, DOC_inflow)

# ── 5. NEC AND NEP CALCULATIONS ───────────────────────────────────────────────
# Tank surface area (planar): 22.5 × 22.5 cm
SurfaceArea <- 22.5 * 22.5  # cm²

Data <- pHSlope2 |>
  ungroup() |>
  filter(!TANK_NUM %in% c("Inflow1", "Inflow2")) |>
  mutate(TANK_NUM = as.numeric(TANK_NUM)) |>
  left_join(TableID |> select(TANK_NUM, INFLOW_TABLE), by = "TANK_NUM") |>
  left_join(InflowData, by = c("DATE", "TIME", "INFLOW_TABLE")) |>
  mutate(
    DATETIME       = ymd_hms(paste(DATE, TIME)),
    deltapH        = pH - pH_inflow,
    deltaDOC       = NPOC_uM - DOC_inflow,
    totalflow      = FLOW_RIGHT + FLOW_LEFT,
    # Residence time (h): volume (10 L = 10,000 mL) / flow (mL min⁻¹) / 60
    residence_time = (1 / totalflow) * (10000 / 60),
    flowrate       = totalflow / 60,
    deltaTA        = TA_inflow - TA,
    deltaDIC       = DIC_inflow - DIC_umol_kg,
    # NEC (mmol CaCO₃ m⁻² h⁻¹): alkalinity anomaly method
    NEC = (deltaTA / 2) * 1.025 * 10 * (1 / residence_time) * (1 / SurfaceArea),
    # NEP (mmol C m⁻² h⁻¹): DIC anomaly method, corrected for CaCO₃ contribution
    NEP = (deltaDIC * 1.025 * 10 * (1 / residence_time) * (1 / SurfaceArea)) - NEC
  )

# ── 6. TANK FLOW / RESIDENCE TIME SUMMARY ────────────────────────────────────

tank_residence_time <- Data |>
  group_by(TANK_NUM) |>
  summarize(
    avg_res_time   = mean(residence_time, na.rm = TRUE),
    res_error      = sd(residence_time, na.rm = TRUE) / sqrt(n()),
    avg_flow_rate  = mean(flowrate, na.rm = TRUE),
    flow_error     = sd(flowrate, na.rm = TRUE) / sqrt(n())
  ) |>
  drop_na()

avg_total_flow <- Data |>
  summarize(
    avg_total_flow = mean(flowrate, na.rm = TRUE),
    flow_error     = sd(flowrate, na.rm = TRUE) / sqrt(n())
  )

# ── 7. DAILY SUMMARIES (DAY and NIGHT) ────────────────────────────────────────
# Daytime = 12:00; Nighttime = 21:00

chem_reframe_DAY <- Data |>
  filter(TIME == hms::as_hms("12:00:00")) |>
  group_by(TREATMENT, DATE, TANK_NUM) |>
  select(DATE, TIME, TANK_NUM, TREATMENT,
         TA, pH, DIC_umol_kg, NPOC_uM,
         deltapH, deltaTA, deltaDOC, deltaDIC, NEC, NEP) |>
  rename(DOC = NPOC_uM)

chem_reframe_NIGHT <- Data |>
  filter(TIME == hms::as_hms("21:00:00")) |>
  group_by(TREATMENT, DATE, TANK_NUM) |>
  select(DATE, TIME, TANK_NUM, TREATMENT,
         TA, pH, DIC_umol_kg, NPOC_uM,
         deltapH, deltaTA, deltaDOC, deltaDIC, NEC, NEP) |>
  rename(DOC = NPOC_uM)

# ── 8. OUTLIER REMOVAL ────────────────────────────────────────────────────────
# Thresholds set by visual inspection of raw tank × date distributions.
# Outliers represent measurement or flow-rate artefacts, not biological signal.

chem_reframe_DAY_clean <- chem_reframe_DAY |>
  mutate(
    deltaDIC     = ifelse(deltaDIC < -250,  NA, deltaDIC),
    deltaDOC     = ifelse(deltaDOC > 200,   NA, deltaDOC),
    deltaTA      = ifelse(deltaTA < -400,   NA, deltaTA),
    DIC_umol_kg  = ifelse(DIC_umol_kg > 2250, NA, DIC_umol_kg),
    DOC          = ifelse(DOC > 350,         NA, DOC),
    NEC          = ifelse(NEC < -2,          NA, NEC),
    NEC          = ifelse(TREATMENT == "Rubble_Dom" & NEC > 2, NA, NEC),
    NEP          = ifelse(NEP < -2.1,        NA, NEP),
    NEP          = ifelse(NEP > 5,           NA, NEP),
    pH           = ifelse(pH < 7.8,          NA, pH),
    TA           = ifelse(TA > 2700,         NA, TA)
  )

chem_reframe_NIGHT_clean <- chem_reframe_NIGHT |>
  mutate(
    deltaDOC    = ifelse(deltaDOC > 300,  NA, deltaDOC),
    deltaTA     = ifelse(deltaTA > 100,   NA, deltaTA),
    DIC_umol_kg = ifelse(DIC_umol_kg < 1900, NA, DIC_umol_kg),
    DOC         = ifelse(DOC > 400,       NA, DOC),
    NEC         = ifelse(NEC > 1,         NA, NEC),
    NEC         = ifelse(NEC < -1,        NA, NEC),
    NEP         = ifelse(NEP > 2,         NA, NEP),
    NEP         = ifelse(TREATMENT == "Control" & NEP > 1, NA, NEP),
    TA          = ifelse(TA < 2250,       NA, TA)
  )

# ── 9. TANK-LEVEL SUMMARY STATISTICS ─────────────────────────────────────────

summarize_chem <- function(df) {
  df |>
    group_by(TREATMENT, TANK_NUM) |>
    summarize(
      TA_mean        = mean(TA,      na.rm = TRUE),
      TA_se          = sd(TA,        na.rm = TRUE) / sqrt(n()),
      deltaTA_mean   = mean(deltaTA, na.rm = TRUE),
      deltaTA_se     = sd(deltaTA,   na.rm = TRUE) / sqrt(n()),
      pH_mean        = mean(pH,      na.rm = TRUE),
      pH_se          = sd(pH,        na.rm = TRUE) / sqrt(n()),
      deltapH_mean   = mean(deltapH, na.rm = TRUE),
      deltapH_se     = sd(deltapH,   na.rm = TRUE) / sqrt(n()),
      DOC_mean       = mean(DOC,     na.rm = TRUE),
      DOC_se         = sd(DOC,       na.rm = TRUE) / sqrt(n()),
      deltaDOC_mean  = mean(deltaDOC,na.rm = TRUE),
      deltaDOC_se    = sd(deltaDOC,  na.rm = TRUE) / sqrt(n()),
      NEC_mean       = mean(NEC,     na.rm = TRUE),
      NEC_se         = sd(NEC,       na.rm = TRUE) / sqrt(n()),
      NEP_mean       = mean(NEP,     na.rm = TRUE),
      NEP_se         = sd(NEP,       na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

chem_summary_data_DAY   <- summarize_chem(chem_reframe_DAY_clean)
chem_summary_data_NIGHT <- summarize_chem(chem_reframe_NIGHT_clean)

# ── 10. NEC ANALYSIS ──────────────────────────────────────────────────────────

treatment_levels <- c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom")
treatment_colors <- c(
  "Control"    = "blue",
  "Algae_Dom"  = "darkgreen",
  "Coral_Dom"  = "coral",
  "Rubble_Dom" = "tan"
)
treatment_labels <- c(
  "Algae_Dom"  = "Macroalgae-Enriched",
  "Control"    = "Control",
  "Coral_Dom"  = "Coral-Enriched",
  "Rubble_Dom" = "CCA-Enriched"
)

# Daytime NEC
NEC_data_day <- chem_reframe_DAY_clean |>
  group_by(TREATMENT, DATE, TANK_NUM) |>
  reframe(NEC_day_mean = mean(NEC, na.rm = TRUE)) |>
  drop_na() |>
  mutate(TREATMENT = factor(TREATMENT, levels = treatment_levels))

NEC_day_mean_plot <- NEC_data_day |>
  ggplot(aes(x = TREATMENT, y = NEC_day_mean, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar",
               fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  labs(x = "",
       y = expression(bold("NCC" ~ (mmol ~ CaCO[3] ~ m^-2 ~ h^-1)))) +
  scale_x_discrete(labels = treatment_labels) +
  scale_color_manual(values = treatment_colors) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

NEC_daytime_model <- lmer(NEC_day_mean ~ TREATMENT + (1 | TANK_NUM), data = NEC_data_day)
check_model(NEC_daytime_model)
summary(NEC_daytime_model)
anova(NEC_daytime_model)
emmeans(NEC_daytime_model, pairwise ~ TREATMENT, adjust = "Tukey")

NEC_daytime_model_lm <- lm(NEC_day_mean ~ TREATMENT, data = NEC_data_day)
HSD.test(NEC_daytime_model_lm, "TREATMENT", console = TRUE)

# Nighttime NEC
NEC_data_night <- chem_reframe_NIGHT_clean |>
  group_by(TREATMENT, DATE, TANK_NUM) |>
  reframe(NEC_night_mean = mean(NEC, na.rm = TRUE)) |>
  drop_na() |>
  mutate(TREATMENT = factor(TREATMENT, levels = treatment_levels))

NEC_night_mean_plot <- NEC_data_night |>
  ggplot(aes(x = TREATMENT, y = NEC_night_mean, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar",
               fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  labs(x = "", y = "") +
  scale_x_discrete(labels = treatment_labels) +
  scale_color_manual(values = treatment_colors) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

NEC_nighttime_model <- lmer(NEC_night_mean ~ TREATMENT + (1 | TANK_NUM), data = NEC_data_night)
check_model(NEC_nighttime_model)
summary(NEC_nighttime_model)
anova(NEC_nighttime_model)
emmeans(NEC_nighttime_model, pairwise ~ TREATMENT, adjust = "Tukey")

# ── 11. NEP ANALYSIS ──────────────────────────────────────────────────────────

# Daytime NEP
NEP_data_day <- chem_reframe_DAY_clean |>
  group_by(TREATMENT, DATE, TANK_NUM) |>
  reframe(NEP_day_mean = mean(NEP, na.rm = TRUE)) |>
  drop_na() |>
  mutate(TREATMENT = factor(TREATMENT, levels = treatment_levels))

NEP_day_mean_plot <- NEP_data_day |>
  ggplot(aes(x = TREATMENT, y = NEP_day_mean, color = TREATMENT)) +
  ggtitle("Daytime Mean") +
  geom_hline(yintercept = 0, lty = 2) +
  ylim(c(-3, 5)) +
  geom_point(size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar",
               fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  labs(x = "",
       y = expression(bold("NCP" ~ (mmol ~ C ~ m^-2 ~ h^-1)))) +
  scale_x_discrete(labels = treatment_labels) +
  scale_color_manual(values = treatment_colors) +
  theme_bw() +
  theme(
    axis.text.x     = element_blank(),
    axis.text.y     = element_text(size = 14),
    axis.title      = element_text(size = 14, face = "bold"),
    plot.title      = element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  )

NEP_day_model <- lmer(NEP_day_mean ~ TREATMENT + (1 | TANK_NUM), data = NEP_data_day)
check_model(NEP_day_model)
summary(NEP_day_model)
anova(NEP_day_model)
emmeans(NEP_day_model, pairwise ~ TREATMENT, adjust = "Tukey")

# Nighttime NEP
NEP_data_night <- chem_reframe_NIGHT_clean |>
  group_by(TREATMENT, DATE, TANK_NUM) |>
  reframe(NEP_night_mean = mean(NEP, na.rm = TRUE)) |>
  drop_na() |>
  mutate(TREATMENT = factor(TREATMENT, levels = treatment_levels))

NEP_night_mean_plot <- NEP_data_night |>
  ggplot(aes(x = TREATMENT, y = NEP_night_mean, color = TREATMENT)) +
  ggtitle("Nighttime Mean") +
  geom_hline(yintercept = 0, lty = 2) +
  ylim(c(-2, 2)) +
  geom_point(size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar",
               fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  labs(x = "", y = "") +
  scale_x_discrete(labels = treatment_labels) +
  scale_color_manual(values = treatment_colors) +
  theme_bw() +
  theme(
    axis.text.x     = element_blank(),
    axis.text.y     = element_text(size = 12),
    axis.title      = element_text(size = 14, face = "bold"),
    plot.title      = element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  )

NEP_night_model <- lmer(NEP_night_mean ~ TREATMENT + (1 | TANK_NUM), data = NEP_data_night)
check_model(NEP_night_model)
summary(NEP_night_model)
anova(NEP_night_model)

# ── 12. FIGURE 2: NEP AND NEC PATCHWORK ───────────────────────────────────────

Fig_2 <- (NEP_day_mean_plot + NEP_night_mean_plot) /
         (NEC_day_mean_plot + NEC_night_mean_plot) +
  plot_annotation(tag_levels = "a")
Fig_2

# ggsave(Fig_2, filename = here("Output", "Fig_2.png"), width = 12, height = 10)

# ── 13. NCP–pH–DOC RELATIONSHIPS (FIG. 3) ────────────────────────────────────
# Combine day and night cleaned data; filter extreme NEP values.

Clean_Chem_all <- chem_reframe_NIGHT_clean |>
  bind_rows(chem_reframe_DAY_clean) |>
  filter(NEP < 4) |>
  mutate(TREATMENT = factor(TREATMENT, levels = treatment_levels))

# Models (ANCOVA: response ~ NEP * treatment, with DATE as block)
mod_DOC_NEP <- lm(DOC ~ NEP * TREATMENT + (1 | DATE), data = Clean_Chem_all)
anova(mod_DOC_NEP)
summary(mod_DOC_NEP)

mod_NEC_pH  <- lm(NEC ~ pH * TREATMENT + (1 | DATE),  data = Clean_Chem_all)
anova(mod_NEC_pH)
summary(mod_NEC_pH)

# Scatter panels — color by treatment, single pooled regression line
plot_NEP_pH <- Clean_Chem_all |>
  ggplot(aes(x = NEP, y = pH)) +
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = expression("NCP (mmol C m"^-2 ~ "hr"^-1 ~ ")"),
    y = expression("pH"[T])
  ) +
  scale_color_manual(
    labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
    values = c("blue", "darkgreen", "coral", "tan")
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

plot_NEP_DOC <- Clean_Chem_all |>
  ggplot(aes(x = NEP, y = DOC)) +
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = expression("NCP (mmol C m"^-2 ~ "hr"^-1 ~ ")"),
    y = expression("DOC (" ~ mu ~ "mol L"^-1 ~ ")")
  ) +
  scale_color_manual(
    labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
    values = c("blue", "darkgreen", "coral", "tan")
  ) +
  theme_bw() +
  theme(
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.text    = element_text(size = 14)
  )

plot_pH_NEC <- Clean_Chem_all |>
  ggplot(aes(x = pH, y = NEC)) +
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = expression("pH"[T]),
    y = expression("NCC" ~ (mmol ~ CaCO[3] ~ m^-2 ~ hr^-1))
  ) +
  scale_color_manual(
    labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
    values = c("blue", "darkgreen", "coral", "tan")
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

plot_NEC_DOC <- Clean_Chem_all |>
  ggplot(aes(x = NEC, y = DOC)) +
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5) +
  labs(
    x = expression("NCC" ~ (mmol ~ CaCO[3] ~ m^-2 ~ hr^-1)),
    y = expression("DOC (" ~ mu ~ "mol L"^-1 ~ ")")
  ) +
  scale_color_manual(
    labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
    values = c("blue", "darkgreen", "coral", "tan")
  ) +
  theme_bw() +
  theme(
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.text    = element_text(size = 14)
  )

Fig_3 <- (plot_NEP_pH + plot_pH_NEC) /
          (plot_NEP_DOC + plot_NEC_DOC) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
Fig_3

# ggsave(Fig_3, filename = here("Output", "Fig_3.png"), width = 12, height = 12)

# ── WRITE OUTPUTS (uncomment to save) ─────────────────────────────────────────
# write_csv(chem_summary_data_DAY,   here("Data", "Chemistry", "chem_summary_data_DAY.csv"))
# write_csv(chem_summary_data_NIGHT, here("Data", "Chemistry", "chem_summary_data_NIGHT.csv"))
# write_csv(tank_residence_time,     here("Data", "Chemistry", "Tank_Residence_Time.csv"))
# write_csv(Data,                    here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))
# ggsave(Fig_2, filename = here("Output", "Fig_2.png"), width = 12, height = 10)
# ggsave(Fig_3, filename = here("Output", "Fig_3.png"), width = 12, height = 12)
