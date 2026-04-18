# =============================================================================
# Standardized Effect Sizes of pH and DOC on Coral Physiology — Revised
# Material Legacies / Indirect Effects Mesocosm Experiment
#
# Fits mixed-effects models for each physiological response variable as a
# function of mean tank pH and DOC (standardized). Reports standardized
# coefficients (effect sizes) with 95% CIs as a forest plot.
#
# KEY CHANGE FROM ORIGINAL: TANK_NUM is now included as a random effect
# alongside GENOTYPE to account for non-independence of coral fragments
# that share the same upstream community tank (Reviewer 2, Comment 2).
# =============================================================================

library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(performance)
library(parameters)
library(patchwork)

# ── 1. DATA LOADING ───────────────────────────────────────────────────────────

metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv")) |>
  mutate(
    TREATMENT = factor(TREATMENT,
                       levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom")),
    # Cap one extreme tissue biomass value identified during EDA
    mean_tissue_biomass = ifelse(mean_tissue_biomass > 30, NA, mean_tissue_biomass)
  )

# ── 2. STANDARDIZE PREDICTORS AND RESPONSES ───────────────────────────────────
# Log-transform skewed biological variables before z-scoring.
# pH and DOC are only z-scored (already on linear scales).

metadata_scaled <- metadata |>
  select(TREATMENT, GENOTYPE, TANK_NUM,
         chla_ug_cm2, endo_per_cm2, mean_tissue_biomass,
         GP, R, pH_mean, DOC_mean) |>
  mutate(
    across(c(chla_ug_cm2, endo_per_cm2, mean_tissue_biomass),
           ~ as.numeric(scale(log(.x)))),
    across(c(GP, R, pH_mean, DOC_mean),
           ~ as.numeric(scale(.x)))
  )

# ── 3. MODEL FITTING ──────────────────────────────────────────────────────────
# Response ~ pH + DOC + (1 | TANK_NUM) + (1 | GENOTYPE)
# Tank random effect accounts for shared water chemistry among fragments from
# the same community tank. Genotype random effect accounts for clone variation.

responses <- list(
  Endo_per_cm2        = "endo_per_cm2",
  Chla_ug_cm2         = "chla_ug_cm2",
  Mean_Tissue_Biomass = "mean_tissue_biomass",
  R                   = "R",
  GP                  = "GP"
)

fit_model <- function(response_var) {
  formula <- as.formula(
    paste0(response_var,
           " ~ pH_mean + DOC_mean + (1 | TANK_NUM) + (1 | GENOTYPE)")
  )
  lmer(formula, data = metadata_scaled, REML = TRUE)
}

models <- map(responses, fit_model)

# Print ANOVA tables for each model
walk2(names(models), models, function(nm, mod) {
  message("\n── ", nm, " ──")
  print(anova(mod))
})

# ── 4. EXTRACT STANDARDIZED COEFFICIENTS ─────────────────────────────────────

effect_table <- imap_dfr(models, function(mod, response_name) {
  tibble(model_parameters(mod)[1:3, ]) |>
    mutate(Phys = response_name)
}) |>
  filter(Parameter != "(Intercept)") |>
  mutate(
    Parameter = case_match(
      Parameter,
      "pH_mean"  ~ "Mean pH",
      "DOC_mean" ~ "Mean DOC"
    ),
    sig = ifelse(p < 0.055, 1, 0.5)
  )

# ── 5. FOREST PLOT ────────────────────────────────────────────────────────────

phys_labels <- c(
  Endo_per_cm2        = expression(atop("Endosymbiont Density",
                                        ~ (cells ~ 10^6 ~ cm^-2))),
  Chla_ug_cm2         = expression(atop("Chlorophyll-a Content",
                                        ~ (mu * g ~ cm^-2))),
  Mean_Tissue_Biomass = expression(atop("Tissue Biomass",
                                        ~ (mg ~ cm^-2))),
  R                   = expression(atop("Respiration Rate",
                                        ~ (mu * mol ~ O[2] ~ cm^-2 ~ hr^-1))),
  GP                  = expression(atop("Gross Photosynthesis",
                                        ~ (mu * mol ~ O[2] ~ cm^-2 ~ hr^-1)))
)

effect_plot <- ggplot(
  effect_table,
  aes(x = Coefficient, y = Phys, color = Parameter, alpha = sig)
) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 3, position = position_dodge(width = 0.4), shape = 19) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0, position = position_dodge(0.4)
  ) +
  scale_alpha(range = c(0.4, 1)) +
  guides(alpha = "none") +
  scale_color_manual(values = c("Mean pH" = "#082a54", "Mean DOC" = "#800000")) +
  scale_y_discrete(
    breaks = names(phys_labels),
    labels = phys_labels
  ) +
  labs(
    x     = "Standardized Effect Size",
    y     = "",
    color = ""
  ) +
  theme_bw() +
  theme(
    axis.title.x  = element_text(size = 14, color = "black"),
    axis.text.x   = element_text(size = 11, color = "black"),
    axis.text.y   = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.text   = element_text(size = 11)
  )

effect_plot

# ── 6. WRITE OUTPUT (uncomment to save) ───────────────────────────────────────
# ggsave(effect_plot, filename = here("Output", "Fig_4.png"), width = 6, height = 5)
# ggsave(effect_plot, filename = here("Revision", "Output", "Fig_4_revised.png"),
#        width = 6, height = 5)
