library(tidyverse)

# ---- 1. ORIGINAL DATA (as you defined it) ----
bca2 <- tribble(
  ~conc_ug_ml, ~PS1,  ~PS2,  ~empty_plasm_r1, ~eGFP_r1, ~aav2_wt_r1, ~aav2_mut_r1, ~empty_plasm_r2, ~eGFP_r2, ~aav2_wt_r2, ~aav2_mut_r2,
  2000,        3.202, 3.352, 3.855,           3.87,     3.015,       3.148,        3.889,          3.608,    3.463,       3.729,
  1500,        2.489, 2.654, 3.789,           3.276,    3.256,       3.018,        3.876,          3.5,      3.307,       3.275,
  1000,        1.837, 1.747, 3.889,           3.289,    3.307,       2.978,        3.862,          3.376,    3.32,        3.274,
  750,         1.546, 1.506, 3.596,           3.479,    3.521,       2.978,        3.662,          3.369,    3.803,       3.717,
  500,         1.045, 1.012, 3.71,            3.471,    3.512,       2.979,        3.829,          3.369,    3.216,       3.231,
  250,         0.621, 0.621, 3.254,           3.532,    3.175,       2.878,        3.812,          3.157,    3.209,       3.082,
  125,         0.396, 0.402, 3.53,            3.266,    3.305,       2.946,        3.755,          3.232,    2.993,       3.048,
  0,           0.117, 0.113, 3.796,           3.338,    3.265,       3.035,        3.8,            3.347,    3.377,       2.73
)
# conc_ug_ml: known standard concentrations only for PS1/PS2. Row 8 is 0 µg/mL blank for PS1/PS2 only.

# ---- 2. DEFINE TRUE BLANK FROM PS1 & PS2 ROW 8 ----
blank_PS <- bca2 %>%
  filter(conc_ug_ml == 0) %>%
  summarise(blank = mean(c(PS1, PS2))) %>%
  pull(blank)

# ---- 3. BLANK-CORRECT PS1/PS2 AND FIT STANDARD CURVE ----
bca_stds <- bca2 %>%
  mutate(
    PS1_bc = PS1 - blank_PS,
    PS2_bc = PS2 - blank_PS,
    PS_mean_bc = rowMeans(across(c(PS1_bc, PS2_bc)))
  )

std_data <- bca_stds %>%
  # here, rows 1–7 are the real standards, row 8 is the 0 µg/mL blank
  select(conc_ug_ml, PS_mean_bc)

fit <- lm(PS_mean_bc ~ conc_ug_ml, data = std_data)
coef(fit)  # check: intercept ~ 0, slope > 0

# ---- 4. FUNCTION TO CONVERT ABSORBANCE TO CONCENTRATION (µg/mL) ----
predict_conc_ug_ml <- function(A_raw) {
  b <- coef(fit)[["(Intercept)"]]
  m <- coef(fit)[["conc_ug_ml"]]
  (A_raw - b) / m
}

# ---- 5. APPLY CURVE TO ALL SAMPLE WELLS (NO EXTRA BLANKS) ----
# We now use raw absorbance in each sample column.
sample_cols <- c("empty_plasm_r1", "eGFP_r1", "aav2_wt_r1", "aav2_mut_r1",
                 "empty_plasm_r2", "eGFP_r2", "aav2_wt_r2", "aav2_mut_r2")

long_results <- bca2 %>%
  select(conc_ug_ml, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample_id",
    values_to = "A_raw"
  ) %>%
  mutate(
    conc_ug_ml_calc = predict_conc_ug_ml(A_raw),
    conc_mg_ml_calc = conc_ug_ml_calc / 1000
  )

long_results
dilution_factor <- 3  # 25 µL sample + 50 µL buffer

long_results <- long_results %>%
  mutate(
    conc_mg_ml_plate    = conc_mg_ml_calc,
    conc_mg_ml_original = conc_mg_ml_calc * dilution_factor
  )


#Simplify the data
summary_table <- long_results %>%
  group_by(sample_id) %>%
  summarise(
    n_wells          = n(),
    mean_mg_ml       = mean(conc_mg_ml_original),
    sd_mg_ml         = sd(conc_mg_ml_original),
    se_mg_ml         = sd_mg_ml / sqrt(n_wells)
  )
summary_table
ggplot(long_results,
       aes(x = sample_id, y = conc_mg_ml_original)) +
  geom_boxplot() +
  ylab("Protein concentration (mg/mL, original sample)") +
  xlab("Sample")
