
summary(aml_data_outcomes_cleaned$survival_time)

aml_data_outcomes_cleaned <- aml_data_outcomes_cleaned %>%
  mutate(survival_time = ifelse(survival_time <= 0, 1, survival_time))

# Bayesian Survival Modelling with brms
aml_bayesian_model <- brm(
  survival_time | cens(1 - vital_status) ~ FLT3_mut_status + age_at_diagnosis + WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis,
  data = aml_data_outcomes_cleaned,
  family = weibull(),
  chains = 4, iter = 2000
  #backend = "cmdstanr" # Speed up compilation
  save_model = "aml_model.stan"
)

## Diag
# Summarise Bayesian Model
summary(aml_bayesian_model)
plot(aml_bayesian_model) # Trace plots & density plots

posterior_summary(aml_bayesian_model)

pp_check(aml_bayesian_model)

pp_check(aml_bayesian_model, type = "hist")

# Lognormal
aml_bayesian_model_lognormal <- brm(
  survival_time | cens(1 - vital_status) ~ FLT3_mut_status + age_at_diagnosis + WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis,
  data = aml_data_outcomes_cleaned,
  family = lognormal(),
  chains = 4, iter = 2000
)

summary(aml_bayesian_model_lognormal)

pp_check(aml_bayesian_model_lognormal)

# Gamma
aml_bayesian_model_gamma <- brm(
  survival_time | cens(1 - vital_status) ~ FLT3_mut_status + age_at_diagnosis + WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis,
  data = aml_data_outcomes_cleaned,
  family = Gamma(link = "log"),
  chains = 4, iter = 2000
)

# Gamma
aml_bayesian_model_gamma <- brm(
  survival_time |  ~ FLT3_mut_status + age_at_diagnosis + WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis,
  data = aml_data_outcomes_cleaned,
  family = Gamma(link = "log"),
  chains = 4, iter = 2000
)

summary(aml_bayesian_model_gamma)

pp_check(aml_bayesian_model_gamma)

write.csv(aml_data_outcomes_cleaned, "aml_data_oucomes_cleaned", row.names = FALSE)

install.packages("writexl")  # Install if not already installed
library(writexl)

write_xlsx(aml_data_outcomes_cleaned, "aml_surv_time.xlsx")
