## Exploratory Data Analysis
# Exploratory Data Analysis
# Histogram of survival times by vital status
ggplot(aml_data_outcomes_cleaned, aes(x = survival_time, fill = factor(vital_status))) +
  geom_histogram(binwidth = 500, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("blue", "red"), labels = c("Alive", "Deceased")) +
  theme_minimal() +
  labs(title = "Distribution of Survival Times by Vital Status",
       x = "Survival Time (Days)", y = "Count", fill = "Vital Status")

# Histogram excluding censored patients
ggplot(aml_data_outcomes_cleaned %>% filter(vital_status == 1), 
       aes(x = survival_time)) +
  geom_histogram(binwidth = 250, fill = "red", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Survival Times (Deceased Patients Only)",
       x = "Survival Time (Days)", y = "Count")

# Kaplan-Meier Survival Curve
km_fit <- survfit(Surv(survival_time, vital_status) ~ FLT3_mut_status, data = aml_data_outcomes_cleaned)
ggsurvplot(km_fit, data = aml_data_outcomes_cleaned, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve by FLT3 Mutation Status")

survdiff(Surv(survival_time, vital_status) ~ FLT3_mut_status, data = aml_data_outcomes_cleaned)

# Stratified Cox Proportional Hazards Model to correct age proportionality issue
cox_model_stratified <- coxph(Surv(survival_time, vital_status) ~ FLT3_mut_status + WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis + strata(age_group),
                              data = aml_data_outcomes_cleaned)
summary(cox_model_stratified)

cox.zph(cox_model)
