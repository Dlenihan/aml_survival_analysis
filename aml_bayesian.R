install.packages("survminer")
install.packages("rstanarm")
install.packages("survminer")
install.packages("brms")
install.packages("survival")
install.packages("dplyr")

search()
ls("package:rstanarm")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(rstanarm)
library(brms)
library(survival)
library(dplyr)

# Load the dataset
aml_data_outcomes <- read.csv("aml_followup.csv")

## Data cleaning
# Convert days_to_death and days_to_last_followup to numeric
aml_data_outcomes$days_to_death <- as.numeric(gsub("[^0-9]", "", aml_data_outcomes$days_to_death))
aml_data_outcomes$days_to_last_followup <- as.numeric(gsub("[^0-9]", "", aml_data_outcomes$days_to_last_followup))

# Compute survival time using last follow-up for censored patients
aml_data_outcomes$survival_time <- ifelse(
   aml_data_outcomes$vital_status == "DECEASED", 
   aml_data_outcomes$days_to_death, 
   aml_data_outcomes$days_to_last_followup
)

# Impute missing survival times using median survival time
median_survival <- median(aml_data_outcomes$survival_time, na.rm = TRUE)
aml_data_outcomes$survival_time[is.na(aml_data_outcomes$survival_time)] <- median_survival

## Troubleshooting
sum(is.na(aml_data_outcomes$survival_time))  # Should return 0 if correctly fixed

aml_data_outcomes %>%
   filter(is.na(survival_time)) %>%
   select(vital_status, days_to_last_followup, days_to_death)

# Encode categorical variables
aml_data_outcomes$sex <- ifelse(aml_data_outcomes$sex == "Male", 0, 1)
aml_data_outcomes$vital_status <- ifelse(aml_data_outcomes$vital_status == "LIVING", 0, 1)

# Convert categorical variables into factors
aml_data_outcomes <- aml_data_outcomes %>%
   mutate_at(vars(risk_cyto_at_diagnosis, risk_molecular_at_diagnosis, FAB), as.factor)

# Remove unnecessary columns
aml_data_outcomes_cleaned <- aml_data_outcomes %>%
   select(-c(patient_ID, race, days_to_birth, days_to_death, days_to_last_followup,
             date_of_form_completion, date_of_initial_pathologic_diagnosis,
             date_of_initial_pathologic_diagnosis.1))

## Exploratory Data Analysis
# Exploratory Data Analysis
# Histogram of survival times by vital status
ggplot(aml_data_outcomes_cleaned, aes(x = survival_time, fill = factor(vital_status))) +
   geom_histogram(binwidth = 500, alpha = 0.7, position = "identity") +
   scale_fill_manual(values = c("blue", "red"), labels = c("Alive", "Deceased")) +
   theme_minimal() +
   labs(title = "Distribution of Survival Times by Vital Status",
        x = "Survival Time (Days)", y = "Count", fill = "Vital Status")

# Kaplan-Meier Survival Curve
km_fit <- survfit(Surv(survival_time, vital_status) ~ FLT3_mut_status, data = aml_data_outcomes_cleaned)
ggsurvplot(km_fit, data = aml_data_outcomes_cleaned, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve by FLT3 Mutation Status")

# Bayesian Survival Modelling with brms
aml_bayesian_model <- brm(
   survival_time | cens(1 - vital_status) ~ FLT3_mut_status + age_at_diagnosis + WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis,
   data = aml_data_outcomes_cleaned,
   family = weibull(),  # Weibull survival model
   chains = 4, iter = 2000
)

# Summarize Bayesian Model
summary(aml_bayesian_model)
plot(aml_bayesian_model)


"""
## Explanatory Markdown Cell

### Overview
This script performs **Bayesian survival analysis** on **Acute Myeloid Leukemia (AML)** patient data in **RStudio**.

### Steps:

1. **Data Cleaning**:
   - Converts `days_to_death` to numeric.
   - Computes `survival_time`.
   - Encodes categorical variables (`sex`, `vital_status`).
   - Converts categorical variables into factors.
   - Drops unnecessary columns.

2. **Exploratory Data Analysis (EDA)**:
   - Plots **histogram of survival times**.
   - Generates a **Kaplan-Meier survival curve** grouped by `FLT3_mut_status`.

3. **Bayesian Survival Modelling**:
   - Uses `stan_surv()` from the **rstanarm** package to fit a **Bayesian survival model**.
   - Defines an **exponential baseline hazard**.
   - Uses Markov Chain Monte Carlo (**MCMC**) sampling.
   - Summarizes and plots posterior distributions.

This approach provides **a probabilistic understanding of AML survival probabilities** and allows us to compare different risk factors.
"""
