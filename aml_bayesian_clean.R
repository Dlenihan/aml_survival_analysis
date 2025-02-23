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

table(aml_data_outcomes_cleaned$vital_status)

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

# Stratify by age group to correct proportional hazards issue
aml_data_outcomes_cleaned <- aml_data_outcomes_cleaned %>%
  mutate(age_group = cut(age_at_diagnosis, breaks=c(0, 50, 60, 70, Inf), labels=c("≤50", "51-60", "61-70", ">70")))
