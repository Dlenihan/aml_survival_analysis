setwd("C:/Users/diarmuidlenihan/aml_work")
getwd()

# Install necessary packages if not already installed
if (!require("survival")) install.packages("survival", dependencies=TRUE)
if (!require("survminer")) install.packages("survminer", dependencies=TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies=TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies=TRUE)
if (!require("ggpubr")) install.packages("ggpubr", dependencies=TRUE)
if (!require("car")) install.packages("car", dependencies=TRUE)
install.packages("forestmodel")
install.packages("broom")

# Load required libraries
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(car)
library(forestmodel)
library(broom)

# Load the dataset (Ensure the path is correct)
aml_data <- read.csv("aml_dataset.csv", stringsAsFactors = TRUE)

# Convert to factors
aml_data$age_group <- as.factor(aml_data$age_group)
aml_data$FLT3_mut_status <- as.factor(aml_data$FLT3_mut_status)
aml_data$risk_cyto_at_diagnosis <- as.factor(aml_data$risk_cyto_at_diagnosis)
aml_data$risk_molecular_at_diagnosis <- as.factor(aml_data$risk_molecular_at_diagnosis)

# Convert survival time to years
aml_data$survival_time_years <- aml_data$survival_time / 365

# Ensure vital_status is a binary event variable (1 = deceased, 0 = alive)
aml_data$vital_status <- as.numeric(aml_data$vital_status == "DECEASED")

# Normality checks with histograms
normality_tests <- function(column, group, group_name) {
  cat("\nShapiro-Wilk normality test for", group_name, "\n")
  
  # Perform the Shapiro-Wilk test
  shapiro_result <- shapiro.test(group[[column]])
  print(shapiro_result)
  
  # Plot histogram
  ggplot(group, aes_string(x = column)) +  # aes_string() avoids tidy eval issues
    geom_histogram(binwidth = 50, fill = "blue", alpha = 0.5, color = "black") +
    labs(title = paste("Histogram of", column, "for", group_name),
         x = column, y = "Frequency") +
    theme_minimal()
}

# Subset data by FLT3 mutation status
wt_group <- subset(aml_data, FLT3_mut_status == 0)
mut_group <- subset(aml_data, FLT3_mut_status == 1)

# Check normality for survival time, BM blasts, and WBC count
normality_tests("survival_time", wt_group, "FLT3 WT")
normality_tests("survival_time", mut_group, "FLT3 Mut")
normality_tests("BM_blast_percentage_at_diagnosis", wt_group, "FLT3 WT")
normality_tests("BM_blast_percentage_at_diagnosis", mut_group, "FLT3 Mut")
normality_tests("WBC_count_at_diagnosis", wt_group, "FLT3 WT")
normality_tests("WBC_count_at_diagnosis", mut_group, "FLT3 Mut")

# Mann-Whitney U test (non-parametric comparisons)
mann_whitney_tests <- function(column) {
  test_result <- wilcox.test(
    aml_data[[column]] ~ aml_data$FLT3_mut_status, 
    data = aml_data, exact = FALSE
  )
  cat("\nMann-Whitney U test for", column, "\n")
  print(test_result)
}

mann_whitney_tests("survival_time")
mann_whitney_tests("BM_blast_percentage_at_diagnosis")
mann_whitney_tests("WBC_count_at_diagnosis")

# Boxplots for survival time, BM blasts, and WBC count
boxplot_graph <- function(column, y_label) {
  ggplot(aml_data, aes(x = as.factor(FLT3_mut_status), y = !!sym(column), fill = as.factor(FLT3_mut_status))) +
    geom_boxplot() +
    scale_fill_manual(values = c("blue", "red"), labels = c("WT", "Mut")) +
    labs(title = paste("Boxplot of", y_label, "by FLT3 Mutation Status"),
         x = "FLT3 Mutation Status", y = y_label) +
    theme_minimal()
}

boxplot_graph("survival_time", "Survival Time (days)")
boxplot_graph("BM_blast_percentage_at_diagnosis", "BM Blast % at Diagnosis")
boxplot_graph("WBC_count_at_diagnosis", "WBC Count at Diagnosis")

# Kaplan-Meier survival analysis
km_fit <- survfit(Surv(survival_time_years, vital_status) ~ FLT3_mut_status, data = aml_data)

# Plot Kaplan-Meier curves
ggsurvplot(km_fit, data = aml_data, 
           conf.int = TRUE, 
           pval = TRUE, 
           risk.table = TRUE,
           legend.title = "FLT3 Mutation Status",
           legend.labs = c("WT", "Mutant"),
           xlab = "Time (Years)", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal())

# Log-rank test
log_rank_test <- survdiff(Surv(survival_time_years, vital_status) ~ FLT3_mut_status, data = aml_data)
cat("\nLog-Rank Test Results:\n")
print(log_rank_test)

# Fit the Cox model using FLT3 mutation status as a predictor
cox_model <- coxph(Surv(survival_time_years, vital_status) ~ FLT3_mut_status, data = aml_data)

# Display model summary
summary(cox_model)

## Adjusting for additional risk factors

# Fit a Cox model adjusting for additional clinical factors
cox_model_adjusted <- coxph(Surv(survival_time_years, vital_status) ~ 
                              FLT3_mut_status + 
                              age_at_diagnosis + 
                              BM_blast_percentage_at_diagnosis + 
                              WBC_count_at_diagnosis + 
                              risk_cyto_at_diagnosis + 
                              risk_molecular_at_diagnosis, 
                            data = aml_data)

# Display adjusted model summary
summary(cox_model_adjusted)

## Forest plot of adjusted model
# Extract summary results
cox_results_adjusted <- tidy(cox_model_adjusted, exponentiate = TRUE, conf.int = TRUE)  # HR = exp(coef)

# Rename columns for clarity
cox_results_adjusted <- cox_results_adjusted %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::rename(
    Variable = term, 
    HR = estimate, 
    CI_Lower = conf.low, 
    CI_Upper = conf.high,
    P_Value = p.value
  )

## Removing risk scores - collinearity

# Modify model
cox_model_adjusted1 <- coxph(Surv(survival_time, vital_status) ~ 
                              FLT3_mut_status + age_at_diagnosis + 
                              WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis,
                            data = aml_data)

# Extract
cox_results_adjusted1 <- tidy(cox_model_adjusted1, exponentiate = TRUE, conf.int = TRUE)

# Rename columns for clarity
cox_results_adjusted1 <- cox_results_adjusted1 %>%
  dplyr::rename(
    Variable = term, 
    HR = estimate,  # Rename estimate → HR
    CI_Lower = conf.low,  # Rename conf.low → CI_Lower
    CI_Upper = conf.high  # Rename conf.high → CI_Upper
  )

# View results
print(cox_results_adjusted1)

# Forest plot
ggplot(cox_results_adjusted1, aes(y = Variable, x = HR, xmin = CI_Lower, xmax = CI_Upper)) +
  geom_point(color = "blue", size = 3) +  # HR point
  geom_errorbarh(height = 0.2, color = "black") +  # CI bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  theme_minimal() +
  labs(title = "Forest Plot of Adjusted Cox Model",
       x = "Hazard Ratio (HR)",
       y = "Covariates") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))

# Create the forest plot with log-transformed HRs
ggplot(cox_results_adjusted1, aes(y = Variable, x = HR, xmin = CI_Lower, xmax = CI_Upper)) +
  geom_point(color = "blue", size = 3) +  # HR points
  geom_errorbarh(height = 0.2, color = "black") +  # Confidence interval bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  scale_x_log10() +  # Log-transform the x-axis
  theme_minimal() +
  labs(title = "Forest Plot of Adjusted Cox Model",
       x = "Hazard Ratio (HR, log scale)",
       y = "Covariates") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))




## Retaining risk_cyto

# Modify model
cox_model_adjusted2 <- coxph(Surv(survival_time, vital_status) ~ 
                               FLT3_mut_status + age_at_diagnosis + 
                               WBC_count_at_diagnosis + BM_blast_percentage_at_diagnosis +risk_cyto_at_diagnosis,
                             data = aml_data)

# Extract
cox_results_adjusted2 <- tidy(cox_model_adjusted2, exponentiate = TRUE, conf.int = TRUE)

# Rename columns for clarity
cox_results_adjusted2 <- cox_results_adjusted2 %>%
  dplyr::rename(
    Variable = term, 
    HR = estimate,  # Rename estimate → HR
    CI_Lower = conf.low,  # Rename conf.low → CI_Lower
    CI_Upper = conf.high  # Rename conf.high → CI_Upper
  )

# View results
print(cox_results_adjusted2)

# Forest plot
ggplot(cox_results_adjusted2, aes(y = Variable, x = HR, xmin = CI_Lower, xmax = CI_Upper)) +
  geom_point(color = "blue", size = 3) +  # HR point
  geom_errorbarh(height = 0.2, color = "black") +  # CI bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  theme_minimal() +
  labs(title = "Forest Plot of Adjusted Cox Model",
       x = "Hazard Ratio (HR)",
       y = "Covariates") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))






# Create a clean dataset
hr_data <- data.frame(
  Variable = c("FLT3 Mutation", "Age at Diagnosis", "WBC Count", 
               "BM Blast %", "Cyto Risk: Intermediate", 
               "Cyto Risk: N.D.", "Cyto Risk: Poor"),
  HR = c(1.62, 1.04, 1.00, 0.999, 2.46, 9.08, 4.36),
  Lower = c(1.04, 1.03, 1.00, 0.989, 1.32, 2.83, 2.21),
  Upper = c(2.52, 1.06, 1.01, 1.01, 4.61, 29.2, 8.59),  # Large CI at 29.2
  p_value = c(0.0319, 0.0000000210, 0.0443, 0.771, 0.00474, 0.000212, 0.0000204)
)

# Order variables for better display
hr_data <- hr_data %>% arrange(HR)
hr_data$Variable <- factor(hr_data$Variable, levels = rev(hr_data$Variable))

# Cap extreme confidence intervals at 15 for display
hr_data <- hr_data %>%
  mutate(Upper_display = ifelse(Upper > 15, 15, Upper))  # Cap at 15

# Create the plot
ggplot(hr_data, aes(x = Variable, y = HR, ymin = Lower, ymax = Upper_display)) +
  geom_pointrange(color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  scale_y_log10(limits = c(0.5, 15)) +  # Adjust axis scale
  labs(title = "Forest Plot of Adjusted Cox Model",
       x = "Covariates",
       y = "Hazard Ratio (HR, log scale)") +
  theme_minimal(base_size = 12) +
  
  # Add annotation for capped CI
  annotate("text", x = "Cyto Risk: N.D.", y = 15, label = ">15", color = "red", size = 4)




# Check multicollinearity using VIF
vif(cox_model_adjusted)

## Checking proportionl hazards assumption

# Test the proportional hazards assumption
ph_test <- cox.zph(cox_model_adjusted)

# Print results
print(ph_test)

# Plot residuals to visually inspect proportionality
ggcoxzph(ph_test)

##Stratifying by risk group
# Filter data for high-risk patients only
aml_high_risk <- subset(aml_data, risk_cyto_at_diagnosis == "High")

# Fit a Cox model for high-risk patients
cox_model_high_risk <- coxph(Surv(survival_time_years, vital_status) ~ FLT3_mut_status, 
                             data = aml_high_risk)
## Removing risk_molecular
cox_model_fixed <- coxph(Surv(survival_time_years, vital_status) ~ 
                           FLT3_mut_status + 
                           age_at_diagnosis + 
                           BM_blast_percentage_at_diagnosis + 
                           WBC_count_at_diagnosis + 
                           risk_cyto_at_diagnosis, 
                         data = aml_data)

summary(cox_model_fixed)

# Display model summary
summary(cox_model_high_risk)

## Checking proportionl hazards assumption

ph_test <- cox.zph(cox_model_fixed)
print(ph_test)
ggcoxzph(ph_test)

## Stratifying by age

# Create age groups (you can adjust breakpoints as needed)
aml_data$age_group <- cut(
  aml_data$age_at_diagnosis, 
  breaks = c(0, 50, 65, Inf), 
  labels = c("≤50", "51-65", ">65")
)

cox_model_stratified <- coxph(Surv(survival_time_years, vital_status) ~ 
                                FLT3_mut_status + 
                                BM_blast_percentage_at_diagnosis + 
                                WBC_count_at_diagnosis + 
                                risk_cyto_at_diagnosis + 
                                strata(age_group), 
                              data = aml_data)

summary(cox_model_stratified)

ph_test_stratified <- cox.zph(cox_model_stratified)
print(ph_test_stratified)
ggcoxzph(ph_test_stratified)  # Visualize assumption check

## Viz

# Convert 'age_group' to a factor (if not already)
aml_data$age_group <- as.factor(aml_data$age_group)

# Create survival object
surv_obj <- survfit(Surv(survival_time_years, vital_status) ~ FLT3_mut_status + strata(age_group), 
                    data = aml_data)

ggsurvplot(surv_obj, 
           data = aml_data, 
           conf.int = TRUE, 
           pval = TRUE, 
           legend.title = "Group",  # More descriptive title
           legend.labs = c("WT ≤50", "Mutant ≤50", 
                           "WT 51-65", "Mutant 51-65", 
                           "WT >65", "Mutant >65"),  # 6 groups
           xlab = "Time (Years)", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal())


saveRDS(aml_data, "processed_aml_data.rds")  # Saves as an RDS file


## Time-dep adjustment for age

cox_model_tvc <- coxph(Surv(survival_time_years, vital_status) ~ 
                         FLT3_mut_status + 
                         tt(age_at_diagnosis) + 
                         BM_blast_percentage_at_diagnosis + 
                         WBC_count_at_diagnosis + 
                         risk_cyto_at_diagnosis, 
                       data = aml_data, tt = function(x, t, ...) x * log(t))

summary(cox_model_tvc)

## Viz adjusted survival curves

# Refit a Cox model without the time-dependent term to extract baseline hazard
cox_model_baseline <- coxph(Surv(survival_time_years, vital_status) ~ 
                              FLT3_mut_status + 
                              BM_blast_percentage_at_diagnosis + 
                              WBC_count_at_diagnosis + 
                              risk_cyto_at_diagnosis, 
                            data = aml_data)

# Extract baseline cumulative hazard
baseline_hazard <- basehaz(cox_model_baseline, centered = FALSE)

# Create new dataset for prediction (100 time points for each FLT3 status)
new_data <- expand.grid(
  FLT3_mut_status = c(0, 1),  # WT and Mutant
  BM_blast_percentage_at_diagnosis = mean(aml_data$BM_blast_percentage_at_diagnosis),
  WBC_count_at_diagnosis = mean(aml_data$WBC_count_at_diagnosis),
  risk_cyto_at_diagnosis = "Intermediate"  # Choose a common risk group
)

# Repeat dataset for multiple time points
time_points <- unique(baseline_hazard$time)  # Extract all time points
new_data <- new_data[rep(seq_len(nrow(new_data)), each = length(time_points)), ]
new_data$time <- rep(time_points, times = 2)  # Repeat time points for each FLT3 status

# Compute the linear predictor (risk score)
new_data$linear_predictor <- predict(cox_model_tvc, newdata = new_data, type = "lp")

# Merge baseline hazard into prediction dataset
new_data <- merge(new_data, baseline_hazard, by.x = "time", by.y = "time", all.x = TRUE)

# Compute survival probability at each time point
new_data$survival_prob <- exp(-new_data$hazard * exp(new_data$linear_predictor))

library(ggplot2)

ggplot(new_data, aes(x = time, y = survival_prob, color = as.factor(FLT3_mut_status))) +
  geom_step() +
  labs(title = "Estimated Survival Probabilities",
       x = "Time (Years)", 
       y = "Survival Probability",
       color = "FLT3 Mutation Status") +
  scale_color_manual(values = c("blue", "red"), labels = c("WT", "Mutant")) +
  theme_minimal()

# Plot adjusted survival curves based on Cox model
ggsurvplot(survfit(cox_model_tvc, data = aml_data), 
           data = aml_data, 
           conf.int = TRUE, 
           pval = TRUE, 
           legend.title = "FLT3 Mutation Status", 
           legend.labs = c("WT", "Mutant"),
           xlab = "Time (Years)", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal())




cox_model_interaction <- coxph(Surv(survival_time_years, vital_status) ~ 
                                 FLT3_mut_status * age_group + 
                                 BM_blast_percentage_at_diagnosis + 
                                 WBC_count_at_diagnosis + 
                                 risk_cyto_at_diagnosis, 
                               data = aml_data)

summary(cox_model_interaction)


# Extract coefficients
coef_main <- coef(cox_model_interaction)["FLT3_mut_status1"]
coef_age51_65 <- coef(cox_model_interaction)["FLT3_mut_status1"] + coef(cox_model_interaction)["FLT3_mut_status1:age_group51-65"]
coef_age65 <- coef(cox_model_interaction)["FLT3_mut_status1"] + coef(cox_model_interaction)["FLT3_mut_status1:age_group>65"]

# Convert to hazard ratios
HR_FLT3_age_50 <- exp(coef_main)
HR_FLT3_age_51_65 <- exp(coef_age51_65)
HR_FLT3_age_65 <- exp(coef_age65)

# Print results
HR_FLT3_age_50
HR_FLT3_age_51_65
HR_FLT3_age_65

# Compute variance for FLT3 mutation in each age group
var_main <- vcov(cox_model_interaction)["FLT3_mut_status1", "FLT3_mut_status1"]
var_age51_65 <- var_main + vcov(cox_model_interaction)["FLT3_mut_status1:age_group51-65", "FLT3_mut_status1:age_group51-65"] + 
  2 * vcov(cox_model_interaction)["FLT3_mut_status1", "FLT3_mut_status1:age_group51-65"]
var_age65 <- var_main + vcov(cox_model_interaction)["FLT3_mut_status1:age_group>65", "FLT3_mut_status1:age_group>65"] + 
  2 * vcov(cox_model_interaction)["FLT3_mut_status1", "FLT3_mut_status1:age_group>65"]

# Compute confidence intervals
CI_FLT3_age_50 <- exp(c(coef_main - 1.96 * sqrt(var_main), coef_main + 1.96 * sqrt(var_main)))
CI_FLT3_age_51_65 <- exp(c(coef_age51_65 - 1.96 * sqrt(var_age51_65), coef_age51_65 + 1.96 * sqrt(var_age51_65)))
CI_FLT3_age_65 <- exp(c(coef_age65 - 1.96 * sqrt(var_age65), coef_age65 + 1.96 * sqrt(var_age65)))

# Print results
CI_FLT3_age_50
CI_FLT3_age_51_65
CI_FLT3_age_65


# Data
age_groups <- c("≤50", "51-65", ">65")
hr_values <- c(1.721, 1.295, 2.712)  # Hazard Ratios
ci_lower <- c(0.806, 0.670, 1.294)  # Lower bounds of 95% CI
ci_upper <- c(3.675, 2.504, 5.686)  # Upper bounds of 95% CI

# Create a dataframe
forest_data <- data.frame(
  Age_Group = age_groups,
  HR = hr_values,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper
)

# Plot with labels
ggplot(forest_data, aes(x = HR, y = Age_Group)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "blue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = sprintf("%.2f", HR), x = HR + 0.3, y = as.numeric(Age_Group) +, size = 4, hjust = 0) +  # Adds HR labels
  labs(title = "Forest Plot of FLT3 Mutation Hazard Ratios by Age Group",
       x = "Hazard Ratio", y = "Age Group") +
  theme_minimal()

