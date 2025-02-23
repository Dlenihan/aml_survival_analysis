# Set wd


install.packages("survminer")
install.packages("rstanarm")
install.packages("survminer")
install.packages("brms")
install.packages("survival")
install.packages("dplyr")
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

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