# Set wd

# Install libraries
install.packages("randomForest") # Creating random forests
install.packages("caret") # General ML
install.packages("ggplot2") # Data viz
install.packages("igraph") # Network analysis
install.packages("ppcor") # Partial correlation analysis
install.packages("infotheo") # Mutual information analysis
install.packages("GGally")  # For visualization
install.packages("ggraph")  # For improved network plotting
install.packages("factoextra")  # For clustering visualization

# Load packages
library(randomForest)
library(caret)
library(ggplot2)
library(igraph)
library(ppcor)
library(infotheo)
library(GGally)
library(ggraph)
library(factoextra)

# Load dataset
aml_data <- read.csv("Acute Myeloid Leukemia Mutation Dataset.csv")


# View
head(aml_data)

# Check structure
str(aml_data)

# Summ stats
summary(aml_data)

# Check for missing vals
colSums(is.na(aml_data))


# Clean data

# Rename columns
colnames(aml_data) <- c("Patient_ID", "FLT3_Mut_Status", "Sex", "Race", "FAB", "Age",
                        "BM_Blast_Percentage", "WBC", "Risk_Cyto", "Risk_Molecular")
# Remove unnecessary columns
aml_data <- aml_data[, -c(11:16)]

# Save cleaned data
write.csv(aml_data, "cleaned_aml_dataset", row.names = FALSE)


