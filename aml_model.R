# Random forest model

# Cleaning


# New DF for conv to factors
aml_data_factors <- aml_data

# Convert categorical variables to factors
aml_data_factors$FLT3_Mut_Status <- as.factor(aml_data_factors$FLT3_Mut_Status)
aml_data_factors$Sex <- as.factor(aml_data_factors$Sex)
aml_data_factors$FAB <- as.factor(aml_data_factors$FAB)
aml_data_factors$Risk_Cyto <- as.factor(aml_data_factors$Risk_Cyto) # Target variable

# Count missing values in each column
colSums(is.na(aml_data_factors))

# Remove rows with missing values (Or use imputation?)
#aml_data <- na.omit(aml_data)

# Splitting

set.seed(123)  # Ensure reproducibility

# Create train-test split
train_index <- createDataPartition(aml_data_factors$Risk_Cyto, p = 0.8, list = FALSE)
train_data <- aml_data_factors[train_index, ]
test_data <- aml_data_factors[-train_index, ]

# Check split sizes
table(train_data$Risk_Cyto)
table(test_data$Risk_Cyto)

# Train rf model
rf_model <- randomForest(Risk_Cyto ~ Age + BM_Blast_Percentage + WBC + FLT3_Mut_Status + Sex + FAB, 
                         data = train_data, ntree = 500, importance = TRUE)

# View model summary
print(rf_model)
