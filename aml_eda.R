# EDA

# Histogram of Age Distribution
hist(aml_data$Age, main="Age Distribution of AML Patients", xlab="Age", col="blue")

# Boxplot for White Blood Cell Count (WBC)
boxplot(aml_data$WBC, main="WBC Count Distribution", ylab="WBC Count", col="red")

# Table of FLT3 Mutation Status
table(aml_data$FLT3_Mut_Status)

# Bar Plot of FLT3 Mutation Status
barplot(table(aml_data$FLT3_Mut_Status), main="FLT3 Mutation Status", 
        names.arg=c("Wild Type", "Mutant"), col=c("green", "orange"))

# Boxplot of BM Blast % by FLT3 Mutation Status
boxplot(BM_Blast_Percentage ~ FLT3_Mut_Status, data=aml_data, 
        main="Bone Marrow Blast % by FLT3 Mutation Status", 
        xlab="FLT3 Mutation Status", ylab="BM Blast %", col=c("purple", "yellow"))

# Compute correlation matrix
cor_matrix <- cor(aml_data[, c("Age", "BM_Blast_Percentage", "WBC")], use="complete.obs")

# Print the correlation matrix
print(cor_matrix)
