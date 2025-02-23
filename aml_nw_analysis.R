# Prep data

# Select key clinical variables (excluding risk scores)
network_data <- aml_data[, c("Age", "BM_Blast_Percentage", "WBC", "Sex", "FAB", "FLT3_Mut_Status")]

# Conv categorical variables to numeric
network_data$Sex <- as.numeric(as.factor(network_data$Sex))
network_data$FAB <- as.numeric(as.factor(network_data$FAB))
network_data$FLT3_Mut_Status <- as.numeric(as.factor(network_data$FLT3_Mut_Status))

# Conv continuous variables to discrete categories (Binning)
network_data$Age <- as.integer(discretize(network_data$Age, disc = "equalfreq", nbins = 4)$X)
network_data$BM_Blast_Percentage <- as.integer(discretize(network_data$BM_Blast_Percentage, disc = "equalfreq", nbins = 4)$X)
network_data$WBC <- as.integer(discretize(network_data$WBC, disc = "equalfreq", nbins = 4)$X)

# Convert back to df
network_data <- as.data.frame(network_data)

# Compute mutual information matrix on discretised data
mi_matrix <- mutinformation(network_data)

# Compute pearson correlation matrix on original data
cor_matrix <- cor(network_data, use = "complete.obs")

# Combine MI and Pearson Correlation*
# - Keeps edges if MI > 0.05 OR correlation > 0.2
adj_matrix <- (mi_matrix > 0.05) | (abs(cor_matrix) > 0.2)

# Convert adjacency matrix to an igraph object
graph_obj <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# Remove self-loops (neat)
graph_obj <- simplify(graph_obj, remove.loops = TRUE)

# Assign node labels (feature names)
V(graph_obj)$name <- colnames(network_data)

# Set node colors
V(graph_obj)$color <- "orange"

# Set edge weights (thicker lines = stronger relationships)
E(graph_obj)$width <- E(graph_obj)$weight * 5

# Compute feature importance (Total MI Contribution)
feature_importance <- colSums(mi_matrix)
print(sort(feature_importance, decreasing = TRUE))

# Better layout for improved structure
plot(graph_obj, 
     layout = layout_with_kk(graph_obj),  # Kamada-Kawai layout
     vertex.size = 15, 
     vertex.label.color = "black", 
     vertex.label.cex = 1.2,
     edge.color = "gray50",
     main = "Final AML Feature Network (Hybrid MI + Correlation)")
