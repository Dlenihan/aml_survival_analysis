# Load AML dataset (assuming it's already in a data frame)
patient_data <- aml_data[, c("Age", "BM_Blast_Percentage", "WBC", "Sex", "FAB", "FLT3_Mut_Status")]

# Convert categorical variables to numeric
patient_data$Sex <- as.numeric(as.factor(patient_data$Sex))
patient_data$FAB <- as.numeric(as.factor(patient_data$FAB))
patient_data$FLT3_Mut_Status <- as.numeric(as.factor(patient_data$FLT3_Mut_Status))

# **Compute Similarity Between Patients**
# Standardize numeric features to give equal importance
patient_data_scaled <- scale(patient_data)

# Compute Euclidean distance between patients
patient_distance <- dist(patient_data_scaled, method = "euclidean")

# Convert distance to similarity (1 - normalized distance)
similarity_matrix <- as.matrix(1 / (1 + patient_distance))

# Set a similarity threshold (retain only meaningful connections)
threshold <- 0.5  # Adjust this based on network density
adj_matrix <- similarity_matrix > threshold

# Convert to an igraph object
patient_graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# Assign patient IDs as node labels
V(patient_graph)$name <- rownames(aml_data)

# Set node color based on FLT3 mutation status
V(patient_graph)$color <- ifelse(patient_data$FLT3_Mut_Status == 1, "red", "blue")

# Set edge weights based on similarity score
E(patient_graph)$width <- E(patient_graph)$weight * 5

# **Detect Patient Clusters**
patient_clusters <- cluster_walktrap(patient_graph)

# **Plot the Patient Similarity Network (Fixed Syntax)**
ggraph(patient_graph, layout = "fr") + 
  geom_edge_link(aes(width = weight), alpha = 0.5) +
  geom_node_point(aes(color = factor(V(patient_graph)$color)), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("blue", "red"), name = "FLT3 Status") +  # Fixed color scale
  ggtitle("AML Patient Similarity Network") +
  theme_void()

# Identify Patient Clusters with Community Detection

# Detect patient communities using Louvain clustering
patient_clusters <- cluster_louvain(patient_graph)

# Assign cluster labels to nodes
V(patient_graph)$cluster <- membership(patient_clusters)

# Visualize clusters
ggraph(patient_graph, layout = "fr") + 
  geom_edge_link(aes(width = weight), alpha = 0.5) +
  geom_node_point(aes(color = factor(V(patient_graph)$cluster)), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = rainbow(length(unique(V(patient_graph)$cluster))), name = "Patient Clusters") +
  ggtitle("AML Patient Similarity Clusters") +
  theme_void()

# Find the Most Important Patients in the Network

# Compute centrality scores
patient_centrality <- betweenness(patient_graph)

# Identify top 10 most central patients
top_hubs <- sort(patient_centrality, decreasing = TRUE)[1:10]
print(top_hubs)

# Adjust Similarity Threshold to Refine Network

threshold <- 0.6  # Try increasing for stricter similarity connections
adj_matrix <- similarity_matrix > threshold

# Characterize Each Cluster (Identify AML Subtypes)

# Add cluster membership to patient dataset
patient_data$Cluster <- membership(patient_clusters)

# Summarize key features for each cluster
summary_by_cluster <- aggregate(. ~ Cluster, data = patient_data, FUN = mean)

# Print summary
print(summary_by_cluster)

# Count the number of unique communities detected by Louvain clustering
num_communities <- length(unique(membership(patient_clusters)))
print(num_communities)

# Compute modularity score of the detected communities
modularity_score <- modularity(patient_clusters)
print(modularity_score)
