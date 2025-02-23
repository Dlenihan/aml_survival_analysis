# **Split the dataset into FLT3+ and FLT3- subgroups**
flt3_positive <- network_data[network_data$FLT3_Mut_Status == 1, ]
flt3_negative <- network_data[network_data$FLT3_Mut_Status == 0, ]

# **Function to Remove Constant Columns**
remove_constant_columns <- function(df) {
  df <- df[, sapply(df, function(x) length(unique(x)) > 1)]  # Keep only features with more than one unique value
  if (ncol(df) < 2) {
    stop("Error: Not enough variables left after filtering.")
  }
  return(df)
}

# Apply the function to FLT3+ and FLT3- groups
flt3_positive <- remove_constant_columns(flt3_positive)
flt3_negative <- remove_constant_columns(flt3_negative)

# **Convert Continuous Variables to Discrete Categories (Binning)**
for (df in list(flt3_positive, flt3_negative)) {
  if ("Age" %in% colnames(df)) df$Age <- as.integer(discretize(df$Age, disc = "equalfreq", nbins = 4)$X)
  if ("BM_Blast_Percentage" %in% colnames(df)) df$BM_Blast_Percentage <- as.integer(discretize(df$BM_Blast_Percentage, disc = "equalfreq", nbins = 4)$X)
  if ("WBC" %in% colnames(df)) df$WBC <- as.integer(discretize(df$WBC, disc = "equalfreq", nbins = 4)$X)
}

# **Function to Build and Plot the Network**
build_network <- function(data, title) {
  
  # **Ensure there are at least two variables before computing correlations**
  if (ncol(data) < 2) {
    print(paste("Skipping network plot for", title, "due to insufficient variables."))
    return(NULL)
  }
  
  # Compute **Mutual Information Matrix** on Discretized Data
  mi_matrix <- mutinformation(data)
  
  # Compute **Pearson Correlation Matrix (if enough variables exist)**
  cor_matrix <- NULL
  if (ncol(data) > 1) {
    cor_matrix <- cor(data, use = "complete.obs")
  }
  
  # **Hybrid Network Approach (MI + Correlation)**
  adj_matrix <- (mi_matrix > 0.05)
  if (!is.null(cor_matrix)) {
    adj_matrix <- adj_matrix | (abs(cor_matrix) > 0.2)
  }
  
  # Convert adjacency matrix to an igraph object
  graph_obj <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  
  # Remove self-loops
  graph_obj <- simplify(graph_obj, remove.loops = TRUE)
  
  # **Fix: Remove edges with NaN weights**
  E(graph_obj)$weight[is.nan(E(graph_obj)$weight)] <- 0
  
  # Assign node labels (feature names)
  V(graph_obj)$name <- colnames(data)
  
  # Set node colors
  V(graph_obj)$color <- "orange"
  
  # Set edge weights (thicker lines for stronger relationships)
  E(graph_obj)$width <- E(graph_obj)$weight * 5
  
  # **Plot the Network**
  plot(graph_obj, 
       layout = layout_with_kk(graph_obj),  # Kamada-Kawai layout
       vertex.size = 15, 
       vertex.label.color = "black", 
       vertex.label.cex = 1.2,
       edge.color = "gray50",
       main = title)
}

# **Plot FLT3+ and FLT3- Networks**
par(mfrow = c(1, 2))  # Arrange plots side by side

build_network(flt3_positive, "FLT3+ AML Feature Network")
build_network(flt3_negative, "FLT3- AML Feature Network")
