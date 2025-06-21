#--------------------------------With simulation_data_Scenario_1.csv----------------------------------


#--------------------------------PCA Analysis----------------------------------

# Load required packages
library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(mixOmics)
library(dplyr)

# Read the CSV file (replace with actual path)
dfTOT <- read.csv("PATH_TO_YOUR_FILE/simulation_data_Scenario_1.csv", header = TRUE, sep = ",")

# Remove the first two columns (e.g., IDs or labels) and scale the data
rmvec <- c(1,2)
X <- scale(dfTOT[-rmvec])
rm(rmvec)

# Run PCA
pca_result <- prcomp(X, center = FALSE, scale. = FALSE)

# Get eigenvalues
eig_val <- get_eigenvalue(pca_result)

# Scree plot of eigenvalues
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "Scree Plot", x = "Principal Components", y = "Explained Variance (%)")

# Get variable information
var <- get_pca_var(pca_result)

# Sort variables by contribution to PC1
var_pca1_sort <- sort(var$contrib[,1], decreasing = TRUE)

# Visualize PCA variables (colored by contribution)
fviz_pca_var(pca_result, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) +
  labs(title = "PCA - Variable Contributions")

#--------------------------------PLS-DA Analysis----------------------------------

# Create a categorical grouping variable based on metabolite production
median_metabolites <- median(dfTOT$Total_Metabolites_Produced)

dfTOT$Production_Group <- ifelse(dfTOT$Total_Metabolites_Produced >= median_metabolites,
                                 "High", "Low")
dfTOT$Production_Group <- as.factor(dfTOT$Production_Group)

# Prepare data for PLS-DA
X <- dfTOT %>%
  dplyr::select(-Step, -Time_stamp, -Total_Metabolites_Produced, -Production_Group)

Y <- dfTOT$Production_Group

# Run PLS-DA with 3 components
plsda_result <- plsda(X, Y, ncomp = 3)

# Plot individuals (samples) in PLS-DA space
plotIndiv(plsda_result, comp = c(1,2), group = Y,
          ind.names = FALSE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA: Individuals (Comp 1 vs 2)")

# Plot variable contributions (loadings)
plotVar(plsda_result, comp = c(1,2),
        title = "PLS-DA: Variable Contributions")
