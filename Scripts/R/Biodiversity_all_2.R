#Beta diversity - All

# read in data

df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

head(df)
nrow(df)
ncol(df)

# Load the vegan package
library(vegan)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

# Check for missing values in the distance matrix
any(is.na(df))

# Remove missing values from the distance matrix
df_clean <- na.omit(df)

# Convert the distance matrix to a distance object
dist_matrix <- dist(df_clean)

# Perform a PCoA analysis
pcoa <- cmdscale(dist_matrix, k = nrow(df_clean) - 1)

# Create a data frame containing the PCoA scores
pcoa_df <- data.frame(Sample = rownames(pcoa), PC1 = pcoa[, 1], PC2 = pcoa[, 2])

# Create a scatter plot of the PCoA scores
plot(pcoa_df$PC1, pcoa_df$PC2, pch = 16, cex = 1.5, col = "black", xlab = "PC1", ylab = "PC2")

# Add labels to the plot
text(pcoa_df$PC1, pcoa_df$PC2, labels = pcoa_df$Sample, pos = 4, cex = 0.7)



#ADD GROUPS TEST 1

# Load the vegan package
library(vegan)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

# Check for missing values in the distance matrix
any(is.na(df))

# Remove missing values from the distance matrix
df_clean <- na.omit(df)

# Convert the distance matrix to a distance object
dist_matrix <- dist(df_clean)

# Perform a PCoA analysis
pcoa <- cmdscale(dist_matrix, k = nrow(df_clean) - 1)

# Define the sample groups
Congo_captive_gorilla <- 1:16
Congo_gorilla <- 17:35
Congo_human <- 36:116
Congo_wild_gorilla <- 117:144
Uganda_HC_AGRF <- 145:149
Uganda_HH_AGRF <- 150:154
Italy_human <- 155:200
US_human_hospital <- 201:300
Mexico_human <- 301:320
Uganda_wild_gorilla <- 321:418
Uganda_habituated_gorilla <- 419:518
Uganda_human_hospital <- 519:616
Uganda_human_community <- 617:713

# Create a vector of colors for each group
group_colors <- c("red", "blue", "green", "orange", "purple", "pink", "yellow", "lightblue", "black", "grey")

# Create a scatter plot of the PCoA scores, colored by group
plot(pcoa_df$PC1, pcoa_df$PC2, pch = 16, cex = 1.5, col = group_colors[pcoa_df$group], 
     xlab = "PC1", ylab = "PC2")

# Add labels to the plot
text(pcoa_df$PC1, pcoa_df$PC2, labels = pcoa_df$Sample, pos = 4, cex = 0.7)

# Add legend to the plot
legend("topright", legend = levels(pcoa_df$group), col = group_colors, pch = 16)


#TEST 2

# Load the vegan package
library(vegan)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

# Check for missing values in the distance matrix
any(is.na(df))

# Remove missing values from the distance matrix
df_clean <- na.omit(df)

# Convert the distance matrix to a distance object
dist_matrix <- dist(df_clean)

# Perform a PCoA analysis
pcoa <- cmdscale(dist_matrix, k = nrow(df_clean) - 1)

# Create a scatter plot of the PCoA scores, colored by group
plot(pcoa_df$PC1, pcoa_df$PC2, pch = 16, cex = 1.5, col = group_colors[pcoa_df$group], 
     xlab = "PC1", ylab = "PC2")

# Add labels to the plot
text(pcoa_df$PC1, pcoa_df$PC2, labels = pcoa_df$Sample, pos = 4, cex = 0.7)

# Add legend to the plot
legend("topright", legend = levels(pcoa_df$group), col = group_colors, pch = 16)


#TEST 3

# Load required libraries
library(ggplot2)
library(vegan)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

# Check for missing values in the distance matrix
any(is.na(df))

# Remove missing values from the distance matrix
df_clean <- na.omit(df)

# Convert the distance matrix to a distance object
dist_matrix <- dist(df_clean)

# Perform a PCoA analysis
pcoa <- cmdscale(dist_matrix, k = nrow(df_clean) - 1)

# Create a ggplot object with PC1 and PC2 on x and y axes, respectively
ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  
  # Add points colored by group
  geom_point(aes(color = group), size = 3.5) +
  
  # Add labels to the points
  geom_text(aes(label = Sample), size = 3, hjust = 0.5, vjust = 0.5) +
  
  # Add dissimilarity values as tooltips
  ggplot2::tooltip_aes("dissimilarity") +
  
  # Add a title to the plot
  ggtitle("PCoA Plot") +
  
  # Customize the legend
  scale_color_manual(values = group_colors, name = "Group") +
  
  # Customize the x and y axis labels
  xlab("PC1") +
  ylab("PC2") +
  
  # Customize the theme of the plot
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right")


# Heatmap:

# Load required libraries
library(ggplot2)
library(reshape2)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

# Check for missing values in the distance matrix
any(is.na(df))

# Remove missing values from the distance matrix
df_clean <- na.omit(df)

# Convert the distance matrix to a distance object
dist_matrix <- dist(df_clean)

# Create a matrix from the distance object
dist_matrix_mat <- as.matrix(dist_matrix)

# Create a rownames column for the matrix
dist_matrix_mat <- cbind(rownames(dist_matrix_mat), dist_matrix_mat)

# Convert the matrix to a long format
dist_matrix_long <- melt(dist_matrix_mat, id.vars = "rownames")

# Rename the columns
colnames(dist_matrix_long) <- c("Sample1", "Sample2", "Distance")

# Convert the Distance column to a numeric data type
dist_matrix_long$Distance <- as.numeric(dist_matrix_long$Distance)

# Create a ggplot object with Sample1 and Sample2 on x and y axes, respectively
ggplot(dist_matrix_long, aes(x = Sample1, y = Sample2, fill = Distance)) +
  
  # Add the heatmap
  geom_tile() +
  
  # Customize the color scale
  scale_fill_gradient(low = "white", high = "red", name = "Distance") +
  
  # Rotate the x-axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  
  # Add a title to the plot
  ggtitle("Distance Matrix Heatmap")


# TEST Heatmap 2

# Load required libraries
library(ggplot2)
library(reshape2)

# Define the groups for each sample
groups <- c(rep("Wild gorilla", 98),
            rep("Habituated gorilla", 100),
            rep("Human - hospital", 98),
            rep("Human - community", 97),
            rep("Cattle", 83),
            rep("Goat", 73))

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE)

# Check for missing values in the distance matrix
any(is.na(df))

# Remove missing values from the distance matrix
df_clean <- na.omit(df)

# Convert the distance matrix to a distance object
dist_matrix <- dist(df_clean)

# Create a matrix from the distance object
dist_matrix_mat <- as.matrix(dist_matrix)

# Create a rownames column for the matrix
dist_matrix_mat <- cbind(rownames(dist_matrix_mat), dist_matrix_mat)

# Convert the matrix to a long format
dist_matrix_long <- melt(dist_matrix_mat, id.vars = "rownames")

# Rename the columns
colnames(dist_matrix_long) <- c("Sample1", "Sample2", "Distance")

# Convert the Distance column to a numeric data type
dist_matrix_long$Distance <- as.numeric(dist_matrix_long$Distance)

# Add a column for the groups
dist_matrix_long$Group <- groups

# Create a ggplot object with Sample1 and Sample2 on x and y axes, respectively
ggplot(dist_matrix_long, aes(x = Sample1, y = Sample2, fill = Distance)) +
  
  # Add the heatmap
  geom_tile() +
  
  # Color the tiles by group
  scale_fill_gradient(low = "white", high = "red", name = "Distance") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(dist_matrix_long$Sample1))) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(. ~ Group) +
  
  # Rotate the x-axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  
  # Add a title to the plot
  ggtitle("Distance Matrix Heatmap")


# TEST 3 PCoA and Heatmap: PERFECT

library(vegan)
library(ggplot2)
library(RColorBrewer)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE, row.names = 1)
dissimilarity <- as.dist(df)

# Compute the PCoA
pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# Define the groups
sample_groups <- rep(NA, nrow(pcoa_df))
sample_groups[1:16] <- "Congo - Captive gorilla"
sample_groups[17:35] <- "Congo - Gorilla"
sample_groups[36:116] <- "Congo - Human"
sample_groups[117:144] <- "Congo - Wild gorilla"
sample_groups[145:149] <- "Uganda - HH - AGRF"
sample_groups[150:154] <- "Uganda - HC - AGRF"
sample_groups[155:200] <- "Italy - Human"
sample_groups[201:300] <- "US - Human - hospital"
sample_groups[301:320] <- "Mexico - Human"
sample_groups[321:418] <- "Uganda - Wild gorilla"
sample_groups[419:518] <- "Uganda - Habituated gorilla"
sample_groups[519:616] <- "Uganda - Human - hospital"
sample_groups[617:713] <- "Uganda - Human - community"

pcoa_df$Group <- sample_groups

# Create the PCoA plot with groups
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "yellow", "gold", "purple", "orange", "black", "pink", "lightblue", "darkgreen", "grey", "black")) +
  ggtitle("PCoA Plot with Groups") +
  xlab("PC1") + ylab("PC2") +
  theme(plot.title = element_text(hjust = 0.5))

# Create a heatmap
heatmap_data <- as.matrix(df)
heatmap_data_scale <- scale(heatmap_data)
heatmap_data_dist <- dist(heatmap_data_scale)
heatmap_colors <- brewer.pal(n = 9, name = "YlOrRd")
heatmap_breaks <- seq(from = min(heatmap_data), to = max(heatmap_data), length.out = length(heatmap_colors) + 1)

ggplot(data.frame(as.table(heatmap_data)), aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = rev(heatmap_colors), breaks = heatmap_breaks, na.value = "white") +
  ggtitle("Distance Matrix Heatmap") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_fixed(ratio = 1)


#TRYING TO FIX THE CODE:

# TEST 3 PCoA and Heatmap: PERFECT

library(vegan)
library(ggplot2)
library(RColorBrewer)

# Try reinstalling the gtable package
if (!requireNamespace("gtable", quietly = TRUE)) {
  install.packages("gtable")
}

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE, row.names = 1)
dissimilarity <- as.dist(df)

# Compute the PCoA
pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# Define the groups
sample_groups <- rep(NA, nrow(pcoa_df))
sample_groups[1:16] <- "Congo - Captive gorilla"
sample_groups[17:35] <- "Congo - Gorilla"
sample_groups[36:116] <- "Congo - Human"
sample_groups[117:144] <- "Congo - Wild gorilla"
sample_groups[145:149] <- "Uganda - HH - AGRF"
sample_groups[150:154] <- "Uganda - HC - AGRF"
sample_groups[155:200] <- "Italy - Human"
sample_groups[201:300] <- "US - Human - hospital"
sample_groups[301:320] <- "Mexico - Human"
sample_groups[321:418] <- "Uganda - Wild gorilla"
sample_groups[419:518] <- "Uganda - Habituated gorilla"
sample_groups[519:616] <- "Uganda - Human - hospital"
sample_groups[617:713] <- "Uganda - Human - community"

pcoa_df$Group <- sample_groups

# Create the PCoA plot with groups
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "darkblue", "darkgrey", "green", "purple", "orange", "black", "pink", "lightblue", "darkgreen", "grey", "black")) +
  ggtitle("PCoA Plot with Groups") +
  xlab("PC1") + ylab("PC2") +
  theme(plot.title = element_text(hjust = 0.5))





# TESTE de novo

library(ggplot2)

# Create a data frame with PCoA scores and sample groups
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2], Group = sample_groups)

# Create a ggplot object with PC1 and PC2 on x and y axes, respectively, colored by group
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  # Add points for each sample
  geom_point(alpha = 0.5) +
  # Add labels for each sample
  geom_text(aes(label = Sample), vjust = -1, hjust = 0.5) +
  # Add a title to the plot
  ggtitle("PCoA Plot with Transparency") +
  # Adjust color scale
  scale_color_discrete(name = "Group") +
  # Adjust transparency
  theme(legend.key.size = unit(0.5, "cm"))

#TEST 5

library(ggplot2)
library(viridis)

# Create a data frame with PC1, PC2, and group information
pcoa_df <- data.frame(Sample = rownames(pcoa$points), 
                      PC1 = pcoa$points[, 1], 
                      PC2 = pcoa$points[, 2], 
                      Group = sample_groups)

# Create a ggplot object with PC1 and PC2 on x and y axes, respectively
ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  
  # Add points to the plot, color-coded by group
  geom_point(aes(color = Group), size = 3, alpha = 0.7) +
  
  # Add color scale to the legend
  scale_color_viridis(discrete = TRUE) +
  
  # Add a title to the plot
  ggtitle("PCoA Plot of Distance Matrix") +
  
  # Add axis labels
  xlab("PC1") + ylab("PC2") +
  
  # Adjust theme
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


# TEST ^^

# Set custom colors for each group
group_colors <- c("Uganda - HH - AGRF" = "darkblue", "Uganda - HC - AGRF" = "#001BB2", "Mexico - Human" = "#E8ED52", "Congo - Captive gorilla" = "#A6FF70", "Congo - Gorilla" = "#8AE552", "Congo - Wild gorilla" = "#48BE00", "Congo - Human" = "#8488BA", "Italy - Human" = "#C36767", "US - Human hospital" = "#004AF5", "Uganda - Wild gorilla" = "#38761d", "Uganda - Habituated gorilla" = "#b6d7a8", 
                  "Uganda - Human hospital" = "#3d85c6", "Uganda - Human community" = "#cfe2f3")

# Change the order of the legend to match the order of the groups
group_order <- c("Congo - Captive gorilla", "Congo - Gorilla", "Congo - Wild gorilla", 
                 "Congo - Human", "Italy - Human", "Mexico - Human", "US - Human hospital", "Uganda - Wild gorilla", 
                 "Uganda - Habituated gorilla", "Uganda - Human hospital", "Uganda - Human community", "Uganda - HH - AGRF", "Uganda - HC - AGRF")

# Create a ggplot object with PC1 and PC2 on x and y axes, respectively
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  # Add points with custom colors and alpha value
  geom_point(alpha = 0.8, size = 3, aes(color = Group)) +
  scale_color_manual(values = group_colors, limits = group_order) +
  # Add a title to the plot
  ggtitle("PCoA - Beta diversity") +
  # Add x and y axis labels
  xlab("PC1") + ylab("PC2") +
  # Adjust the theme of the plot
  theme_bw()


#TEST#$%^^

# TEST 3 PCoA and Heatmap: PERFECT

library(vegan)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_2_clean.csv", header = TRUE, row.names = 1)
dissimilarity <- as.dist(df)

# Compute the PCoA
pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# Define the groups
sample_groups <- rep(NA, nrow(pcoa_df))
sample_groups[1:16] <- "Congo - Captive gorilla"
sample_groups[17:35] <- "Congo - Gorilla"
sample_groups[36:116] <- "Congo - Human"
sample_groups[117:144] <- "Congo - Wild gorilla"
sample_groups[145:149] <- "Uganda - HH - AGRF"
sample_groups[150:154] <- "Uganda - HC - AGRF"
sample_groups[155:200] <- "Italy - Human"
sample_groups[201:300] <- "US - Human - hospital"
sample_groups[301:320] <- "Mexico - Human"
sample_groups[321:418] <- "Uganda - Wild gorilla"
sample_groups[419:518] <- "Uganda - Habituated gorilla"
sample_groups[519:616] <- "Uganda - Human - hospital"
sample_groups[617:713] <- "Uganda - Human - community"

pcoa_df$Group <- sample_groups

# Create the PCoA plot with groups
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = viridis(n = length(unique(sample_groups)), option = "D")) +
  ggtitle("PCoA Plot with Groups") +
  xlab("PC1") + ylab("PC2") +
  theme(plot.title = element_text(hjust = 0.5))


#TESTTTTT

# Load required packages
library(vegan)
library(ggplot2)
library(viridis)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_all_4_clean.csv", header = TRUE, row.names = 1)
dissimilarity <- as.dist(df)

# Compute the PCoA
pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# Define the groups
sample_groups <- rep(NA, nrow(pcoa_df))
sample_groups[1:16] <- "Congo - Captive gorilla"
sample_groups[17:35] <- "Congo - Gorilla"
sample_groups[36:116] <- "Congo - Human"
sample_groups[117:144] <- "Congo - Wild gorilla"
sample_groups[145:149] <- "Uganda - GH - AGRF"
sample_groups[150:154] <- "Uganda - GW - AGRF"
sample_groups[155:159] <- "Uganda - HC - AGRF"
sample_groups[160:164] <- "Uganda - HH - AGRF"
sample_groups[165:210] <- "Italy - Human"
sample_groups[211:310] <- "US - Human - hospital"
sample_groups[311:330] <- "Mexico - Human"
sample_groups[331:428] <- "Uganda - Wild gorilla"
sample_groups[429:528] <- "Uganda - Habituated gorilla"
sample_groups[529:626] <- "Uganda - Human - hospital"
sample_groups[627:723] <- "Uganda - Human - community"

pcoa_df$Group <- sample_groups

# Set custom colors for each group
group_colors <- c("Uganda - HH - AGRF" = "black", "Uganda - HC - AGRF" = "red", "Mexico - Human" = "#E8ED52", "Congo - Captive gorilla" = "#A6FF70", "Congo - Gorilla" = "#8AE552", "Congo - Wild gorilla" = "#48BE00", "Congo - Human" = "#8488BA", "Italy - Human" = "#C36767", "US - Human - hospital" = "lightgrey", "Uganda - Wild gorilla" = "#38761d", "Uganda - Habituated gorilla" = "#b6d7a8", 
                  "Uganda - Human - hospital" = "#3d85c6", "Uganda - Human - community" = "#cfe2f3", "Uganda - GH - AGRF" = "#00F0DA", "Uganda - GW - AGRF" = "#00BFAD")

# Change the order of the legend to match the order of the groups
group_order <- c("Congo - Captive gorilla", "Congo - Gorilla", "Congo - Wild gorilla", 
                 "Congo - Human", "Italy - Human", "Mexico - Human", "US - Human - hospital", "Uganda - Wild gorilla", 
                 "Uganda - Habituated gorilla", "Uganda - Human - hospital", "Uganda - Human - community", "Uganda - HH - AGRF", "Uganda - HC - AGRF", "Uganda - GH - AGRF", "Uganda - GW - AGRF")

# Create the PCoA plot with groups and custom colors
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors) +
  ggtitle("PCoA Plot with Groups") +
  xlab("PC1") + ylab("PC2") +
  theme(plot.title = element_text(hjust = 0.5))

# Create the PCoA plot with groups and custom colors, with transparency
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) + # Add alpha parameter for transparency
  scale_color_manual(values = group_colors) +
  ggtitle("PCoA Plot with Groups") +
  xlab("PC1") + ylab("PC2") +
  theme(plot.title = element_text(hjust = 0.5))

# Create a ggplot object with PC1 and PC2 on x and y axes, respectively
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  # Add points with custom colors and alpha value
  geom_point(alpha = 0.6, size = 3, aes(color = Group)) +
  scale_color_manual(values = group_colors, limits = group_order) +
  # Add a title to the plot
  ggtitle("PCoA - Beta diversity") +
  # Add x and y axis labels
  xlab("PC1") + ylab("PC2") +
  # Adjust the theme of the plot
  theme_bw()
  

