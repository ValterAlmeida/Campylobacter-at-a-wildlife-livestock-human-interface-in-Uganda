#TESTTTTT

# Load required packages
library(vegan)
library(ggplot2)
library(viridis)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_3_clean.csv", header = TRUE, row.names = 1)
dissimilarity <- as.dist(df)

# Compute the PCoA
pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# Define the groups
sample_groups <- rep(NA, nrow(pcoa_df))
sample_groups[01:81] <- "Congo - Human"
sample_groups[82:86] <- "Uganda - HC - AGRF"
sample_groups[87:91] <- "Uganda - HH - AGRF"
sample_groups[92:137] <- "Italy - Human"
sample_groups[138:237] <- "US - Human - hospital"
sample_groups[238:257] <- "Mexico - Human"

pcoa_df$Group <- sample_groups

# Set custom colors for each group
group_colors <- c("Uganda - HH - AGRF" = "#3d85c6", "Uganda - HC - AGRF" = "#cfe2f3", "Mexico - Human" = "#E8ED52", "Congo - Human" = "#8488BA", "Italy - Human" = "#C36767", "US - Human - hospital" = "lightgrey")

# Change the order of the legend to match the order of the groups
group_order <- c("Congo - Human", "Italy - Human", "Mexico - Human", "US - Human - hospital", "Uganda - HH - AGRF", "Uganda - HC - AGRF")

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

