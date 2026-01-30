# TEST 3 PCoA and Heatmap: PERFECT

library(vegan)
library(ggplot2)
library(RColorBrewer)

# Read in the distance matrix
df <- read.csv("/Users/valter/Downloads/Kraken2/beta_diversity_uganda_AGRF_clean.csv", header = TRUE, row.names = 1)
dissimilarity <- as.dist(df)

# Compute the PCoA
pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# Define the groups
sample_groups <- rep(NA, nrow(pcoa_df))
sample_groups[1:5] <- "Wild gorilla"
sample_groups[6:10] <- "Habituated gorilla"
sample_groups[11:15] <- "Human - hospital"
sample_groups[16:20] <- "Human - community"
sample_groups[21:25] <- "Cattle"
sample_groups[26:30] <- "Goat"

pcoa_df$Group <- sample_groups

# Create the PCoA plot with groups
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "black")) +
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

# Set custom colors for each group
group_colors <- c("Wild gorilla" = "#38761d", "Habituated gorilla" = "#b6d7a8", 
                  "Human - hospital" = "#3d85c6", "Human - community" = "#cfe2f3", 
                  "Cattle" = "#e69138", "Goat" = "#ffd966")

# Change the order of the legend to match the order of the groups
group_order <- c("Wild gorilla", "Habituated gorilla", "Human - hospital", 
                 "Human - community", "Cattle", "Goat")

# Create a ggplot object with PC1 and PC2 on x and y axes, respectively
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  # Add points with custom colors and alpha value
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = group_colors, limits = group_order) +
  # Add a title to the plot
  ggtitle("PCoA - Beta diversity") +
  # Add x and y axis labels
  xlab("PC1") + ylab("PC2") +
  # Adjust the theme of the plot
  theme_bw()
