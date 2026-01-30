#Completeness - Treponema

# Load ggplot2 package
library(ggplot2)

# Read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam_trep.csv", header = TRUE)

# Define custom color palette
my_colors <- c("#ffa600", "#ff6361", "#bc5090", "#58508d", "#003f5c")

# Create bubble plot with custom colors
ggplot(df, aes(x = Completeness, y = Contamination, size = N50, color = Completeness)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(7, 10)) +
  scale_color_gradient(low = "#009AE0", high = my_colors[5], guide = "colorbar") +
  labs(x = "Completeness", y = "Contamination", size = "N50", color = "Completeness") +
  theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        legend.key.size = unit(1.2, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 5)) +
  scale_x_continuous(limits = c(0, 100))

