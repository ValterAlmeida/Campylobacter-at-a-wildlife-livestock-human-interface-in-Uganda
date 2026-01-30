#Completeness

#Comple_contam_after.csv

#Completeness and contamination

#Strain heterogeneity

head(df)
ncol(df)
nrow(df)

library(ggplot2)

# read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam_campy.csv", header = TRUE)

# create scatterplot
ggplot(df, aes(x=Completeness, y=Contamination, color=Contamination)) + 
  geom_point(size=3, alpha=0.6) +
  scale_color_gradient(low = "#0072B2", high = "#FF4500") +
  xlab("Completeness (%)") + 
  ylab("Contamination (%)") +
  ggtitle("Mags - Completeness and Contamination") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 10, b = 20, l = 0, r = 0, "pt"),
        axis.line = element_line(linetype = "solid", color = "black", 
                                 lineend = "butt", linewidth = 0.5),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10),
        plot.tag = element_text(size = 10, hjust = 1))

#Test2

library(ggplot2)

# read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam_campy.csv", header = TRUE)

# create scatterplot
ggplot(df, aes(x=Completeness, y=Contamination, color=Contamination)) + 
  geom_point(size=3, alpha=0.6) +
  scale_color_gradient(low = "#0072B2", high = "#FF4500") +
  xlab("Completeness (%)") + 
  ylab("Contamination (%)") +
  ggtitle("Mags - Completeness and Contamination") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 10, b = 20, l = 0, r = 0, "pt"),
        axis.line = element_line(linetype = "solid", color = "black", 
                                 lineend = "butt", linewidth = 0.5),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10),
        plot.tag = element_text(size = 10, hjust = 1)) +
  scale_y_continuous(limits = c(0, 100))

#Test3

ggplot(df, aes(x=Completeness, y=Contamination, color=Completeness)) + 
  geom_point(size=3, alpha=0.6) +
  scale_color_gradient(low = "#0072B2", high = "#FFD700") +
  xlab("Completeness (%)") + 
  ylab("Contamination (%)") +
  ggtitle("Mags - Completeness and Contamination") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 10, b = 20, l = 0, r = 0, "pt"),
        axis.line = element_line(linetype = "solid", color = "black", 
                                 lineend = "butt", linewidth = 0.5),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10),
        plot.tag = element_text(size = 10, hjust = 1))
scale_y_continuous(limits = c(0, 100))

#Test4

ggplot(df, aes(x=Completeness, y=Contamination, color=Completeness)) + 
  geom_point(size=3, alpha=0.6) +
  scale_color_gradient(low = "#E18936", high = "#0072B2") + 
  xlab("Completeness (%)") + 
  ylab("Contamination (%)") +
  ggtitle("MAGs - Completeness and Contamination") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + # add this line
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 10, b = 20, l = 0, r = 0, "pt"),
        axis.line = element_line(linetype = "solid", color = "black", 
                                 lineend = "butt", linewidth = 0.5),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10),
        plot.tag = element_text(size = 10, hjust = 1))

# Bublle plot

# Load ggplot2 package
library(ggplot2)

# Read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam_campy.csv", header = TRUE)

# Create bubble plot
ggplot(df, aes(x = Completeness, y = Contamination, size = N50)) +
  geom_point(alpha = 0.5) +
  scale_size(range = c(2, 10)) +
  labs(x = "Completeness", y = "Contamination", size = "N50")


# TEST COLOURS

# Load ggplot2 package
library(ggplot2)

# Read in data
df <- read.csv("/path/to/your/data.csv", header = TRUE)

# Define custom color palette
my_colors <- c("#ffa600", "#ff6361", "#bc5090", "#58508d", "#003f5c")

# Create bubble plot with custom colors
ggplot(df, aes(x = Completeness, y = Contamination, size = N50, color = N50)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(2, 10)) +
  scale_color_gradient(low = "#60C8F8", high = my_colors[5], guide = "colorbar") +
  labs(x = "Completeness", y = "Contamination", size = "N50", color = "N50") +
  theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        legend.key.size = unit(1.2, "cm"), legend.position = "bottom")

#Test

# Load ggplot2 package
library(ggplot2)

# Read in data
df <- read.csv("/path/to/your/data.csv", header = TRUE)

# Define custom color palette
my_colors <- c("#ffa600", "#ff6361", "#60C8F8", "#58508d", "#003f5c")

# Create bubble plot with custom colors
ggplot(df, aes(x = Completeness, y = Contamination, size = N50, color = Contamination)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(2, 10)) +
  scale_color_gradient(low = "#003f5c", high = my_colors[3], guide = "colorbar") +
  labs(x = "Completeness", y = "Contamination", size = "N50", color = "Contamination") +
  theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        legend.key.size = unit(1.2, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 5)) +
  scale_x_continuous(limits = c(0, 100))


#Test 2 (The best one)

# Load ggplot2 package
library(ggplot2)

# Read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam_campy.csv", header = TRUE)

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

















