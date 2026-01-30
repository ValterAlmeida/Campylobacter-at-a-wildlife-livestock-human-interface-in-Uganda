# read in data

df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_1.csv", header = TRUE)


#TEST 3

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_1.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Create boxplot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_boxplot(alpha = 0.3, colour = "black") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank()) +
  labs(x = "Group", y = "Shannon's diversity")


#TEST 4 (BEST)

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_1.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Run ANOVA test
anova_result <- aov(Alpha_diversity ~ Group, data = df)

# Display ANOVA test results
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result

# Create boxplot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_boxplot(alpha = 0.7, colour = "black") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.4, color = "black") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "Group", y = "Shannon's diversity")


#NEEEW

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_1.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Run ANOVA test
anova_result <- aov(Alpha_diversity ~ Group, data = df)

# Display ANOVA test results
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result

# Create boxplot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_boxplot(alpha = 0.7, colour = "black") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.4, color = "black", size = 4) + # increased size parameter
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "Group", y = "Shannon's diversity")


#NEEEW colours

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_1.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Run ANOVA test
anova_result <- aov(Alpha_diversity ~ Group, data = df)

# Display ANOVA test results
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result

# Create boxplot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_boxplot(alpha = 0.7, colour = "black") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.4, color = "blue", size = 5) + # changed color parameter to red
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "Group", y = "Shannon's diversity")


#TEST NEWW AGAIN

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_1.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Run ANOVA test
anova_result <- aov(Alpha_diversity ~ Group, data = df)

# Display ANOVA test results
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result

# Create boxplot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_boxplot(alpha = 0.7, colour = "black") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5, size = 4) + # changed color parameter to black and size to 5
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) + # added a new scale_color_manual function with my_colors vector
  theme_classic() +
  labs(x = "Group", y = "Shannon's diversity")




#TEST 5

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_AGRF.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Run ANOVA test
anova_result <- aov(Alpha_diversity ~ Group, data = df)

# Display ANOVA test results
summary(anova_result)

# Compute Tukey test results
tukey_result <- TukeyHSD(anova_result)

# Print Tukey test results to console
tukey_summary <- summary(tukey_result)
tukey_summary

# Create violin plot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "Group", y = "Shannon's diversity")

#TEST 6 - Violin

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_AGRF.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Run ANOVA test
anova_result <- aov(Alpha_diversity ~ Group, data = df)

# Display ANOVA test results
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result

# Create violin plot using ggplot2 with average bars
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_violin(alpha = 0.7, colour = "black") +
  geom_jitter(width = 0.18, height = 0, alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, size = 0.6, colour = "#e06666") +
  theme_classic() +
  labs(x = "Group", y = "Shannon's diversity")

#TEST 7 for the significance letters

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/alpha_diversity_Uganda_AGRF.csv", header = TRUE)

# Load required packages
library(ggplot2)
library(dplyr)

# Define the order of the groups
group_order <- c("Wild_gorilla", "Habituated_gorilla", "Human_hospital", "Human_community", "Cattle", "Goat")

# Reorder the levels of the Group variable
df$Group <- factor(df$Group, levels = group_order)

# Define colours for the boxplot and points
my_colors <- c("#38761d", "#d9ead3", "#3d85c6", "#cfe2f3", "#f6b26b", "#ffd966")

# Create a fancy boxplot using ggplot2
ggplot(df, aes(x = Group, y = Alpha_diversity, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(x = "Group", y = "Shannon's diversity")



