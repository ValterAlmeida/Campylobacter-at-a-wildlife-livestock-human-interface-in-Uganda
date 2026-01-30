#Microbiome_comparison_all

df <- read.csv("/Users/valter/Downloads/R_Microbiome/Comparisons_Orion_3_stats.csv",
               header = TRUE)

library(ggplot2)
library(reshape2)
library(tidyr) # load tidyr package

head(df)
ncol(df)
nrow(df)

# Convert data from wide to long format
df_long <- reshape2::melt(df, id.vars = c("Origin", "Host"))

# Define custom color palette
my_colors <- c("#1f78b4", "#6D058E", "#e31a1c", "#fdbf6f", "#33a02c", "#000000", "grey")

# Reorder levels of the Host variable based on their row order in the original data frame
df_long$Host <- factor(df_long$Host, levels = unique(df$Host))

# Create stacked bar plot without spaces between bars
ggplot(df_long, aes(x = Host, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.9, position = "fill") +
  scale_fill_manual(values = my_colors, name = "Phylum") +
  xlab("Host") +
  ylab("Relative Abundance") +
  ggtitle("Human - Microbiome Composition") +
  theme_classic() +
  labs(fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1.2),
        axis.title.x = element_text(margin = margin(t = 10)))



#STATS:

#Wild gorillas x habituateds gorillas:

#Firmicutes

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Wild_gorilla" & df_long$variable == "Firmicutes"]
group2 <- df_long$value[df_long$Host == "Habituated_gorilla" & df_long$variable == "Firmicutes"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Actinobacteria

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Wild_gorilla" & df_long$variable == "Actinobacteria"]
group2 <- df_long$value[df_long$Host == "Habituated_gorilla" & df_long$variable == "Actinobacteria"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Proteobacteria

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Wild_gorilla" & df_long$variable == "Proteobacteria"]
group2 <- df_long$value[df_long$Host == "Habituated_gorilla" & df_long$variable == "Proteobacteria"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Bacteroidetes

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Wild_gorilla" & df_long$variable == "Bacteroidetes"]
group2 <- df_long$value[df_long$Host == "Habituated_gorilla" & df_long$variable == "Bacteroidetes"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Spirochaetes

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Wild_gorilla" & df_long$variable == "Spirochaetes"]
group2 <- df_long$value[df_long$Host == "Habituated_gorilla" & df_long$variable == "Spirochaetes"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Verrucomicrobia

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Wild_gorilla" & df_long$variable == "Verrucomicrobia"]
group2 <- df_long$value[df_long$Host == "Habituated_gorilla" & df_long$variable == "Verrucomicrobia"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Humans in hospital x Humans in the community:

#Firmicutes

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Human_hospital" & df_long$variable == "Firmicutes"]
group2 <- df_long$value[df_long$Host == "Human_community" & df_long$variable == "Firmicutes"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Actinobacteria

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Human_hospital" & df_long$variable == "Actinobacteria"]
group2 <- df_long$value[df_long$Host == "Human_community" & df_long$variable == "Actinobacteria"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Proteobacteria

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Human_hospital" & df_long$variable == "Proteobacteria"]
group2 <- df_long$value[df_long$Host == "Human_community" & df_long$variable == "Proteobacteria"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Bacteroidetes

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Human_hospital" & df_long$variable == "Bacteroidetes"]
group2 <- df_long$value[df_long$Host == "Human_community" & df_long$variable == "Bacteroidetes"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Spirochaetes

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Human_hospital" & df_long$variable == "Spirochaetes"]
group2 <- df_long$value[df_long$Host == "Human_community" & df_long$variable == "Spirochaetes"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

#Verrucomicrobia

# Subset the data for the two groups of interest
group1 <- df_long$value[df_long$Host == "Human_hospital" & df_long$variable == "Verrucomicrobia"]
group2 <- df_long$value[df_long$Host == "Human_community" & df_long$variable == "Verrucomicrobia"]

# Perform the Wilcoxon/Mann-Whitney U test
result <- wilcox.test(group1, group2)

# Print the test results
cat("Wilcoxon/Mann-Whitney U test:\n")
cat(paste("W =", round(result$statistic, 3), "\n"))
p_value <- formatC(result$p.value, digits = 5, format = "f")
cat(paste("p-value <", p_value, "\n"))

