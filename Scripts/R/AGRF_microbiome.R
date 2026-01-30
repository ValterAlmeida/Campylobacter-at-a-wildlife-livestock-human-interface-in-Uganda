#AGRF - Microbiome

df <- read.csv("/Users/valter/Downloads//Kraken2/AGRF_kraken2_3.csv",
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
my_colors <- c("#1f78b4", "#6D058E", "#e31a1c", "#fdbf6f", "#33a02c", "#000000", "#a6cee3")

# Reorder levels of the Host variable based on their row order in the original data frame
df_long$Host <- factor(df_long$Host, levels = unique(df$Host))

# Create stacked bar plot without spaces between bars
ggplot(df_long, aes(x = Host, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8, position = "fill") +
  scale_fill_manual(values = my_colors, name = "Phylum") +
  xlab("Host") +
  ylab("Relative Abundance") +
  ggtitle("Microbiome Composition by Host (pooled) - Uganda") +
  theme_classic() +
  labs(fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1.2),
        axis.title.x = element_text(margin = margin(t = 10)))
