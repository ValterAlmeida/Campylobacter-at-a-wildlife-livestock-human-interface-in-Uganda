#Bar plots - Genera

#gw_agrf_genera_copy.csv

df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/gw_agrf_genera_copy.csv",
               header = TRUE)

#Colours

# Load the ggplot2 library
library(ggplot2)

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/gw_agrf_genera_copy.csv", header = TRUE)

# Subset the data to only include the first 15 columns
df_subset <- df[, 2:16]

# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Define a new color palette
new_colors <- c("#FF7F00", "#E41A1C", "#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FFC20A", "#F781BF", "#A65628", "#F0E442", "#000000", "#999999", "#A6CEE3", "#1F78B4", "#B2DF8A")

# Create the bar plot with the new color palette
ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = column_colors) + # Set the new color palette
  ylim(0, 0.5) + # Set the y-axis limits
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic

#Habituated gorilla

df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/gh_agrf_genera_copy.csv",
               header = TRUE)

#Colours

# Load the ggplot2 library
library(ggplot2)

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/gh_agrf_genera_copy.csv", header = TRUE)

# Subset the data to only include the first 15 columns
df_subset <- df[, 2:16]

# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Define a new color palette
column_names <- c("Olsenella", "Parolsenella", "Collinsella", "Thermophilibacter", "Streptomyces", "Prevotella", "Parafannyhessea", "Faecalibacterium", "Blautia", "Clostridium", "Pseudomonas", "Bacteroides", "Paenibacillus", "Vescimonas", "Bifidobacterium")
new_colors <- c("#F781BF", "#F0E442", "#984EA3", "#B2DF8A", "#1F78B4", "#000000", "#A65628", "#FFC20A", "#377EB8", "#4DAF4A", "#999999", "#FFFF33", "#5684FF", "#2AFF00", "#00FFF3")

# Create the bar plot with the new color palette
ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = new_colors) + # Set the new color palette
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic

#TEST - Matching colours

# Define a new color palette
column_names <- c("Olsenella", "Parolsenella", "Collinsella", "Thermophilibacter", "Streptomyces", "Prevotella", "Parafannyhessea", "Faecalibacterium", "Blautia", "Clostridium", "Pseudomonas", "Bacteroides", "Paenibacillus", "Vescimonas", "Bifidobacterium")
new_colors <- c("#F781BF", "#F0E442", "#984EA3", "#B2DF8A", "#1F78B4", "#000000", "#A65628", "#FFC20A", "#377EB8", "#4DAF4A", "#999999", "#FFFF33", "#5684FF", "#2AFF00", "#00FFF3")

# Define a named vector of column names and colors
column_colors <- c(Olsenella = new_colors[1],
                   Parolsenella = new_colors[2],
                   Collinsella = new_colors[3],
                   Thermophilibacter = new_colors[4],
                   Streptomyces = new_colors[5],
                   Prevotella = new_colors[6],
                   Parafannyhessea = new_colors[7],
                   Faecalibacterium = new_colors[8],
                   Blautia = new_colors[9],
                   Clostridium = new_colors[10],
                   Pseudomonas = new_colors[11],
                   Bacteroides = new_colors[12],
                   Paenibacillus = new_colors[13],
                   Vescimonas = new_colors[14],
                   Bifidobacterium = new_colors[15])

# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Create the bar plot with the new color palette
ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = column_colors) + # Set the new color palette
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic


#Test variation:

# Define a new color palette
new_colors <- c("#F781BF", "#F0E442", "#984EA3", "#B2DF8A", "#1F78B4", "#000000", "#A65628", "#FFC20A", "#377EB8", "#4DAF4A", "#999999", "#FFFF33", "#5684FF", "#2AFF00", "#00FFF3")

# Define a named vector of column names and colors
column_colors <- c(Olsenella = new_colors[1],
                   Parolsenella = new_colors[2],
                   Collinsella = new_colors[3],
                   Thermophilibacter = new_colors[4],
                   Streptomyces = new_colors[5],
                   Prevotella = new_colors[6],
                   Parafannyhessea = new_colors[7],
                   Faecalibacterium = new_colors[8],
                   Blautia = new_colors[9],
                   Clostridium = new_colors[10],
                   Pseudomonas = new_colors[11],
                   Bacteroides = new_colors[12],
                   Paenibacillus = new_colors[13],
                   Vescimonas = new_colors[14],
                   Bifidobacterium = new_colors[15])

# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Create the bar plot with the new color palette
ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = column_colors) + # Set the new color palette
  ylim(0, 0.5) + # Set the y-axis limits
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic


# Human

df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/gw_agrf_genera_copy.csv",
               header = TRUE)

#Colours

# Load the ggplot2 library
library(ggplot2)

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/hh_agrf_genera_copy.csv", header = TRUE)

# Subset the data to only include the first 15 columns
df_subset <- df[, 2:16]
df_subset_complete <- tidyr::complete(df_subset, fill = list(value = 0))


# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Define a new color palette
new_colors <- c("#FF7F00", "#E41A1C", "#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FFC20A", "#F781BF", "#A65628", "#F0E442", "#000000", "#999999", "#A6CEE3", "#1F78B4", "#B2DF8A")

# Define a named vector of column names and colors
column_colors <- c(Prevotella = new_colors[1],
                   Faecalibacterium = new_colors[2],
                   Treponema = new_colors[3],
                   Bifidobacterium = new_colors[4],
                   Blautia = new_colors[5],
                   Bacteroides = new_colors[6],
                   Ruminococcus = new_colors[7],
                   Coprococcus = new_colors[8],
                   Collinsella = new_colors[9],
                   Roseburia = new_colors[10],
                   Mediterraneibacter = new_colors[11],
                   Clostridium = new_colors[12],
                   Anaerobutyricum = new_colors[13],
                   Anaerostipes = new_colors[14],
                   Dorea = new_colors[15])

ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = column_colors) + # Set the new color palette
  ylim(0, 2.5) + # Set the y-axis limits
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic

#TEST

#Human - community

# Load the ggplot2 library
library(ggplot2)

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/hc_agrf_genera_copy.csv", header = TRUE)

# Subset the data to only include the first 15 columns
df_subset <- df[, 2:16]
df_subset_complete <- tidyr::complete(df_subset, fill = list(value = 0))


# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Define a new color palette
new_colors <- c("#FF7F00", "#E41A1C", "#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FFC20A", "#F781BF", "#A65628", "#F0E442", "#000000", "#922B21", "#A6CEE3", "#900C3F", "#2B8884")

# Define a named vector of column names and colors
column_colors <- c(Prevotella = new_colors[1],
                   Faecalibacterium = new_colors[2],
                   Treponema = new_colors[3],
                   Bifidobacterium = new_colors[4],
                   Blautia = new_colors[5],
                   Bacteroides = new_colors[6],
                   Ruminococcus = new_colors[7],
                   Coprococcus = new_colors[8],
                   Collinsella = new_colors[9],
                   Roseburia = new_colors[10],
                   Mediterraneibacter = new_colors[11],
                   Lachnospira = new_colors[12],
                   Anaerobutyricum = new_colors[13],
                   Streptococcus = new_colors[14],
                   Phocaeicola = new_colors[15])

ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = column_colors) + # Set the new color palette
  ylim(0, 2.5) + # Set the y-axis limits
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic


#Cattle

# Load the ggplot2 library
library(ggplot2)

# Read in the data
df <- read.csv("/Users/valter/Downloads/Kraken2/Genera/ca_agrf_genera_copy.csv", header = TRUE)

# Subset the data to only include the first 15 columns
df_subset <- df[, 2:16]

# Convert the data to long format
df_long <- tidyr::gather(df_subset, key = "variable", value = "value")

# Define a new color palette
new_colors <- c("#FF7F00", "#E41A1C", "#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FFC20A", "#F781BF", "#A65628", "#F0E442", "#000000", "#922B21", "#A6CEE3", "#900C3F", "#2B8884")

ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = column_colors) + # Set the new color palette
  ylim(0, 0.3) + # Set the y-axis limits
  theme_minimal() +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  coord_flip() +
  theme(text = element_text(face = "italic")) # Set the font to italic
