library(mixOmics) # import the mixOmics library

## ==================================
#  Section 1: PCA 
## ==================================

df <- read.csv("/Users/valter/Downloads/R_Microbiome/microbiomes_orion.csv",
               header = TRUE)

library(ggplot2)
library(tidyverse)

plot(df[,2:6])

ggplot(data = df) +

  geom_point(mapping = aes(x = Host_label, y = Firmicutes)) +
  geom_point(mapping = aes(x = Host_label, y = Actinobacteria))



# Microbiome - Youtube

library(tidyverse)
library(readxl)

metadata <- read.csv("/Users/valter/Downloads/R_Microbiome/microbiomes_orion.csv", na="NA") %>%
  select(Host_label, Firmicutes) %>%
  drop_na(Host_label)

otu_counts <- read_t("raw_data/schubert.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")

taxonomy <- read_tsv("raw_data/schubert.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")

otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon")



# MixOmics:
head(df)
ncol(df)
nrow(df)

microbiome.pca = pca(df[,2:6], ncomp = 5, logratio = 'CLR') # undergo PCA with 10 components
plot(microbiome.pca) # plot explained variance

## ---- fig.cap = "FIGURE 2: Sample plots from PCA performed on the Koren OTU data. Samples are projected into the space spanned by the first two components."----
df$Origin
df$Host
# cbdind(df$Host_label, df$Origin) 
plotIndiv(microbiome.pca, # plot samples projected onto PCs
          ind.names = FALSE, # not showing sample names
          group = df$Host_label, # color according to Y 
          legend = TRUE,
          title = 'PCA Components')


## Question: now, try without a logratio transformation, what do you observe?

## ==================================
#  Section 2: PLS-DA first pass 
## ==================================
# here we evaluate the performance of the PLSDA with repeated cross validation,
# and choose the number of components for the remainder of the analysis
# we choose a large-ish number of component
## -----------------------------------------------------------------------------------------------------------
basic.koren.plsda = plsda(df[,2:6], df$Host_label, logratio = 'CLR', 
                          ncomp = 5)

plotIndiv(basic.koren.plsda, ind.names = FALSE, legend = TRUE)

## ==================================
#  Section 3: sparse PLS-DA first pass 
## ==================================
# Before we embark on tuning the optimal number of variables to select per component
# lets first do a first pass

## -----------------------------------------------------------------------------------------------------------
basic.koren.splsda = splsda(df[,2:6], df$Host_label, logratio = 'CLR', ncomp = 5, keepX = c(5, 5, 5)) # select 50 variables per component

plotIndiv(basic.koren.splsda, ind.names = FALSE, legend = TRUE)

#Working

#Testing:

#Testing:

## ==================================
#  Section 3: sparse PLS-DA with tuning
## ==================================
# we first set up a grid of values of keepX we want to evaluate
# then, using repeated cross-validation we evaluate the performance for this grid of values
# It may take a few minutes to run, you can skip this part and look on the website what the output
# looks like

## ---- fig.cap = "FIGURE 4: Tuning keepX for the sPLS-DA performed on the Koren OTU data. Each coloured line represents the balanced error rate (y-axis) per component across all tested keepX values (x-axis) with the standard deviation based on the repeated cross-validation folds."----
grid.keepX = c(seq(5,150, 5))

koren.tune.splsda = tune.splsda(df[,2:6], df$Host_label,
                                ncomp = optimal.ncomp, # use optimal component number
                                logratio = 'CLR', # transform data to euclidean space
                                test.keepX = grid.keepX,
                                validation = c('Mfold'),
                                folds = 5, nrepeat = 10, # use repeated CV
                                dist = 'max.dist', # maximum distance as metric
                                progressBar = FALSE)

# results are plotted here
plot(koren.tune.splsda) # plot this tuning




