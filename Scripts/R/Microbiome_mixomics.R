
library(mixOmics) # import the mixOmics library

## ==================================
#  Section 1: PCA 
## ==================================

df <- read.csv("/Users/valter/Downloads/R_Microbiome/microbioma_r.csv",
               header = TRUE)

head(df)
ncol(df)
nrow(df)

df

pca

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

#Working

## ==================================
#  Section 3: sparse PLS-DA first pass 
## ==================================
# Before we embark on tuning the optimal number of variables to select per component
# lets first do a first pass

## -----------------------------------------------------------------------------------------------------------
basic.koren.splsda = splsda(df[,2:6], df$Host_label, logratio = 'CLR', ncomp = 5, keepX = c(5, 5, 5)) # select 50 variables per component

plotIndiv(basic.koren.splsda, ind.names = FALSE, legend = TRUE)

#Working#

#____________________________________

#Testing:

## ==================================
#  Section 3: sparse PLS-DA with tuning
## ==================================
# we first set up a grid of values of keepX we want to evaluate
# then, using repeated cross-validation we evaluate the performance for this grid of values
# It may take a few minutes to run, you can skip this part and look on the website what the output
# looks like

## ---- fig.cap = "FIGURE 4: Tuning keepX for the sPLS-DA performed on the Koren OTU data. 
#Each coloured line represents the balanced error rate (y-axis) per component 
#across all tested keepX values (x-axis) with the standard deviation 
#based on the repeated cross-validation folds."----

grid.keepX = c(seq(5,150, 5))
??splsda

??tune.splsda

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








#Renata:
#---------------------------

library(tidyverse)
library(factoextra)
library(FactoMineR)
library(wesanderson)

gpca <- PCA(df[,2:6], graph = FALSE)


pal <- wesanderson::wes_palette("FantasticFox1", 
                                length(unique(df$Host_label)), type = "continuous")[1:length(unique(df$Host_label))]

plot(gpca)

fviz_pca_ind(gpca,
             col.var="contrib",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(df$Host_label), # color by groups
             palette = pal,
             addEllipses = FALSE, # Concentration ellipses
             legend.title = "Groups") + geom_text(aes(label = df$Host_label))

fviz_pca_var(gpca, col.var="contrib")+
  scale_color_gradient2(low="gray", mid="blue",
                        high="red", midpoint=96) +
  theme_minimal()

fviz_pca_biplot(gpca,
                label='var',
                habillage = factor(df$Host_label), # color by groups
                palette = pal,
                #addEllipses = TRUE, # Concentration ellipses
                legend.title = "Groups")

