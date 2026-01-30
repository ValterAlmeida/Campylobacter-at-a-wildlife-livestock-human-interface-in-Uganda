#Completeness and contamination

#All genomes

df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam.csv",
               header = TRUE)
head(df)
ncol(df)
nrow(df)

library(ggplot2)

# read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam.csv", header = TRUE)

# create scatterplot
ggplot(df, aes(x=Completeness, y=Contamination, color=Contamination)) + 
  geom_point(size=3, alpha=0.6) +
  scale_color_gradient(low = "#0072B2", high = "#FF4500") +
  xlab("Completeness (%)") + 
  ylab("Contamination (%)") +
  ggtitle("MAGs - Completeness and Contamination") +
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

#Test to put in the bottom:

library(ggplot2)

# read in data
df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam.csv", header = TRUE)

# create scatterplot

  ggplot(df, aes(x=Completeness, y=Contamination, color=Completeness)) + 
  geom_point(size=3, alpha=0.6) +
  scale_color_gradient(low = "#0072B2", high = "#FF4500") + 
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

  # Test 3
  
  library(ggplot2)
  
  # read in data
  df <- read.csv("/Users/valter/Downloads/R_Microbiome/MAGs/Comple_contam.csv", header = TRUE)
  
  # create scatterplot
  ggplot(df, aes(x=Completeness, y=Contamination, color=Contamination)) + 
    geom_point(size=3, alpha=0.6) +
    scale_color_gradient(low = "#0072B2", high = "#FF4500") +
    xlab("Completeness (%)") + 
    ylab("Contamination (%)") +
    ggtitle("MAGs - Completeness and Contamination") +
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
  
  
  


