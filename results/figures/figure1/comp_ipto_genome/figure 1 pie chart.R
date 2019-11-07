##produce a gene enrichment pie chart with the following data structure in a .csv file
## (Pathway, Genes, Type) where Pathway is the ontilogical pathway from kegg, Genes is the number of genes with that ontology, and Type being the ontological catgory (subjective)

#choose the file to upload, use .csv format
df <- read.csv(file(file.choose()), header=TRUE, sep = ",")

#load in the packages needed
library(tidyverse)
library(ggplot2)
library(ggrepel)
#create a factor for the Type object
dfType<-factor(df$Type)

#create sums for the genes in the Type factor using aggragate
sum_data<-aggregate(df$Genes, list(df$Type), sum)

#change teh column names to prevent weirdness with x call
colnames(sum_data) <- c('Category', 'Genes')
sum_data$id <- paste(sum_data$Category, sum_data$Genes, sep=" ")


ggplot(data= sum_data, aes(x = "", y = Genes, fill =  reorder(id, -Genes))) + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  #geom_text(aes(label = paste0(Genes, " ", Category), x = 1.63), position = position_stack(vjust = 0.5)) + 
  coord_polar(theta = "y") +
  scale_fill_discrete(name = "Ontological Groups") +
  ggtitle("Kegg Pto Genome Ontological Groupings") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom',
)
