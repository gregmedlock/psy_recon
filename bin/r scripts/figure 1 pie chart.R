##produce a gene enrichment pie chart with the following data structure in a .csv file
## (Pathway, Genes, Type) where Pathway is the ontilogical pathway from kegg, Genes is the number of genes with that ontology, and Type being the ontological catgory (subjective)

#choose the file to upload, use .csv format
df <- read.csv(file(file.choose()), header=TRUE, sep = ",")