#import kegg rest api
library(KEGGREST)
library(tidyverse)

filevar <- 'complete_100_e3'

filepath <- '../ensemble_essentialiity/'
filesuf <- 'uniquefor_r.csv'

filename <- paste(filepath,filevar, filesuf, sep = "")
table <-  read.csv(filename)

tablen<-table[,2:2]
View(tablen)
filesuf_save <- '_orthologies.csv'
filename <- paste(filepath + filevar + filesuf_save, sep = "")
for (val in tablen) {
  n <- keggLink("pathway", val)
  for (paths in n){
   q <- keggList(paths)
   #print (as.character(val))
   #print (as.character(q))
   m = data.frame(as.character(val),as.character(q))
   write.table(m, file=filename,
               append=TRUE,
               col.names = FALSE,
               sep = ',')
  }
}
print ("done")


       