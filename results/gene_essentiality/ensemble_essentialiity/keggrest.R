#import kegg rest api
library(KEGGREST)
library(tidyverse)
filevar <- 'val_100_e3'

filesuf <- 'uniquefor_r.csv'
filename <-"val_100_e3uniquefor_r.csv"
filename <- paste(filevar, filesuf, sep = "")
table <-  read.csv(file.choose())

tablen<-table[,2:2]
View(tablen)
filesuf_save <- '_orthologies.csv'
filename <- paste(filevar, filesuf_save, sep = "")
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