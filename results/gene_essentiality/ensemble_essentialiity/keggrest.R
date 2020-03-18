#import kegg rest api
library(KEGGREST)
library(tidyverse)
filevar <- 'suc_mm_100_for_analysis'

filesuf <- 'uniquefor_r.csv'
filename <-"gl_mixed.csv"
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


#choose file and save file name
library(KEGGREST)
library(tidyverse)

filename <-"gl_mixed.csv"
table <-  read.csv(file.choose())

tablen<-table[,2:2]
View(tablen)

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