#import kegg rest api
library(KEGGREST)
library(tidyverse)
filevar <- 'glucose_mm_100_for_analysis'

filesuf <- 'uniquefor_r.csv'
filename <-"../glucose_mm_100_for_analysisuniquefor_r.csv"
filename <- paste(filevar, filesuf, sep = "")


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

filevar <- 'gaba_mm_100_for_analysis'
filesuf <- 'uniquefor_r.csv'
filename <-"gaba_mm_100_for_analysisuniquefor_r.csv"
filename <- paste(filevar, filesuf, sep = "")
table <-  read.csv("gaba_mm_100_for_analysisuniquefor_r.csv")

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

filevar <- 'asp_mm_100_for_analysis'
filesuf <- 'uniquefor_r.csv'
filename <-"asp_mm_100_for_analysisuniquefor_r.csv"
filename <- paste(filevar, filesuf, sep = "")
table <-  read.csv("asp_mm_100_for_analysisuniquefor_r.csv")

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

filevar <- 'cit_mm_100_for_analysis'
filesuf <- 'uniquefor_r.csv'
filename <-"cit_mm_100_for_analysisuniquefor_r.csv"
filename <- paste(filevar, filesuf, sep = "")
table <-  read.csv("cit_mm_100_for_analysisuniquefor_r.csv")

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

filevar <- 'suc_mm_100_for_analysis'
filesuf <- 'uniquefor_r.csv'
filename <-"suc_mm_100_for_analysisuniquefor_r.csv"
filename <- paste(filevar, filesuf, sep = "")
table <-  read.csv("suc_mm_100_for_analysisuniquefor_r.csv")

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