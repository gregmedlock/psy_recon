#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 10:09:51 2018

@author: pjt7wd
"""

import os
import pandas as pd
import cobra
from cobra import Reaction, Metabolite, Gene

#read the two library files for the lists, as well as the reaction list to call reactions
supp_pao1 = pd.read_excel("/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/SupplementaryData4_mPAO1.xlsx")
orthologs = pd.read_csv('/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/PA01_PSY_orthologs.csv')
rxn_list = pd.read_csv('/Users/pjt7wd/Desktop/PST_updates_list.csv')
print (rxn_list)

#create reaction list to call in the supp
rxns = []
for x in range(0,len(rxn_list)):
    reaction_row = rxn_list.iloc[x,]
    ids = reaction_row['0']
    rxns.append(ids)
print (rxns)

#now we have the reactions extracted into a list to compare with the pa01 library, which must be similarly compared
#when comparing the reaction, create a dictionary with the reaction id as the key and the genes as the value
pao1_genes = {}

#extract the reactions from pao1 and the corresponding genes to create a dictionary
for x in range(0,len(supp_pao1)):
    reaction_row = supp_pao1.iloc[x,]
    if "rxn" in reaction_row["Abbreviation"]:
        rxn_abb = reaction_row['Abbreviation'] + "_c0"
        rxn_genes = reaction_row['Genes']
    print (rxn_genes)
    print (rxn_abb)
    #populate the dictionary with the keys (rxn abbreviations) and values (gene ids)
    pao1_genes[rxn_abb] = rxn_genes
    print (pao1_genes)    

#now compare the dictionary keys (rxns) to the list of rxns in the PST update list of rxns, rxn_list
matching_rxns = {}
for k in rxns:
    if k in pao1_genes:
        print (k, pao1_genes[k])
        str_pao1_genes = str(pao1_genes[k])
        split_genes = str.split(str_pao1_genes, ' ')
        print (split_genes)
        matching_rxns[k] = split_genes        
        print (matching_rxns)

#this has the format that I want, but now add to a pandas dataframe
df_matching_rxns = pd.DataFrame.from_dict(matching_rxns, orient='index')
print (df_matching_rxns)

#dataframe is now in the format to be used. it includes the reaction ids and the associated genes
#now I must match the asociated genes with the orthologous pairs
#first create dictionary with the orthologous pairs
ortholog_dict = {}
print (orthologs)
for x in range(0, len(orthologs)):
    pair_row = orthologs.iloc[x,]
    query = pair_row['Locus Tag (Query)']
    result = pair_row['Locus Tag (Hit)']
    ortholog_dict[query]=result
    print (ortholog_dict)

#the ortholog dict is assembled, with each query value (all Pao1 genes and PSPTO genes were surveyed)
#this dict will now be iterated over to match with the missing reactions PAo1 genes
rxn_PST = {}
for k,v in matching_rxns.items():
    for x in matching_rxns[v]:
        if x in ortholog_dict:
            print (k ,v, ortholog_dict[x])
            rxn_PST[k] = ortholog_dict[x]
            print (rxn_PST)
        else:
            pass
df_match_orthologs = pd.DataFrame.from_dict(rxn_PST, orient='index') 
print(df_match_orthologs)   
df_match_orthologs.to_csv("/Users/pjt7wd/Desktop/PST_gpr.csv")
#the above code will only produce matches for all singly annotated reactions with one gene annoted
#plural genes have not been matched this way
for k in matching_rxns.items():
    if k in ortholog_dict:
        print (k, v, ortholog_dict[k])
    else:
        print (k ,v, 'no')
        pass
     #   df_matching_rxns[]
        