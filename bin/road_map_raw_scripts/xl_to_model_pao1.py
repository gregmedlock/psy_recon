import cobra
import os
import pandas as pd
from cobra import Reaction, Metabolite, Model
import numpy

#generate the list of reactions to be used for generating the model, from the curation file
xl_file = pd.read_excel('/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/SupplementaryData4_mPAO1.xlsx')

#at this point we have a grab bag of reactions ready to be inserted into the model, which we have yet to create
universal = cobra.io.load_json_model('/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/seed_universal.json') #load the grab bag
universal.reactions
for x in universal.reactions:
    x.id = x.id + "0" #this is to match the formating of the rxn id in the PST model
for y in universal.metabolites:
    y.id = y.id + "0" #formatting again to match previous model generated
universal.metabolites
#reaction ids now match those found in pst

PST = Model('PST') #empty model
pst_nan=[]
#add reactions from the grab bag
for i in range(0,len(xl_file)):
    reaction_row = xl_file.iloc[i,]
    reaction_id = reaction_row['Abbreviation'] + '_c0'
 #   print (reaction_id)
    reaction_to_add = universal.reactions.get_by_id(reaction_id).copy()
    reaction_to_add.subsystem = reaction_row ['Subsystems'] #add subsystem for the reaction 
 #   pst_nan.append(reaction_to_add.subsystem)   
    reaction_to_add.gene_reaction_rule = reaction_row['GPR'] #add gpr for the reaction
    PST.add_reaction(reaction_to_add)
    print (PST.reactions.get_by_id(reaction_id))
#print (pst_nan)
#print ([g.id for g in PST.genes])
#print ([m.id for m in PST.metabolites])
gen = [rxn.gene_reaction_rule for rxn in PST.reactions]
#print (PST.reactions[0].gene_reaction_rule)
#write the model
#for x in gen:
 #   if type(x) is not str:
  #      print(x)
#print (gen)
cobra.io.save_json_model(PST,"pao1.json")
#do the reactions have gprs and do the contain wht you expect

