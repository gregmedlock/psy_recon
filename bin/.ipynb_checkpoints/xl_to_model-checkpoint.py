import cobra
import os
import pandas as pd
from cobra import Reaction, Metabolite, Model

#generate the list of reactions to be used for generating the model, from the curation file
xl_file = pd.read_excel('../data/PST_feeder.xlsx')

#at this point we have a grab bag of reactions ready to be inserted into the model, which we have yet to create
universal = cobra.io.load_json_model('../data/seed_universal.json') #load the grab bag
universal.reactions
for x in universal.reactions:
    x.id = x.id + "0" #this is to match the formating of the rxn id in the PST model
for y in universal.metabolites:
    y.id = y.id + "0" #formatting again to match previous model generated
universal.metabolites
#reaction ids now match those found in pst

PST = Model('PST') #empty model

#add reactions from the grab bag
for i in range(0,len(xl_file)):
    reaction_row = xl_file.iloc[i,]
    reaction_id = reaction_row['Reaction ID']
    reaction_to_add = universal.reactions.get_by_id(reaction_id).copy()
    reaction_to_add.subsystem = reaction_row ['Subsystems'] #add subsystem for the reaction 
    reaction_to_add.gene_reaction_rule = reaction_row['GPR'] #add gpr for the reaction
    PST.add_reaction(reaction_to_add)

#write the model
cobra.io.save_json_model(PST,"../results/pst_feeder.json")

