#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 16:07:54 2018

@author: pjt7wd
"""

import os
import pandas as pd
import cobra
from cobra import Reaction, Metabolite, Gene

#import the draft model
PST = cobra.io.read_sbml_model('/Users/pjt7wd/Desktop/Psy_modeling/models/01_31_18/PSY_02.sbml')

#generate initial empty list
PST_reactions = []
#generate the reaction list to be added to teh pandas df
for x in PST.reactions:
    PST_reactions.append(x.id)
print (PST_reactions)

#create dataframe
df = pd.DataFrame(PST_reactions)
print (df)

#convert df to excel
df.to_csv("/Users/pjt7wd/Desktop/PST_reaction_draft.csv")
