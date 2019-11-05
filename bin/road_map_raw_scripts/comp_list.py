#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 16:30:40 2018

@author: pjt7wd
"""

import os
import pandas as pd
import cobra
from cobra import Reaction, Metabolite, Gene

#import the excel file for PST_curation
xl_file = pd.read_excel('/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/PST_curation.xlsx')

#import the reaction list from the draft
rxn_list = pd.read_csv('/Users/pjt7wd/Desktop/PST_reaction_draft.csv')

#create new list blank
not_in_curation = []
just_ids = []
draft_rxns = []
#compare the dataframes

for x in range(0, len(xl_file)):
    reaction_row = xl_file.iloc[x,]
    n_id = reaction_row['Reaction ID']
    just_ids.append(n_id)
print (just_ids)

for x in range(0, len(rxn_list)):
    reaction_row_m = rxn_list.iloc[x,]
    m_id = reaction_row_m['id']
    draft_rxns.append(m_id)
print (draft_rxns)

for x in draft_rxns:
    if x not in just_ids:
        not_in_curation.append(x)
    else:
        pass
print (not_in_curation)

#df = pd.DataFrame(not_in_curation)
#print (df)
#df.to_csv("/Users/pjt7wd/Desktop/PST_updates_list.csv")