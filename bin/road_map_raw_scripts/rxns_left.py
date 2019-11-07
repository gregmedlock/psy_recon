#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 09:32:28 2018

@author: pjt7wd
"""

#compare the list of annotated vs unannotated pao1 matches to generate new list of unnamed rxns
import pandas as pd

unannotated_rxns = pd.read_excel('/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/reactions_left.xlsx')
matched_rxns = pd.read_csv("/Users/pjt7wd/Desktop/Psy_modeling/psy_recon/data/matched_rxns.xlsx")

print(matched_rxns)
