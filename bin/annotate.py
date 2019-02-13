import pandas as pd
import cobra

features = pd.read_csv('../data/pst_kegg.tsv', sep = '^')
features['Kegg Ontology'], features['RefSeq'] = features['KO'].str.split('|').str

# Separate the KO number and text description into separate columns.
features['KO'] = features['KO'].str[:6]
features['Kegg Ontology'] = features['Kegg Ontology'].str[6:]

# Replace missing annotations with consistent nomenclature
features['KO'] = features['KO'].str.replace('no KO ','none')
features.loc[features['KO'] == 'none','Kegg Ontology'] = 'none'

# Remove the RefSeq tag in the RefSeq column
features['RefSeq'] = features['RefSeq'].str[10:]

# extract EC numbers and then transfer to reactions
gene_to_EC = {}
for geneid in features['kegg id']:
    if features.loc[features['kegg id'] == geneid, 'Kegg Ontology'].str.find('EC:').values[0] > -1:
        # Extract the EC number from the Kegg Ontology column
        gene_to_EC[geneid] = features.loc[features['kegg id'] == geneid, 'Kegg Ontology'].str.split('EC:').values[0][1].split(']')[0]
    else:
        gene_to_EC[geneid] = 'none'


gene_to_EC = pd.Series(gene_to_EC, name = 'EC')
gene_to_EC.index.name = 'kegg id'
gene_to_EC = gene_to_EC.reset_index()

features = features.merge(gene_to_EC, on='kegg id')

#import reconstruction to be annotated
psy = cobra.io.read_sbml_model("../results/reconstructions/suffix_fixed_no_biomass.xml")

# fix seed reaction/metabolite identifiers by removing '0' suffix if it exists
for metabolite in psy.metabolites:
    if metabolite.id.endswith('0'):
        metabolite.id = metabolite.id[:-1]

for reaction in psy.reactions:
    if reaction.id.endswith('0'):
        reaction.id = reaction.id[:-1]

#For each gene that is in the model, add the annotated kegg gene and refseq function name
for index, row in features.iterrows():
    gene_id = row['kegg id']
    if gene_id in psy.genes:
        gene_obj = psy.genes.get_by_id(gene_id)
        gene_obj.annotation['kegg.genes'] = 'pst:' + gene_id

        # Assign the function as the name for the gene
        gene_obj.name = row['RefSeq'].split(';')[0]

        # Add the EC annotation to the reactions for each gene
        for reaction in gene_obj.reactions:
            reaction_obj = psy.reactions.get_by_id(reaction.id)
            reaction_obj.annotation['ec-code'] = features.loc[features['kegg id'] == gene_id, 'EC']

# Next, load the modelseed reaction and compound aliases to annotate reactions and metabolites.
seed_rxn_aliases = pd.read_csv('../data/modelseed_data/Reactions_Aliases.tsv', sep = '\t')
seed_cpd_aliases = pd.read_csv('../data/modelseed_data/Compounds_Aliases.tsv', sep = '\t')

# Replace the source IDs to be consistent with the identifiers memote is looking for
# These are the MIRIAM compliant versions of the resources, available at identifiers.org
# see the memote annotations.py file for the regular expressions expected for
# identifiers from each resource.
seed_rxn_aliases.loc[seed_rxn_aliases['Source'] == 'KEGG', 'Source'] = 'kegg.reaction'
seed_rxn_aliases.loc[seed_rxn_aliases['Source'] == 'BiGG', 'Source'] = 'bigg.reaction'
seed_rxn_aliases.loc[seed_rxn_aliases['Source'] == 'MetaCyc', 'Source'] = 'biocyc'

seed_cpd_aliases.loc[seed_cpd_aliases['Source'] == 'KEGG', 'Source'] = 'kegg.compound'
seed_cpd_aliases.loc[seed_cpd_aliases['Source'] == 'BiGG', 'Source'] = 'bigg.compound'
seed_cpd_aliases.loc[seed_cpd_aliases['Source'] == 'MetaCyc', 'Source'] = 'biocyc'


# Get and add all of the reaction annotations
for reaction in psy.reactions:
    if not reaction.id.startswith('EX_'): # don't look for exchange reactions
        # get the reaction ID without compartment suffix e.g. '_c'
        reaction_baseid = reaction.id.split('_')[0]

        annotation_dict = {}
        annotation = seed_rxn_aliases.loc[seed_rxn_aliases['MS ID'] == reaction_baseid]
        for source in annotation['Source'].unique():
            db_annotation = annotation.loc[annotation['Source'] == source]
            annotation_id = db_annotation['External ID'].values[0]
            annotation_dict[source] = annotation_id

        # if there were no annotations, this reaction either has no alias
        # or is only in modelseed.
        if (reaction.id == (reaction_baseid + '_c')) or (reaction.id == (reaction_baseid + '_e')):
                # if the ID is a standard modelseed format, add a seed annotation.
                # otherwise, this might be a custom object (which should not have an annotation)
                annotation_dict['seed.reaction'] = reaction_baseid

        # add each annotation in a way that maintains any existing annotations
        for key in annotation_dict:
            reaction.annotation[key] = annotation_dict[key]


# Get and add all of the metabolite annotations
for metabolite in psy.metabolites:
    # get the reaction ID without compartment suffix e.g. '_c'
    metabolite_baseid = metabolite.id.split('_')[0]

    annotation_dict = {}
    annotation = seed_cpd_aliases.loc[seed_cpd_aliases['MS ID'] == metabolite_baseid]
    for source in annotation['Source'].unique():
        db_annotation = annotation.loc[annotation['Source'] == source]
        annotation_id = db_annotation['External ID'].values[0]
        annotation_dict[source] = annotation_id

    if (metabolite.id == (metabolite_baseid + '_c')) or (metabolite.id == (metabolite_baseid + '_e')):
                # if the ID is a standard modelseed format, add a seed annotation.
                # otherwise, this might be a custom object (which should not have an annotation)
        annotation_dict['seed.compound'] = metabolite_baseid

    # add each annotation in a way that maintains any existing annotations
    for key in annotation_dict:
        metabolite.annotation[key] = annotation_dict[key]


# Add inchi keys for all metabolites using the ModelSEED biochemistry files.
seed_cpd_structures = pd.read_csv('../data/modelseed_data/ModelSEED_Structures.txt', sep = '\t')

# Add the inchi and inchikey annotations for metabolites
for metabolite in psy.metabolites:
    metabolite_baseid = metabolite.id.split('_')[0]
    if metabolite_baseid in seed_cpd_structures['ID'].tolist():
        annotations = seed_cpd_structures.loc[seed_cpd_structures['ID'] == metabolite_baseid]
        if 'InChI' in annotations['Type'].tolist():
            inchi = seed_cpd_structures.loc[(seed_cpd_structures['ID'] == metabolite_baseid) &
                                        (seed_cpd_structures['Type'] == 'InChI'),'Structure'].values[0]
            metabolite.annotation['inchi'] = inchi

        if 'InChIKey' in annotations['Type'].tolist():
            inchikey = seed_cpd_structures.loc[(seed_cpd_structures['ID'] == metabolite_baseid) &
                                        (seed_cpd_structures['Type'] == 'InChIKey'),'Structure'].values[0]
            metabolite.annotation['inchikey'] = inchikey

# Add SBO terms

# The expected SBO terms are as follows:
# Metabolite SBO:0000247
# Metabolic Reaction SBO:0000176
# Transport Reaction SBO:0000185
# Exchange Reaction SBO:0000627
# Demand Reaction SBO:0000628
# Sink Reactions SBO:0000632
# Gene SBO:0000243
# Biomass Reactions SBO:0000629

# For all metabolites, add the metabolite annotation
for metabolite in psy.metabolites:
    metabolite.annotation['sbo'] = 'SBO:0000247'

# for all reactions, add the reaction annotation.
# where applicable, add biomass, transport, and exchange
# terms as well.
for reaction in psy.reactions:
    annotations = []
    if reaction.id.startswith('EX_'):
        annotations.append('SBO:0000627')
    elif len(set([met.compartment for met in reaction.metabolites])) > 1:
        annotations.append('SBO:0000185')
    elif reaction.id.lower().find('bio') > -1:
        annotations.append('SBO:0000629')
    else:
        annotations.append('SBO:0000176')

    if len(annotations) > 1:
        reaction.annotation['sbo'] = annotations
        print(reaction.id + ' has more than one SBO annotation. Memote will not like this.')
    else:
        # save as a single element to cope with the memote issues with lists of
        # sbo annotations.
        reaction.annotation['sbo'] = annotations[0]

    # add gene annotations
for gene in psy.genes:
    gene.annotation['sbo'] = 'SBO:0000243'

# save the model with added annotations

# open the exchange reactions before saving
for reaction in psy.reactions:
    if reaction.id.startswith('EX'):
        reaction.lower_bound = -1000
cobra.io.write_sbml_model(psy,'../results/reconstructions/v4_with_all_annotations.xml')
