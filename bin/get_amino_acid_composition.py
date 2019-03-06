import collections
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json

all_aas = collections.defaultdict(int)
for rec in SeqIO.parse("../data/psy_dc3000_ref_genome.fasta", "fasta"):
    x = ProteinAnalysis(str(rec.seq))
    for aa, count in x.count_amino_acids().items():
        all_aas[aa] += count

# Normalize to total count
total_aas = sum(all_aas.values())
normalized_aas = {key:float(value)/float(total_aas) for
                        key,value in all_aas.items()}

# convert to modelseed identifiers
letter_to_seed = {'G': 'cpd00033_c',
'P': 'cpd00129_c',
'A': 'cpd00035_c',
'V': 'cpd00156_c',
'L': 'cpd00107_c',
'I': 'cpd00322_c',
'M': 'cpd00060_c',
'C': 'cpd00084_c',
'F': 'cpd00066_c',
'Y': 'cpd00069_c',
'W': 'cpd00065_c',
'H': 'cpd00119_c',
'K': 'cpd00039_c',
'R': 'cpd00051_c',
'Q': 'cpd00053_c',
'N': 'cpd00132_c',
'E': 'cpd00023_c',
'D': 'cpd00041_c',
'S': 'cpd00054_c',
'T': 'cpd00161_c'}

normalized_aas_seed = {letter_to_seed[letter]:value for
                letter,value in normalized_aas.items()}

# save the dictionary
with open('../results/biomass/amino_acid_composition.json', 'w') as fp:
    json.dump(normalized_aas_seed, fp)
