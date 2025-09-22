import json
import pickle

import geneakit as gen
import pandas as pd

# Get paths to data files
with open("../paths.json", 'r') as file:
    paths = json.load(file)

# Associations between CARTaGENE identifiers and BALSAC identifiers
with open(paths['cartagene_balsac_matches'], 'r') as file:
    association_df = pd.read_csv(file, sep=' ')

# A dictionary to convert CARTaGENE identifiers into BALSAC identifiers
cag_to_balsac = dict()
for index, (balsac, cartagene) in association_df.iterrows():
    cag_to_balsac[f'{cartagene}_{cartagene}'] = int(balsac)

# Load the realised kinship matrix
ibd_matrix = pd.read_csv(paths['wd'] + "results/realised_kinship", index_col=0, sep='\t')

# Extract the individuals who are in both genetic (realised)
# and genealogical (expected) sets
codes = [code for code in sorted(ibd_matrix.index) if code in cag_to_balsac.keys()]
ibd_matrix = ibd_matrix.loc[codes, codes]
inds = [cag_to_balsac[code] for code in codes]
ibd_matrix.index = inds
ibd_matrix.columns = inds
ibd_matrix.sort_index(axis=0, inplace=True)
ibd_matrix.sort_index(axis=1, inplace=True)
inds = sorted(inds)

# Load the genealogy
ped = gen.genealogy(paths['balsac_genealogy'])
iso_ped = gen.branching(ped, pro=inds)
pro = gen.pro(iso_ped)

# Compute the expected kinship
phi = gen.phi(ped, pro=pro, verbose=True)

# Save the output
with open(paths['wd'] + "results/pickles/cartagene_expected_kinship.pkl", 'wb') as file:
    pickle.dump(phi, file)