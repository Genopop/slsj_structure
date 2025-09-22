import json
import pickle
import geneakit as gen

# Get paths to data files
with open("../paths.json", 'r') as file:
    paths = json.load(file)

# Load the BALSAC genealogy
ped = gen.genealogy(paths['balsac_genealogy'])

# Load the probands who married in Saguenay–Lac-Saint-Jean between 1931 and 1960
with open(paths['wd'] + "results/pickles/balsac_probands.pkl", 'rb') as file:
    pro = pickle.load(file)

# Compute the expected kinship
phi = gen.phi(ped, pro=pro, verbose=True)

# Save the output as a pickle
with open(paths['wd'] + "results/pickles/balsac_kinship.pkl", 'wb') as file:
    pickle.dump(phi, file)