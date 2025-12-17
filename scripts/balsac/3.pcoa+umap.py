import json
import numpy as np
import pickle
from skbio.stats.ordination import pcoa
import umap

# Get paths to data files
with open("../paths.json", 'r') as file:
    paths = json.load(file)

# Load the expected kinship of BALSAC probands
with open(paths['wd'] + "results/pickles/balsac_kinship.pkl", 'rb') as file:
    phi = pickle.load(file)

# Remove individuals who have kinship with less than half of other individuals
counts = phi.astype(bool).sum(axis=1)
inds_to_keep = sorted(counts[counts > counts.shape[0] / 2].index)
phi = phi.loc[inds_to_keep, inds_to_keep]

# Remove siblings while keeping one individual per family
indices_to_drop = np.where(np.triu(phi.values, k=1) >= 0.25)
ids_to_drop = set(phi.columns[indices_to_drop[1]])
print(f"Removed {len(ids_to_drop)} siblings.")

# Remove the siblings from the kinship matrix
pruned = phi.drop(list(ids_to_drop), axis=0).drop(list(ids_to_drop), axis=1)

# Set the diagonal to 1 for MDS and UMAP
np.fill_diagonal(pruned.values, 1)

# Save the IDs of the remaining individuals
with open(paths['wd'] + "results/pickles/balsac_pro_0.5_nosibs.pkl", 'wb') as file:
    pickle.dump(list(pruned.columns), file)

# Compute PCoA for UMAP initialization
pcoa_results = pcoa(1 - pruned, number_of_dimensions=1000, method='eigh', seed=0, warn_neg_eigval=False, inplace=True)
init_embedding = pcoa_results.samples.values[:, :2]

# Save the PCoA output
with open(paths['wd'] + "results/pickles/balsac_pcoa.pkl", 'wb') as file:
    pickle.dump(pcoa_results, file)

# Compute the UMAP projection
emb = umap.UMAP(
    metric='precomputed',
    n_neighbors=15,
    min_dist=0.5,
    verbose=True,
    random_state=0,
    init=init_embedding
).fit_transform(1 - pruned)

# Save the UMAP output
with open(paths['wd'] + "results/pickles/balsac_umap.pkl", 'wb') as file:
    pickle.dump(emb, file)

# Compute the UMAP projection from the PCoA
emb = umap.UMAP(
    n_neighbors=15,
    min_dist=0.5,
    verbose=True,
    random_state=0,
    init=init_embedding
).fit_transform(pcoa_results.samples.values)

# Save the UMAP output
with open(paths['wd'] + "results/pickles/balsac_umap_from_pcoa.pkl", 'wb') as file:
    pickle.dump(emb, file)
