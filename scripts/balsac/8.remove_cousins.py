import json
import logging
import pickle
from pathlib import Path
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

paths_file = Path("../paths.json")

with open(paths_file) as f:
    paths = json.load(f)
    
def load_kinship():
    logger.info(f"Loading expected kinship matrix for probands...")
    with open(Path(paths["wd"]) / "results/pickles/balsac_kinship.pkl", "rb") as f:
        return pickle.load(f)

def remove_siblings_and_cousins(phi):
    logger.info(f"Removing siblings and cousins from probands...")
    pro = phi.columns.tolist()
    indices_to_drop = np.where(np.triu(phi.values, k=1) >= 0.125)
    ids_to_drop = set(phi.columns[indices_to_drop[1]])
    filtered_probands = [pid for pid in pro if pid not in ids_to_drop]
    filtered_kinship = phi.loc[filtered_probands, filtered_probands]
    logger.info(f"Retained {len(filtered_probands)} probands after removing siblings and cousins.")
    return filtered_probands, filtered_kinship

def main():
    phi = load_kinship()
    filtered_probands, filtered_kinship = remove_siblings_and_cousins(phi)
    with open(Path(paths["wd"]) / "results/pickles/balsac_probands_nocousins.pkl", "wb") as f:
        pickle.dump(filtered_probands, f)
    with open(Path(paths["wd"]) / "results/pickles/balsac_kinship_nocousins.pkl", "wb") as f:
        pickle.dump(filtered_kinship, f)

if __name__ == "__main__":
    main()