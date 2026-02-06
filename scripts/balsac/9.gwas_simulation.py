#!/usr/bin/env python3
"""
GWAS Stratification Simulation with Genealogical Pedigrees
===========================================================

Supplementary code for:
"Fine-scale structure of a whole regional population through genetics and genealogies"

Description
-----------
This script demonstrates how cryptic population stratification, driven by specific
founder effects (e.g., Charlevoix founders in the SLSJ region), inflates GWAS 
test statistics. It compares standard correction methods (Genetic PCA) against 
Genealogical PCoA.

Methodology:
1.  **Pedigree Simulation**: Uses `msprime` FixedPedigree model to simulate ancestry 
    through the real BALSAC genealogy.
2.  **Recapitation**: Uses `stdpopsim` (OutOfAfrica_3G09) to provide realistic 
    founder diversity at the top of the pedigree.
3.  **Phenotype**: Simulates a trait correlated with genealogical structure
    (Charlevoix genetic contribution).
4.  **GWAS**: Performs association testing using linear regression.

Usage
-----
    python 9.gwas_simulation.py --config config.json
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set

import matplotlib.pyplot as plt
import msprime
import numpy as np
import pandas as pd
import seaborn as sns
import stdpopsim
import tskit
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.feature_selection import f_regression

# Check for critical dependencies
try:
    import allel
    import geneakit as gk
    from skbio import DistanceMatrix
    from skbio.stats.ordination import pcoa
except ImportError as e:
    sys.exit(f"Critical dependency missing: {e}. Please install requirements.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class SimulationConfig:
    """Simulation parameters and file paths."""
    # Paths
    paths_file: Path = Path("../paths.json")
    output_dir: Path = Path("results")
    
    # Biological parameters (human genome average)
    sequence_length: float = 1e8           # 100 Mb
    mutation_rate: float = 1.2e-8          # Kong et al. 2012, https://doi.org/10.1038/nature11396
    recombination_rate: float = 1.1e-8     # Kong et al. 2002, https://doi.org/10.1038/ng917
    
    # Study design
    sample_size: int = 10000
    seed: int = 2026
    
    # Quality control
    maf_threshold_gwas: float = 0.01
    maf_threshold_pca: float = 0.05
    
    # LD pruning (scikit-allel)
    ld_window: int = 50
    ld_step: int = 5
    ld_threshold: float = 0.2
    
    # Correction & phenotype
    n_pcs1: int = 10
    n_pcs2: int = 100
    h2: float = 0.1                         # Effect size of the stratification

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "SimulationConfig":
        """Initialize config from file and CLI arguments."""
        conf = cls()
        if args.config and Path(args.config).exists():
            with open(args.config) as f:
                data = json.load(f)
            for k, v in data.items():
                if hasattr(conf, k):
                    setattr(conf, k, Path(v) if "path" in k else v)
        conf.seed = args.seed
        return conf


def get_geographic_mappings() -> Dict[int, str]:
    """Returns mapping from integer region codes to region names."""
    return {
        25705: "Abitibi", 25706: "Bas-Saint-Laurent", 25707: "Beauce",
        25708: "Bois-Francs", 25709: "Charlevoix", 25710: "Côte-de-Beaupré",
        25711: "Côte-du-Sud", 25712: "Côte-Nord", 25713: "Estrie",
        25714: "Gaspésie", 25715: "Île-de-Montréal", 25716: "Îles-de-la-Madeleine",
        25717: "Lanaudière", 25718: "Laurentides", 25719: "Mauricie",
        25720: "Outaouais", 25721: "Agglomération de Québec", 25722: "Région de Québec",
        25723: "Nord du Québec", 25724: "Richelieu", 25725: "Rive Nord-Ouest de Montréal",
        25726: "Rive Sud de Montréal", 25727: "Saguenay–Lac-Saint-Jean",
        25728: "Témiscamingue", 27118: "Côte-de-Beaupré", 27119: "Portneuf",
        27120: "Lévis-Lotbinière"
    }


# =============================================================================
# Data Management
# =============================================================================

class DataLoader:
    """Handles robust loading of external genealogical and demographic data."""
    
    def __init__(self, paths_file: Path):
        with open(paths_file) as f:
            self.paths = json.load(f)
        self.region_map = get_geographic_mappings()

    def load_genealogy(self) -> Any:
        path = self.paths.get("balsac_genealogy")
        logger.info(f"Loading genealogy from {path}...")
        return gk.genealogy(path)

    def load_probands(self) -> List[int]:
        path = Path(self.paths["wd"]) / "results/pickles/balsac_probands_nocousins.pkl"
        logger.info(f"Loading probands from {path}...")
        with open(path, "rb") as f:
            probands = pickle.load(f)
            return [int(x) for x in probands]

    def load_region_map_full(self) -> Dict[int, str]:
        """
        Loads region codes for ALL individuals.
        Returns: Dict {individual_id: region_name}
        """
        path = self.paths["demography_information"]
        logger.info(f"Loading full regional demographics from {path}...")
        
        # Load ID (col 0) and parents' region (col 17)
        df = pd.read_csv(
            path, sep=";", encoding="cp1252", na_values="NA", 
            usecols=[0, 17], low_memory=False
        )
        df.columns = ["ind_id", "region_code"]
        
        # QC: Ensure IDs are integers
        df.dropna(subset=["ind_id"], inplace=True)
        df["ind_id"] = df["ind_id"].astype(int)
        
        # QC: Ensure Regions are integers
        df["region_code"] = (
            pd.to_numeric(df["region_code"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        
        # Map to names
        df["region_name"] = df["region_code"].map(self.region_map).fillna("Unknown")
        
        return dict(zip(df["ind_id"], df["region_name"]))


# =============================================================================
# Simulation Engine
# =============================================================================

class GenealogySimulator:
    """Manages the hybrid msprime (FixedPedigree) + stdpopsim (recapitation) pipeline."""

    def __init__(self, ped: Any, config: SimulationConfig):
        self.ped = ped
        self.cfg = config

    def build_pedigree(self, sample_ids: List[int]) -> Tuple[msprime.Pedigree, Dict[int, int], List[int]]:
        """
        Converts BALSAC genealogy to msprime format.
        
        Returns:
            - msprime.Pedigree object
            - Dictionary mapping BALSAC ID -> msprime ID
            - List of ALL involved individuals sorted topologically (ancestors -> descendants)
        """
        logger.info(f"Tracing ancestry for {len(sample_ids)} probands...")
        
        # 1. Identify all ancestors
        ancestors = set()
        stack = list(sample_ids)
        while stack:
            ind_id = stack.pop()
            if ind_id == 0 or ind_id in ancestors:
                continue
            ancestors.add(ind_id)
            
            try:
                ind = self.ped[ind_id]
                if ind.father.ind != 0: stack.append(ind.father.ind)
                if ind.mother.ind != 0: stack.append(ind.mother.ind)
            except KeyError:
                continue

        # 2. Compute topological depth (generations)
        # Generations are assigned so that all parents are strictly older than offspring,
        # satisfying msprime FixedPedigree temporal constraints; absolute times are arbitrary.
        # We start everyone at 0 and bubble up.
        generations = {pid: 0 for pid in ancestors}
        
        changed = True
        all_involved = list(ancestors) # Note: ancestors set includes probands in this logic
        
        while changed:
            changed = False
            for ind_id in all_involved:
                try:
                    ind = self.ped[ind_id]
                except KeyError:
                    continue
                
                t_child = generations[ind_id]
                
                for parent_id in [ind.father.ind, ind.mother.ind]:
                    if parent_id != 0 and parent_id in generations:
                        if generations[parent_id] <= t_child:
                            generations[parent_id] = t_child + 1
                            changed = True
        
        # 3. Build msprime pedigree
        # Sorting High -> Low ensures parents are added before children
        sorted_ids = sorted(all_involved, key=lambda x: generations[x], reverse=True)
        
        pb = msprime.PedigreeBuilder()
        balsac_to_msprime = {}
        
        for ind_id in sorted_ids:
            try:
                ind = self.ped[ind_id]
            except KeyError:
                continue
            
            p1 = balsac_to_msprime.get(ind.father.ind, -1)
            p2 = balsac_to_msprime.get(ind.mother.ind, -1)
            
            mid = pb.add_individual(
                time=generations[ind_id],
                parents=[p1, p2],
                is_sample=(ind_id in sample_ids)
            )
            balsac_to_msprime[ind_id] = mid

        pedigree = pb.finalise()
        pedigree.sequence_length = self.cfg.sequence_length
        
        return pedigree, balsac_to_msprime, sorted_ids

    def run(self, pedigree: msprime.Pedigree) -> tskit.TreeSequence:
        """Executes the simulation pipeline."""
        
        # 1. Forward Time Simulation
        logger.info(f"Simulating fixed pedigree ({len(pedigree.individuals)} individuals)...")
        ts_ped = msprime.sim_ancestry(
            initial_state=pedigree,
            model="fixed_pedigree",
            recombination_rate=self.cfg.recombination_rate,
            random_seed=self.cfg.seed
        )

        # 2. Recapitation (coalescent simulation for deep history)
        logger.info("Recapitating with OutOfAfrica_3G09...")
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        
        # Align pedigree populations to the demographic model (CEU)
        tables = ts_ped.dump_tables()
        tables.populations.clear()
        for pop in model.model.populations:
            tables.populations.add_row(metadata=pop.asdict())
        
        # Map all pedigree nodes to the CEU population ID
        ceu_id = [p.name for p in model.model.populations].index("CEU")
        tables.nodes.population = np.full_like(tables.nodes.population, ceu_id)
        
        ts_recap = msprime.sim_ancestry(
            initial_state=tables.tree_sequence(),
            demography=model.model,
            recombination_rate=self.cfg.recombination_rate,
            random_seed=self.cfg.seed
        )

        # 3. Add mutations
        logger.info("Adding mutations...")
        ts_final = msprime.sim_mutations(
            ts_recap,
            rate=self.cfg.mutation_rate,
            random_seed=self.cfg.seed
        )
        return ts_final

    @staticmethod
    def compute_charlevoix_gc(
        ped: Any,
        sorted_ids: List[int],
        charlevoix_ancestors: Set[int],
    ) -> Dict[int, float]:
        """
        Computes expected genetic contribution from Charlevoix founders.
        Requires `sorted_ids` to be topologically sorted (Ancestors -> Descendants).
        """
        contribution: Dict[int, float] = {}
        
        for ind_id in sorted_ids:
            # Base case: Founder from region
            if ind_id in charlevoix_ancestors:
                contribution[ind_id] = 1.0
                continue
                
            # Recursive case: Average of parents
            try:
                ind = ped[ind_id]
                father_contribution = contribution.get(ind.father.ind, 0.0)
                mother_contribution = contribution.get(ind.mother.ind, 0.0)
                contribution[ind_id] = 0.5 * father_contribution + 0.5 * mother_contribution
            except KeyError:
                contribution[ind_id] = 0.0
        
        return contribution


# =============================================================================
# GWAS Engine
# =============================================================================

class GwasRunner:
    """
    Handles genotype processing and association testing using 
    scikit-learn and scipy.stats standard functions.
    """

    @staticmethod
    def extract_genotypes(ts: tskit.TreeSequence, 
                         balsac_map: Dict[int, int], 
                         sample_ids: List[int]) -> np.ndarray:
        ms_ids = [balsac_map[bid] for bid in sample_ids if bid in balsac_map]
        ts_nodes_flat = []
        for mid in ms_ids:
            ts_nodes_flat.extend(ts.individual(mid).nodes)
        H = ts.genotype_matrix(samples=ts_nodes_flat)
        G = H[:, 0::2] + H[:, 1::2]
        # Return as (samples, variants) for scikit-learn compatibility
        return G.T

    @staticmethod
    def filter_genotypes(G: np.ndarray, maf_thresh: float) -> Tuple[np.ndarray, np.ndarray]:
        ac = G.sum(axis=0)
        an = G.shape[0] * 2
        af = ac / an
        maf = np.minimum(af, 1 - af)
        mask = maf > maf_thresh
        return G[:, mask], maf[mask]

    @staticmethod
    def compute_pca(G: np.ndarray, n_components: int) -> np.ndarray:
        # Standardize: (g - mean) / std
        G_std = (G - np.mean(G, axis=0)) / np.std(G, axis=0)
        pca = PCA(n_components=n_components, svd_solver="full")
        return pca.fit_transform(G_std)

    @staticmethod
    def run_association(G: np.ndarray, y: np.ndarray, C: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Runs GWAS using standard scikit-learn functions.
        
        Methodology:
        1. If covariates (C) exist, regress them out of both Genotypes (G) and Phenotype (y)
           using sklearn.LinearRegression.
        2. Test the association between the residuals using sklearn.f_regression.
        """
        # Ensure y is 1D
        y = y.ravel()

        # 1. Handle Covariates via Residualization
        if C is not None:
            # Fit covariate model for phenotype
            cov_model_y = LinearRegression(fit_intercept=True, copy_X=False)
            cov_model_y.fit(C, y)
            y_resid = y - cov_model_y.predict(C)

            # Fit covariate model for genotypes (Multi-output regression)
            cov_model_G = LinearRegression(fit_intercept=True, copy_X=False)
            cov_model_G.fit(C, G)
            G_resid = G - cov_model_G.predict(C)
        else:
            y_resid = y
            G_resid = G

        # 2. Run Vectorized Association (F-test)
        F_scores, p_values = f_regression(G_resid, y_resid)
        
        # 3. Degrees of Freedom Correction
        p_values = np.clip(p_values, 1e-300, 1.0)
        neg_log10_p = -np.log10(p_values)
        
        return np.nan_to_num(neg_log10_p, nan=0.0)

    @staticmethod
    def calculate_lambda_gc(neg_log10_p: np.ndarray) -> float:
        """Calculates Genomic Inflation Factor (Lambda GC)."""
        # Example: https://github.com/pgxcentre/lambda/blob/master/compute_lambda.py
        if len(neg_log10_p) == 0: return 1.0
        p_vals = 10**(-neg_log10_p)
        chi2 = stats.chi2.ppf(1 - p_vals, df=1)
        return np.median(chi2) / stats.chi2.ppf(0.5, df=1)


# =============================================================================
# Visualization
# =============================================================================

def plot_qq(results: Dict[str, np.ndarray], output_path: Path):
    """
    Generates QQ plots with 95% Confidence Intervals.
    """
    sns.set_theme(style="ticks", font_scale=1.1, rc={"font.family": "sans-serif"})
    n_panels = len(results)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.5 * n_panels, 4.5), constrained_layout=True)
    if n_panels == 1: axes = [axes]

    # Find the global max for axis limits except for when model is "Uncorrected"
    global_max = 0
    for model, logp in results.items():
        if model == "Uncorrected":
            continue
        n = len(logp)
        observed = np.sort(logp)[::-1]
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))
        current_max = max(np.max(expected), np.max(observed))
        if current_max > global_max:
            global_max = current_max

    for i, (ax, (model, logp)) in enumerate(zip(axes, results.items())):
        # Add bold panel label (A, B, C...) to top left
        ax.text(0.00, 1.05, chr(65 + i), transform=ax.transAxes, 
                size=16, weight='bold')

        # 1. Prepare Data
        n = len(logp)
        observed = np.sort(logp)[::-1]
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))
        
        # 2. Subsampling for Plotting Efficiency (if N > 10k)
        if n > 10000:
            indices = np.unique(np.concatenate([
                np.arange(0, 1000), 
                np.linspace(1000, n-1, 5000).astype(int)
            ]))
        else:
            indices = np.arange(n)

        exp_plot = expected[indices]
        obs_plot = observed[indices]
        
        # 3. 95% Confidence Interval Calculation
        k = indices + 1
        ci_upper = -np.log10(stats.beta.ppf(0.025, k, n - k + 1))
        ci_lower = -np.log10(stats.beta.ppf(0.975, k, n - k + 1))
        
        # 4. Plotting
        ax.fill_between(exp_plot, ci_lower, ci_upper, color="gray", alpha=0.25, lw=0, label="95% CI")
        ax.scatter(exp_plot, obs_plot, c="#2171b5", s=12, alpha=0.8, edgecolors='none', rasterized=True)

        max_val = max(np.max(exp_plot), np.max(obs_plot)) * 1.05
        if model != "Uncorrected":
            max_val = max(max_val, global_max * 1.05)
        ax.plot([0, max_val], [0, max_val], "r--", lw=1.5)
        
        gc = GwasRunner.calculate_lambda_gc(logp)
        logger.info(f"Model: {model} | Lambda GC: {gc:.3f}")
        ax.set_title(f"{model}\n$\lambda_{{GC}} = {gc:.3f}$")
        ax.set_xlabel(r"Expected $-\log_{10}(P)$")
        if ax == axes[0]: ax.set_ylabel(r"Observed $-\log_{10}(P)$")
        
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        sns.despine(ax=ax)

    plt.savefig(output_path, dpi=300)
    logger.info(f"QQ plot saved to {output_path}")


# =============================================================================
# Main Pipeline
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="GWAS Simulation with Pedigree contribution")
    parser.add_argument("--config", "-c", help="Path to JSON configuration")
    parser.add_argument("--seed", "-s", type=int, default=2026, help="RNG seed")
    args = parser.parse_args()
    
    cfg = SimulationConfig.from_args(args)
    np.random.seed(cfg.seed)
    logger.info(f"Random seed set to: {cfg.seed}")
    
    # --- 1. Load data ---
    loader = DataLoader(cfg.paths_file)
    ped = loader.load_genealogy()
    all_probands = loader.load_probands()
    full_region_map = loader.load_region_map_full()
    
    # Sampling is deterministic given cfg.seed
    sample_ids = np.random.choice(all_probands, size=cfg.sample_size, replace=False).tolist()
    
    # --- 2. Run simulation ---
    sim = GenealogySimulator(ped, cfg)
    pedigree_ms, id_map, sorted_ids = sim.build_pedigree(sample_ids)
    ts = sim.run(pedigree_ms)
    
    # --- 3. Process genotypes ---
    G_raw = GwasRunner.extract_genotypes(ts, id_map, sample_ids)
    G_gwas, _ = GwasRunner.filter_genotypes(G_raw, cfg.maf_threshold_gwas)
    G_pca_pre, _ = GwasRunner.filter_genotypes(G_raw, cfg.maf_threshold_pca)
    
    logger.info(f"LD pruning for PCA (input variants: {G_pca_pre.shape[1]})...")
    gn = np.ascontiguousarray(G_pca_pre.T)
    is_unlinked = allel.locate_unlinked(
        gn, size=cfg.ld_window, step=cfg.ld_step, threshold=cfg.ld_threshold
    )
    G_pca = G_pca_pre[:, is_unlinked]
    logger.info(f"Retained {G_pca.shape[1]} variants for PCA.")
    
    # --- 4. Phenotype construction ---
    logger.info("Computing Charlevoix genetic contributions...")
    
    charlevoix_ancestors = {
        pid for pid in sorted_ids 
        if full_region_map.get(pid) == "Charlevoix"
    }
    logger.info(f"Identified {len(charlevoix_ancestors)} Charlevoix ancestors in the pedigree.")
    
    contribution_dict = sim.compute_charlevoix_gc(ped, sorted_ids, charlevoix_ancestors)
    
    # Extract genetic contribution for samples involved in simulation
    valid_sample_ids = [bid for bid in sample_ids if bid in id_map]
    gc_vector = np.array([contribution_dict.get(bid, 0.0) for bid in valid_sample_ids])
    
    # QC: GC stats
    gc_mean, gc_var = np.mean(gc_vector), np.var(gc_vector)
    logger.info(f"GC stats | Mean: {gc_mean:.4f} | Var: {gc_var:.6f}")
    
    if gc_var < 1e-9:
        logger.warning("WARNING: GC variance is effectively zero. Phenotype will be pure noise.")

    # 1. Structure (standardized)
    structure_z = (gc_vector - np.mean(gc_vector)) / np.std(gc_vector)
    
    # 2. Noise (standardized)
    noise_raw = np.random.normal(0, 1, len(valid_sample_ids))
    noise_z = (noise_raw - np.mean(noise_raw)) / np.std(noise_raw)

    # 3. Combine
    h2 = cfg.h2
    pheno = (structure_z * np.sqrt(h2)) + (noise_z * np.sqrt(1 - h2))
    
    # Calculate regression statistics
    slope, intercept, r_value, p_value, std_err = stats.linregress(structure_z, pheno)
    logger.info(f"Phenotype construction | Pearson R^2: {r_value**2:.4f} | Slope: {slope:.4f} | p-value: {p_value:.2e}")
    
    # --- 5. Corrections ---
    logger.info(f"Computing {cfg.n_pcs2} genetic PCs...")
    pcs2 = GwasRunner.compute_pca(G_pca, cfg.n_pcs2)

    logger.info(f"Computing 2 genealogical PCoA coordinates...")
    try:
        if len(valid_sample_ids) < 5000:
            phi = gk.phi(ped, pro=valid_sample_ids, verbose=False)
        else:
            with open(cfg.paths_file) as f:
                paths = json.load(f)
            with open(Path(paths["wd"]) / "results/pickles/balsac_kinship_nocousins.pkl", "rb") as f:
                phi = pickle.load(f)
                phi = phi.loc[valid_sample_ids, valid_sample_ids]
        dist = 1 - phi.values
        np.fill_diagonal(dist, 0)
        dist = (dist + dist.T) / 2
        dist[dist < 0] = 0
        
        pcoa_res = pcoa(DistanceMatrix(dist, ids=[str(i) for i in valid_sample_ids]))
        pcoa_coords = pcoa_res.samples.iloc[:, :2].values

        slope, intercept, r_value, p_value, std_err = stats.linregress(structure_z, pcoa_coords[:, 0])
        logger.info(f"Charlevoix GC and PCoA 1 | Pearson R^2: {r_value**2:.4f} | Slope: {slope:.4f} | p-value: {p_value:.2e}")
        slope, intercept, r_value, p_value, std_err = stats.linregress(structure_z, pcs2[:, 0])
        logger.info(f"Charlevoix GC and PCA 1 | Pearson R^2: {r_value**2:.4f} | Slope: {slope:.4f} | p-value: {p_value:.2e}")

    except Exception as e:
        logger.error(f"PCoA failed: {e}. Skipping PCoA model.")
        pcoa_coords = None
    
    # --- 6. Run GWAS ---
    logger.info("Running GWAS Models...")
    results = {}
    
    # Baseline (null hypothesis)
    y_null = np.random.normal(0, 1, len(valid_sample_ids))
    results["Baseline (null)"] = GwasRunner.run_association(G_gwas, y_null)
    
    # Models
    results["Uncorrected"] = GwasRunner.run_association(G_gwas, pheno)
    results["Genetic PCA (10 PCs)"] = GwasRunner.run_association(G_gwas, pheno, pcs2[:, :cfg.n_pcs1])
    results["Genetic PCA (100 PCs)"] = GwasRunner.run_association(G_gwas, pheno, pcs2)
    
    if pcoa_coords is not None:
        results["Genealogical PCoA"] = GwasRunner.run_association(G_gwas, pheno, pcoa_coords)
    
    # --- 7. Visualization & export ---
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
        
    # Plotting QQ
    plot_qq(results, cfg.output_dir / "7.GWASSimulationQQ.png")
        
    logger.info("Simulation pipeline completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as e:
        logger.exception("An error occurred during execution.")
        sys.exit(1)