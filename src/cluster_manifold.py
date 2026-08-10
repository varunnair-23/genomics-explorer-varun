from pathlib import Path
import logging
import yaml

import scanpy as sc
import pandas as pd

def load_config(config_path):
    with open(config_path, "r") as file:
        return yaml.safe_load(file)

def setup_logging(log_dir):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("manifold_pipeline")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] - %(message)s")

    file_handler = logging.FileHandler(log_dir / "manifold_pipeline.log")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def main():
    config = load_config("config/query_config.yaml")

    logger = setup_logging(config["paths"]["logs"])

    processed_dir = Path(config["paths"]["processed_data"])
    results_dir = Path(config["paths"]["results"])

    processed_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    manifold_config = config["manifold"]

    n_top_genes = manifold_config["n_top_genes"]
    n_pcs = manifold_config["n_pcs"]
    n_neighbors = manifold_config["n_neighbors"]
    leiden_resolution = manifold_config["leiden_resolution"]
    random_state = manifold_config["random_state"]

    input_path = processed_dir / "census_subset.h5ad"
    clustered_path = processed_dir / "clustered_census.h5ad"
    umap_csv_path = processed_dir / "umap_coordinates.csv"
    markers_csv_path = processed_dir / "cluster_markers.csv"

    logger.info("Starting manifold clustering pipeline.")
    logger.info("Reading AnnData from %s", input_path)

    adata = sc.read_h5ad(input_path)

    logger.info("Loaded AnnData: %s", adata)

    if "feature_name" in adata.var.columns:
        adata.var["gene_symbol"] = adata.var["feature_name"].astype(str)
        adata.var_names = adata.var["gene_symbol"]
        adata.var_names_make_unique()
        adata.var.index.name = None
        logger.info("Set AnnData var_names to gene symbols from feature_name.")
    else:
        raise ValueError("feature_name column not found in adata.var. Cannot run marker gene analysis.")

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes
    )

    adata = adata[:, adata.var["highly_variable"]].copy()

    logger.info("AnnData after HVG filtering: %s", adata)

    logger.info("Scaling expression matrix.")

    sc.pp.scale(
        adata,
        max_value=10
    )

    logger.info("Running PCA.")

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        random_state=random_state
    )

    logger.info("Saving PCA elbow plot.")

    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=50,
        log=True,
        show=False,
        save="_midway_elbow.png"
    )

    logger.info("Building k-nearest neighbor graph with n_neighbors=%s and n_pcs=%s.", n_neighbors, n_pcs)

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        random_state=random_state
    )

    logger.info("Running Leiden clustering with resolution=%s.", leiden_resolution)

    sc.tl.leiden(
        adata,
        resolution=leiden_resolution,
        random_state=random_state,
        key_added="leiden"
    )

    cluster_count = adata.obs["leiden"].nunique()
    logger.info("Detected %s Leiden clusters.", cluster_count)

    logger.info("Running UMAP.")

    sc.tl.umap(
        adata,
        random_state=random_state
    )

    logger.info("Running marker gene ranking with Wilcoxon test.")

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )

    markers = sc.get.rank_genes_groups_df(
        adata,
        group=None
    )

    markers.to_csv(markers_csv_path, index=False)

    logger.info("Saved marker genes to %s", markers_csv_path)

    umap_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "UMAP_1": adata.obsm["X_umap"][:, 0],
        "UMAP_2": adata.obsm["X_umap"][:, 1],
        "leiden": adata.obs["leiden"].astype(str).values,
        "total_counts": adata.obs["total_counts"].values,
        "n_genes_by_counts": adata.obs["n_genes_by_counts"].values,
        "pct_counts_mt": adata.obs["pct_counts_mt"].values,
    })

    if "cell_type" in adata.obs.columns:
        umap_df["cell_type"] = adata.obs["cell_type"].values

    umap_df.to_csv(umap_csv_path, index=False)

    logger.info("Saved UMAP coordinates to %s", umap_csv_path)

    adata.write_h5ad(clustered_path)

    logger.info("Saved clustered AnnData to %s", clustered_path)
    logger.info("Manifold clustering pipeline completed successfully.")


if __name__ == "__main__":
    main()