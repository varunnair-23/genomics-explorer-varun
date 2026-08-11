"""PCA / Leiden / UMAP manifold for chemical perturbation cells."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_loader import load_yaml, resolve_path, setup_logging


def main() -> None:
    config = load_yaml("config/query_config.yaml")
    logger = setup_logging(
        config["paths"]["logs"],
        "manifold_pipeline",
        "manifold_pipeline.log",
    )

    processed_dir = resolve_path(config["paths"]["processed_data"])
    results_dir = resolve_path(config["paths"]["results"])
    results_dir.mkdir(parents=True, exist_ok=True)

    manifold = config["manifold"]
    input_path = processed_dir / "perturbation_qc.h5ad"
    clustered_path = processed_dir / "clustered_perturbation.h5ad"
    umap_csv_path = processed_dir / "umap_coordinates.csv"

    if not input_path.exists():
        raise FileNotFoundError(
            f"Missing {input_path}. Run src/ingest_perturbation.py first."
        )

    logger.info("Reading AnnData from %s", input_path)
    adata = sc.read_h5ad(input_path)
    logger.info("Loaded AnnData: %s", adata)

    # 1) Find informative genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=manifold["n_top_genes"],
        flavor="seurat",
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    logger.info("AnnData after HVG filtering: %s", adata)

    # 2) Scale + PCA
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(
        adata,
        svd_solver="arpack",
        random_state=manifold["random_state"],
    )

    # 3) Neighbor graph in PCA space
    sc.pp.neighbors(
        adata,
        n_neighbors=manifold["n_neighbors"],
        n_pcs=manifold["n_pcs"],
        random_state=manifold["random_state"],
    )

    # 4) Leiden clusters (on PCA graph, not UMAP)
    sc.tl.leiden(
        adata,
        resolution=manifold["leiden_resolution"],
        random_state=manifold["random_state"],
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
    )

    # 5) UMAP for visualization
    sc.tl.umap(adata, random_state=manifold["random_state"])

    umap_df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "UMAP_1": adata.obsm["X_umap"][:, 0],
            "UMAP_2": adata.obsm["X_umap"][:, 1],
            "leiden": adata.obs["leiden"].astype(str).values,
            "treatment_group": adata.obs["treatment_group"].astype(str).values,
            "perturbation": adata.obs["perturbation"].astype(str).values,
            "dose_value": adata.obs["dose_value_num"].values,
            "percent_mito": adata.obs["percent_mito"].values,
            "ngenes": adata.obs["ngenes"].values,
            "cell_line": adata.obs["cell_line"].astype(str).values,
            "disease": adata.obs["disease"].astype(str).values,
        }
    )
    umap_df.to_csv(umap_csv_path, index=False)
    adata.write_h5ad(clustered_path)

    logger.info("Saved UMAP coordinates to %s", umap_csv_path)
    logger.info("Saved clustered AnnData to %s", clustered_path)
    logger.info("Manifold clustering pipeline completed successfully.")


if __name__ == "__main__":
    main()