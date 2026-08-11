"""Ingest sci-Plex2 chemical perturbation data, QC-filter, and normalize."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Allow: python src/ingest_perturbation.py  (from repo root)
sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_loader import (
    load_yaml,
    resolve_path,
    setup_logging,
    validate_perturbation_adata,
)


def _to_float_dose(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str), errors="coerce")


def subsample_by_group(adata, group_col: str, max_cells_per_group: int, random_state: int):
    if max_cells_per_group <= 0:
        return adata

    keep_indices = []
    rng = np.random.default_rng(random_state)

    for _, idx in adata.obs.groupby(group_col, observed=True).groups.items():
        idx = list(idx)
        if len(idx) > max_cells_per_group:
            idx = list(rng.choice(idx, size=max_cells_per_group, replace=False))
        keep_indices.extend(idx)

    return adata[keep_indices].copy()


def main() -> None:
    config = load_yaml("config/query_config.yaml")
    logger = setup_logging(
        config["paths"]["logs"],
        "ingest_pipeline",
        "ingest_pipeline.log",
    )

    raw_dir = resolve_path(config["paths"]["raw_data"])
    processed_dir = resolve_path(config["paths"]["processed_data"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    dataset = config["dataset"]
    qc = config["qc"]

    input_path = raw_dir / dataset["filename"]
    if not input_path.exists():
        raise FileNotFoundError(
            f"Missing raw dataset: {input_path}\n"
            f"Run: python src/download_data.py\n"
            f"Expected URL: {dataset['download_url']}"
        )

    logger.info("Reading perturbation AnnData from %s", input_path)
    adata = sc.read_h5ad(input_path)
    validate_perturbation_adata(adata)

    logger.info("Loaded AnnData: %s", adata)
    logger.info("Disease labels: %s", sorted(adata.obs["disease"].astype(str).unique()))
    logger.info("Cell lines: %s", sorted(adata.obs["cell_line"].astype(str).unique()))

    focus = set(dataset["focus_perturbations"])
    adata = adata[adata.obs[dataset["perturbation_column"]].isin(focus)].copy()
    logger.info("Cells after focusing on %s: %s", sorted(focus), adata.n_obs)

    adata.obs["dose_value_num"] = _to_float_dose(adata.obs[dataset["dose_column"]])
    adata.obs["perturbation"] = adata.obs[dataset["perturbation_column"]].astype(str)
    adata.obs["is_control"] = adata.obs["perturbation"] == dataset["control_label"]

    min_dose = float(dataset["min_dose_for_treatment"])
    keep_mask = adata.obs["is_control"] | (
        (~adata.obs["is_control"]) & (adata.obs["dose_value_num"] >= min_dose)
    )
    adata = adata[keep_mask].copy()
    logger.info(
        "Cells after dose filter (control or dose >= %s): %s",
        min_dose,
        adata.n_obs,
    )

    # Unified label used by later DE / plotting
    adata.obs["treatment_group"] = np.where(
        adata.obs["is_control"],
        "control",
        adata.obs["perturbation"].astype(str),
    )

    before_qc = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=qc["min_genes_per_cell"])
    adata = adata[adata.obs["ngenes"] <= qc["max_genes_per_cell"]].copy()
    adata = adata[adata.obs["percent_mito"] <= qc["max_mito_percent"]].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    after_qc = adata.n_obs
    logger.info("Cells before QC: %s | after QC: %s", before_qc, after_qc)

    if after_qc < 100:
        raise ValueError("Too few cells remain after QC to continue.")

    adata = subsample_by_group(
        adata,
        group_col="treatment_group",
        max_cells_per_group=int(dataset["max_cells_per_group"]),
        random_state=int(dataset["random_state"]),
    )
    logger.info("Cells after stratified subsample: %s", adata.n_obs)
    logger.info(
        "Group sizes:\n%s",
        adata.obs["treatment_group"].value_counts().to_string(),
    )

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=qc["normalize_target_sum"])
    sc.pp.log1p(adata)

    h5ad_path = processed_dir / "perturbation_qc.h5ad"
    metadata_path = processed_dir / "cell_metadata_qc.csv"

    adata.write_h5ad(h5ad_path)
    adata.obs.to_csv(metadata_path)

    logger.info("Saved cleaned AnnData to %s", h5ad_path)
    logger.info("Saved cell metadata to %s", metadata_path)
    logger.info("Ingestion completed successfully.")


if __name__ == "__main__":
    main()