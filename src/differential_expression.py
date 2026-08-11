"""Differential expression: each drug treatment versus vehicle control."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_loader import load_yaml, resolve_path, setup_logging


def main() -> None:
    config = load_yaml("config/query_config.yaml")
    compounds_config = load_yaml("config/compounds.yaml")
    logger = setup_logging(
        config["paths"]["logs"],
        "de_pipeline",
        "differential_expression.log",
    )

    de_cfg = config["differential_expression"]
    processed_dir = resolve_path(config["paths"]["processed_data"])
    input_path = processed_dir / "perturbation_qc.h5ad"
    output_path = resolve_path(de_cfg["output_file"])

    if not input_path.exists():
        raise FileNotFoundError(
            f"Missing {input_path}. Run src/ingest_perturbation.py first."
        )

    # Map dataset labels (Nutlin) -> display names (Nutlin-3A)
    label_to_name = {
        c["perturbation_label"]: c["name"] for c in compounds_config["compounds"]
    }

    logger.info("Reading AnnData from %s", input_path)
    adata = sc.read_h5ad(input_path)
    adata.var_names_make_unique()

    groups = [
        g
        for g in adata.obs["treatment_group"].astype(str).unique()
        if g != "control"
    ]
    if not groups:
        raise ValueError("No treatment groups found for differential expression.")

    all_rows = []

    for group in sorted(groups):
        subset = adata[adata.obs["treatment_group"].isin(["control", group])].copy()
        counts = subset.obs["treatment_group"].value_counts()

        if counts.get("control", 0) < de_cfg["min_cells_per_group"]:
            logger.warning("Skipping %s: too few control cells.", group)
            continue
        if counts.get(group, 0) < de_cfg["min_cells_per_group"]:
            logger.warning("Skipping %s: too few treated cells.", group)
            continue

        logger.info(
            "Running DE for %s vs control (%s vs %s cells).",
            group,
            counts.get(group, 0),
            counts.get("control", 0),
        )

        sc.tl.rank_genes_groups(
            subset,
            groupby="treatment_group",
            groups=[group],
            reference="control",
            method=de_cfg["method"],
            use_raw=False,
        )

        result = sc.get.rank_genes_groups_df(subset, group=group)
        result["perturbation_label"] = group
        result["compound"] = label_to_name.get(group, group)
        result["direction"] = np.where(result["logfoldchanges"] >= 0, "up", "down")

        n_sig = result[
            (result["pvals_adj"] < de_cfg["adjusted_pvalue_cutoff"])
            & (result["logfoldchanges"].abs() >= de_cfg["min_logfoldchange"])
        ].shape[0]

        logger.info(
            "%s significant DEGs for %s (adj.p < %s, |logFC| >= %s).",
            n_sig,
            group,
            de_cfg["adjusted_pvalue_cutoff"],
            de_cfg["min_logfoldchange"],
        )
        all_rows.append(result)

    if not all_rows:
        raise ValueError("No differential expression results were generated.")

    deg_df = pd.concat(all_rows, ignore_index=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    deg_df.to_csv(output_path, index=False)

    sig_path = processed_dir / "drug_response_degs_significant.csv"
    sig_df = deg_df[
        (deg_df["pvals_adj"] < de_cfg["adjusted_pvalue_cutoff"])
        & (deg_df["logfoldchanges"].abs() >= de_cfg["min_logfoldchange"])
    ].copy()
    sig_df.to_csv(sig_path, index=False)

    logger.info("Saved full DE table to %s", output_path)
    logger.info("Saved significant DE table to %s", sig_path)
    logger.info("Differential expression completed successfully.")


if __name__ == "__main__":
    main()