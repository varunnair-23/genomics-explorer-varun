"""Generate UMAP showcase plots colored by treatment and QC metrics."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_loader import load_yaml, resolve_path, setup_logging


def main() -> None:
    config = load_yaml("config/query_config.yaml")
    logger = setup_logging(config["paths"]["logs"], "plot_umap", "plot_umap.log")

    umap_path = resolve_path("data/processed/umap_coordinates.csv")
    results_dir = resolve_path(config["paths"]["results"])
    results_dir.mkdir(parents=True, exist_ok=True)

    if not umap_path.exists():
        raise FileNotFoundError(
            f"Missing {umap_path}. Run src/cluster_manifold.py first."
        )

    umap_df = pd.read_csv(umap_path)
    required = {"UMAP_1", "UMAP_2", "treatment_group", "percent_mito"}
    missing = required - set(umap_df.columns)
    if missing:
        raise ValueError(f"UMAP file missing columns: {missing}")

    sns.set_theme(style="white")
    palette = {
        "control": "#4C4C4C",
        "Nutlin": "#C0392B",
        "SAHA": "#2471A3",
        "BMS": "#1E8449",
    }

    fig, ax = plt.subplots(figsize=(10, 7))
    for group, frame in umap_df.groupby("treatment_group"):
        ax.scatter(
            frame["UMAP_1"],
            frame["UMAP_2"],
            s=8,
            alpha=0.75,
            label=group,
            c=palette.get(str(group), None),
        )
    ax.set_title("sci-Plex2 A549 UMAP by chemical perturbation")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(title="Treatment", markerscale=2, frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(left=True, bottom=True)
    treatment_path = results_dir / "umap_by_treatment.png"
    fig.tight_layout()
    fig.savefig(treatment_path, dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 7))
    points = ax.scatter(
        umap_df["UMAP_1"],
        umap_df["UMAP_2"],
        c=umap_df["percent_mito"],
        s=8,
        cmap="viridis",
        alpha=0.8,
    )
    cbar = fig.colorbar(points, ax=ax)
    cbar.set_label("Mitochondrial percent")
    ax.set_title("UMAP QC overlay: mitochondrial content")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(left=True, bottom=True)
    qc_path = results_dir / "umap_qc_mito.png"
    fig.tight_layout()
    fig.savefig(qc_path, dpi=300)
    plt.close(fig)

    logger.info("Saved %s", treatment_path)
    logger.info("Saved %s", qc_path)
    logger.info("Plotting completed successfully.")


if __name__ == "__main__":
    main()