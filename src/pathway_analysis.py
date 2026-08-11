"""KEGG pathway enrichment for drug-response DEGs via Enrichr."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_loader import load_yaml, resolve_path, setup_logging


ENRICHR_ADDLIST_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"


def submit_gene_list_to_enrichr(genes: list[str], description: str) -> int:
    if len(genes) < 5:
        raise ValueError("Enrichr requires at least 5 genes.")

    payload = {
        "list": "\n".join(genes),
        "description": description,
    }
    response = requests.post(ENRICHR_ADDLIST_URL, files=payload, timeout=60)
    response.raise_for_status()
    result = response.json()
    if "userListId" not in result:
        raise ValueError(f"Enrichr did not return a userListId: {result}")
    return result["userListId"]


def fetch_enrichment_results(user_list_id: int, library: str):
    response = requests.get(
        ENRICHR_ENRICH_URL,
        params={"userListId": user_list_id, "backgroundType": library},
        timeout=60,
    )
    response.raise_for_status()
    results = response.json()
    if library not in results:
        raise ValueError(f"No results returned for library: {library}")
    return results[library]


def parse_enrichr_results(raw_results, compound: str, direction: str):
    rows = []
    for result in raw_results:
        rows.append(
            {
                "compound": compound,
                "direction": direction,
                "rank": result[0],
                "pathway": result[1],
                "p_value": result[2],
                "z_score": result[3],
                "combined_score": result[4],
                "overlapping_genes": ";".join(result[5]),
                "adjusted_p_value": result[6],
            }
        )
    return rows


def is_ribosomal_gene(symbol: str) -> bool:
    gene = symbol.upper()
    return (
        gene.startswith("RPS")
        or gene.startswith("RPL")
        or gene.startswith("MRPS")
        or gene.startswith("MRPL")
    )


def main() -> None:
    config = load_yaml("config/query_config.yaml")
    logger = setup_logging(
        config["paths"]["logs"],
        "pathway_pipeline",
        "pathway_pipeline.log",
    )

    pathway_cfg = config["pathway"]
    de_cfg = config["differential_expression"]

    deg_file = resolve_path(pathway_cfg["deg_file"])
    output_file = resolve_path(pathway_cfg["output_file"])

    if not deg_file.exists():
        raise FileNotFoundError(
            f"DEG file not found: {deg_file}. Run src/differential_expression.py first."
        )

    logger.info("Reading DEGs from %s", deg_file)
    degs = pd.read_csv(deg_file)

    required = {
        pathway_cfg["gene_column"],
        pathway_cfg["group_column"],
        pathway_cfg["direction_column"],
        "pvals_adj",
        "logfoldchanges",
    }
    missing = required - set(degs.columns)
    if missing:
        raise ValueError(f"DEG file is missing columns: {missing}")

    significant = degs[
        (degs["pvals_adj"] < pathway_cfg["adjusted_pvalue_cutoff"])
        & (degs["logfoldchanges"].abs() >= de_cfg["min_logfoldchange"])
    ].copy()

    logger.info("Significant DEGs available for enrichment: %s", significant.shape[0])

    all_results = []
    group_col = pathway_cfg["group_column"]
    direction_col = pathway_cfg["direction_column"]
    gene_col = pathway_cfg["gene_column"]

    for (compound, direction), group_df in significant.groupby(
        [group_col, direction_col], dropna=False
    ):
        ranked = group_df.sort_values("scores", key=lambda s: s.abs(), ascending=False)

        top_genes = []
        seen = set()
        for gene in ranked[gene_col].dropna().astype(str).tolist():
            if is_ribosomal_gene(gene):
                continue
            if gene in seen:
                continue
            seen.add(gene)
            top_genes.append(gene)
            if len(top_genes) >= pathway_cfg["top_genes_per_group"]:
                break

        if len(top_genes) < 5:
            logger.warning(
                "Skipping %s / %s: fewer than 5 significant genes.",
                compound,
                direction,
            )
            continue

        logger.info(
            "Submitting %s %s-regulated genes for %s to Enrichr.",
            len(top_genes),
            direction,
            compound,
        )

        user_list_id = submit_gene_list_to_enrichr(
            genes=top_genes,
            description=f"{compound} {direction}-regulated DEGs vs control",
        )
        raw_results = fetch_enrichment_results(
            user_list_id=user_list_id,
            library=pathway_cfg["enrichr_library"],
        )
        all_results.extend(
            parse_enrichr_results(raw_results, compound=compound, direction=direction)
        )
        time.sleep(1)  # be polite to the public API

    if not all_results:
        raise ValueError("No pathway enrichment results were generated.")

    enrichment_df = pd.DataFrame(all_results)
    enrichment_df = enrichment_df.sort_values(
        ["compound", "direction", "adjusted_p_value", "combined_score"],
        ascending=[True, True, True, False],
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    enrichment_df.to_csv(output_file, index=False)

    sig_path = resolve_path("data/processed/pathway_enrichment_significant.csv")
    enrichment_df[enrichment_df["adjusted_p_value"] < 0.05].to_csv(sig_path, index=False)

    logger.info("Saved pathway enrichment results to %s", output_file)
    logger.info("Saved significant pathways to %s", sig_path)
    logger.info("Pathway enrichment pipeline completed successfully.")


if __name__ == "__main__":
    main()