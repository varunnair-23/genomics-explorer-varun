from pathlib import Path
import logging
import time

import pandas as pd
import requests
import yaml


ENRICHR_ADDLIST_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"


def load_config(config_path):
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def setup_logging(log_dir):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("pathway_pipeline")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] - %(message)s")

    file_handler = logging.FileHandler(log_dir / "pathway_pipeline.log")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def submit_gene_list_to_enrichr(genes, description):
    payload = {
        "list": "\n".join(genes),
        "description": description,
    }

    response = requests.post(
        ENRICHR_ADDLIST_URL,
        files=payload,
        timeout=30,
    )

    response.raise_for_status()

    result = response.json()

    if "userListId" not in result:
        raise ValueError(f"Enrichr did not return a userListId: {result}")

    return result["userListId"]


def fetch_enrichment_results(user_list_id, library):
    response = requests.get(
        ENRICHR_ENRICH_URL,
        params={
            "userListId": user_list_id,
            "backgroundType": library,
        },
        timeout=30,
    )

    response.raise_for_status()

    results = response.json()

    if library not in results:
        raise ValueError(f"No results returned for library: {library}")

    return results[library]


def parse_enrichr_results(raw_results, cluster_id):
    rows = []

    for result in raw_results:
        rows.append({
            "cluster": cluster_id,
            "rank": result[0],
            "pathway": result[1],
            "p_value": result[2],
            "z_score": result[3],
            "combined_score": result[4],
            "overlapping_genes": ";".join(result[5]),
            "adjusted_p_value": result[6],
        })

    return rows


def main():
    config = load_config("config/query_config.yaml")

    logger = setup_logging(config["paths"]["logs"])

    pathway_config = config["pathway"]

    marker_file = Path(pathway_config["marker_file"])
    output_file = Path(pathway_config["output_file"])

    gene_column = pathway_config["gene_column"]
    cluster_column = pathway_config["cluster_column"]
    adjusted_pvalue_cutoff = pathway_config["adjusted_pvalue_cutoff"]
    top_genes_per_cluster = pathway_config["top_genes_per_cluster"]
    enrichr_library = pathway_config["enrichr_library"]

    logger.info("Starting pathway enrichment pipeline.")
    logger.info("Reading marker genes from %s", marker_file)

    if not marker_file.exists():
        raise FileNotFoundError(f"Marker gene file not found: {marker_file}")

    markers = pd.read_csv(marker_file)

    required_columns = {gene_column, cluster_column, "pvals_adj"}
    missing_columns = required_columns - set(markers.columns)

    if missing_columns:
        raise ValueError(f"Marker file is missing columns: {missing_columns}")

    significant_markers = markers[
        markers["pvals_adj"] < adjusted_pvalue_cutoff
    ].copy()

    logger.info("Total marker genes: %s", markers.shape[0])
    logger.info("Significant marker genes: %s", significant_markers.shape[0])

    all_results = []

    for cluster_id, cluster_df in significant_markers.groupby(cluster_column):
        top_genes = (
            cluster_df
            .sort_values("scores", ascending=False)
            .head(top_genes_per_cluster)[gene_column]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        if len(top_genes) < 5:
            logger.warning(
                "Skipping cluster %s because it has fewer than 5 significant genes.",
                cluster_id
            )
            continue

        logger.info(
            "Submitting cluster %s with %s genes to Enrichr.",
            cluster_id,
            len(top_genes)
        )

        user_list_id = submit_gene_list_to_enrichr(
            genes=top_genes,
            description=f"Leiden cluster {cluster_id} marker genes"
        )

        raw_results = fetch_enrichment_results(
            user_list_id=user_list_id,
            library=enrichr_library
        )

        parsed_rows = parse_enrichr_results(
            raw_results=raw_results,
            cluster_id=cluster_id
        )

        all_results.extend(parsed_rows)

        time.sleep(1)

    if not all_results:
        raise ValueError("No pathway enrichment results were generated.")

    enrichment_df = pd.DataFrame(all_results)

    enrichment_df = enrichment_df.sort_values(
        ["cluster", "adjusted_p_value", "combined_score"],
        ascending=[True, True, False]
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    enrichment_df.to_csv(output_file, index=False)

    logger.info("Saved pathway enrichment results to %s", output_file)
    logger.info("Pathway enrichment pipeline completed successfully.")


if __name__ == "__main__":
    main()