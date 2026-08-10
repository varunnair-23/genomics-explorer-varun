from pathlib import Path
import logging
import yaml

import scanpy as sc
import cellxgene_census


def load_config(config_path):
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def setup_logging(log_dir):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / "pipeline.log"

    logger = logging.getLogger("census_pipeline")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] - %(message)s"
    )

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def build_value_filter(census_config):
    tissue = census_config["tissue_general"]
    disease = census_config["disease"]
    assay = census_config["assay"]

    value_filter = (
        f"tissue_general == '{tissue}' "
        f"and disease == '{disease}' "
        f"and assay == \"{assay}\" "
        f"and is_primary_data == True"
    )

    return value_filter


def main():
    config = load_config("config/query_config.yaml")

    logger = setup_logging(config["paths"]["logs"])

    logger.info("Starting CELLxGENE Census ingestion pipeline.")
    logger.info("Project: %s", config["project"]["name"])
    logger.info("Phase: %s", config["project"]["phase"])

    processed_dir = Path(config["paths"]["processed_data"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    census_config = config["census"]
    qc_config = config["qc"]

    census_version = census_config["census_version"]
    organism_key = census_config["organism_key"]
    organism = census_config["organism"]
    measurement_name = census_config["measurement_name"]
    max_cells = census_config["max_cells"]

    value_filter = build_value_filter(census_config)

    logger.info("Using Census version: %s", census_version)
    logger.info("Metadata filter: %s", value_filter)

    with cellxgene_census.open_soma(census_version=census_version) as census:
        logger.info("Opened CELLxGENE Census.")

        obs = (
            census["census_data"][organism_key]
            .obs
            .read(
                value_filter=value_filter,
                column_names=[
                    "soma_joinid",
                    "cell_type",
                    "tissue_general",
                    "tissue",
                    "disease",
                    "assay",
                    "is_primary_data",
                ],
            )
            .concat()
            .to_pandas()
        )

        logger.info("Matched cells before limiting: %s", obs.shape[0])

        if obs.empty:
            raise ValueError("No cells matched the Census query filter.")

        selected_obs = obs.head(max_cells)
        selected_ids = selected_obs["soma_joinid"].tolist()

        logger.info("Selected %s cells for AnnData streaming.", len(selected_ids))

        id_filter = "soma_joinid in [" + ", ".join(map(str, selected_ids)) + "]"

        adata = cellxgene_census.get_anndata(
            census=census,
            organism=organism,
            measurement_name=measurement_name,
            obs_value_filter=id_filter,
            obs_column_names=[
                "soma_joinid",
                "cell_type",
                "tissue_general",
                "tissue",
                "disease",
                "assay",
                "is_primary_data",
            ],
            var_column_names=[
                "feature_id",
                "feature_name",
            ],
        )

    logger.info("Streamed AnnData object: %s", adata)

    adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True
    )

    logger.info("Calculated QC metrics.")

    min_genes = qc_config["min_genes_per_cell"]
    max_mito = qc_config["max_mito_percent"]
    normalize_target_sum = qc_config["normalize_target_sum"]

    before_qc = adata.n_obs

    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] <= max_mito].copy()

    after_qc = adata.n_obs

    logger.info("Cells before QC filtering: %s", before_qc)
    logger.info("Cells after QC filtering: %s", after_qc)

    sc.pp.normalize_total(
        adata,
        target_sum=normalize_target_sum
    )

    sc.pp.log1p(adata)

    logger.info("Applied total-count normalization and log1p transformation.")

    h5ad_path = processed_dir / "census_subset.h5ad"
    metadata_path = processed_dir / "cell_metadata_qc.csv"

    adata.write_h5ad(h5ad_path)
    adata.obs.to_csv(metadata_path)

    logger.info("Saved cleaned AnnData file to %s", h5ad_path)
    logger.info("Saved cell metadata CSV to %s", metadata_path)
    logger.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()