from pathlib import Path
import logging
import yaml

import scanpy as sc

try:
    import cellxgene_census
    census_available = True
except ImportError:
    census_available = False

def main():
    with open("config/settings.yml", "r") as file:
        config = yaml.safe_load(file)

        log_dir = Path(config["paths"]["logs"])
        log_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(
            filename=log_dir / "check_single_cell_setup.log",
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s"
        )

        logging.info("Starting single-cell setup check")

        print("Project:", config["project"]["name"])
        print("Phase:", config["project"]["phase"])
        print("Scanpy version:", sc.__version__)

        if census_available:
            print("cellxgene census is available")
            logging.info("cellxgene census is available")
        else:
            print("cellxgene census is not available")
            logging.warning("cellxgene census is not available")

        logging.info("Single-cell setup check completed successfully")

if __name__ == "__main__":
            main()