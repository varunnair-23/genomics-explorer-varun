from pathlib import Path
import logging
import yaml

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors


def load_yaml(path):
    with open(path, "r") as file:
        return yaml.safe_load(file)


def setup_logging(log_dir):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("chemical_structure_pipeline")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] - %(message)s")

    file_handler = logging.FileHandler(log_dir / "chemical_structure_pipeline.log")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def smiles_to_molecule(smiles, compound_name):
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError(f"Invalid SMILES string for {compound_name}")

    return molecule


def calculate_descriptors(molecule):
    return {
        "molecular_weight": Descriptors.MolWt(molecule),
        "logp": Descriptors.MolLogP(molecule),
        "h_donors": Descriptors.NumHDonors(molecule),
        "h_acceptors": Descriptors.NumHAcceptors(molecule),
        "rotatable_bonds": Descriptors.NumRotatableBonds(molecule),
    }


def save_structure_image(molecule, output_path):
    Draw.MolToFile(
        molecule,
        str(output_path),
        size=(600, 400)
    )


def main():
    project_config = load_yaml("config/query_config.yaml")
    compound_config = load_yaml("config/compounds.yaml")

    logger = setup_logging(project_config["paths"]["logs"])

    output_dir = Path(project_config["paths"]["results"]) / "chemical_structures"
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    logger.info("Starting chemical structure generation pipeline.")

    for compound in compound_config["compounds"]:
        name = compound["name"]
        smiles = compound["smiles"]

        logger.info("Processing compound: %s", name)

        molecule = smiles_to_molecule(smiles, name)

        descriptors = calculate_descriptors(molecule)

        image_filename = name.lower().replace(" ", "_") + ".png"
        image_path = output_dir / image_filename

        save_structure_image(molecule, image_path)

        summary_rows.append({
            "name": name,
            "class": compound["class"],
            "target": compound["target"],
            "pathway": compound["pathway"],
            "smiles": smiles,
            "structure_image": str(image_path),
            **descriptors,
        })

        logger.info("Saved structure image to %s", image_path)

    summary_df = pd.DataFrame(summary_rows)

    summary_path = Path(project_config["paths"]["processed_data"]) / "compound_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    logger.info("Saved compound summary to %s", summary_path)
    logger.info("Chemical structure generation pipeline completed successfully.")


if __name__ == "__main__":
    main()