"""Translate compound SMILES into 2D structures and physicochemical descriptors."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdMolDescriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_loader import load_yaml, resolve_path, setup_logging


def smiles_to_molecule(smiles: str, compound_name: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES string for {compound_name}: {smiles}")
    return molecule


def calculate_descriptors(molecule):
    return {
        "molecular_weight": round(Descriptors.MolWt(molecule), 3),
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "logp": round(Descriptors.MolLogP(molecule), 3),
        "h_donors": int(Descriptors.NumHDonors(molecule)),
        "h_acceptors": int(Descriptors.NumHAcceptors(molecule)),
        "rotatable_bonds": int(Descriptors.NumRotatableBonds(molecule)),
    }


def validate_against_expected(compound: dict, descriptors: dict, logger) -> None:
    expected_mw = compound.get("expected_mw")
    expected_formula = compound.get("expected_formula")

    if expected_formula and descriptors["formula"] != expected_formula:
        raise ValueError(
            f"{compound['name']} formula mismatch: got {descriptors['formula']}, "
            f"expected {expected_formula}. Check the SMILES string."
        )

    if expected_mw is not None:
        if abs(descriptors["molecular_weight"] - float(expected_mw)) > 1.5:
            raise ValueError(
                f"{compound['name']} MW mismatch: got {descriptors['molecular_weight']}, "
                f"expected ~{expected_mw}. Check the SMILES string."
            )

    logger.info(
        "Validated %s (%s, MW=%.2f).",
        compound["name"],
        descriptors["formula"],
        descriptors["molecular_weight"],
    )


def main() -> None:
    project_config = load_yaml("config/query_config.yaml")
    compound_config = load_yaml("config/compounds.yaml")
    logger = setup_logging(
        project_config["paths"]["logs"],
        "chemical_structure_pipeline",
        "chemical_structure_pipeline.log",
    )

    output_dir = resolve_path(project_config["paths"]["results"]) / "chemical_structures"
    output_dir.mkdir(parents=True, exist_ok=True)
    processed_dir = resolve_path(project_config["paths"]["processed_data"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    logger.info("Starting chemical structure generation pipeline.")

    for compound in compound_config["compounds"]:
        name = compound["name"]
        smiles = compound["smiles"]
        logger.info("Processing compound: %s", name)

        molecule = smiles_to_molecule(smiles, name)
        descriptors = calculate_descriptors(molecule)
        validate_against_expected(compound, descriptors, logger)

        image_filename = name.lower().replace(" ", "_").replace("-", "_") + ".png"
        image_path = output_dir / image_filename
        Draw.MolToFile(molecule, str(image_path), size=(600, 400))

        summary_rows.append(
            {
                "name": name,
                "perturbation_label": compound["perturbation_label"],
                "class": compound["class"],
                "target": compound["target"],
                "pathway": compound["pathway"],
                "smiles": smiles,
                "structure_image": str(
                    Path("results") / "chemical_structures" / image_filename
                ),
                "pubchem_cid": compound.get("pubchem_cid", ""),
                "marker_genes": ";".join(compound.get("marker_genes", [])),
                **descriptors,
            }
        )
        logger.info("Saved structure image to %s", image_path)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = processed_dir / "compound_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    logger.info("Saved compound summary to %s", summary_path)
    logger.info("Chemical structure generation pipeline completed successfully.")


if __name__ == "__main__":
    main()