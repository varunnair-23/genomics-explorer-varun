"""Shared loading helpers for the chemical genomics pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


REQUIRED_OBS_COLUMNS = {
    "perturbation",
    "dose_value",
    "cell_line",
    "disease",
    "percent_mito",
    "ngenes",
}


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return data


def project_root() -> Path:
    # src/ is one level below the repo root
    return Path(__file__).resolve().parents[1]


def resolve_path(relative_path: str | Path) -> Path:
    path = Path(relative_path)
    if path.is_absolute():
        return path
    return project_root() / path


def validate_perturbation_adata(adata) -> None:
    missing = REQUIRED_OBS_COLUMNS - set(adata.obs.columns)
    if missing:
        raise ValueError(
            "AnnData is missing required perturbation metadata columns: "
            + ", ".join(sorted(missing))
        )
    if adata.n_obs == 0:
        raise ValueError("AnnData contains zero cells.")
    if adata.n_vars == 0:
        raise ValueError("AnnData contains zero genes.")


def setup_logging(log_dir: str | Path, logger_name: str, log_filename: str):
    import logging

    log_dir = resolve_path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] - %(message)s")

    file_handler = logging.FileHandler(log_dir / log_filename)
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger