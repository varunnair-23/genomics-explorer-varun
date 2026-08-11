"""Download the sci-Plex2 AnnData (.h5ad) file into data/raw/."""

from __future__ import annotations

import logging
import urllib.request
from pathlib import Path

import yaml


def load_config(config_path: str = "config/query_config.yaml") -> dict:
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with path.open("r") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return data


def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("download_data")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] - %(message)s")
    file_handler = logging.FileHandler(log_dir / "download_data.log")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def main() -> None:
    config = load_config()
    logger = setup_logging(Path(config["paths"]["logs"]))

    raw_dir = Path(config["paths"]["raw_data"])
    raw_dir.mkdir(parents=True, exist_ok=True)

    filename = config["dataset"]["filename"]
    url = config["dataset"]["download_url"]
    destination = raw_dir / filename

    # Skip re-download if a reasonably large file already exists (~145 MB expected)
    if destination.exists() and destination.stat().st_size > 1_000_000:
        logger.info("Dataset already present at %s", destination)
        logger.info("File size: %s bytes", destination.stat().st_size)
        return

    logger.info("Downloading sci-Plex2 dataset")
    logger.info("URL: %s", url)
    logger.info("Saving to: %s", destination)

    try:
        urllib.request.urlretrieve(url, destination)
    except Exception as exc:
        raise RuntimeError(f"Download failed: {exc}") from exc

    size = destination.stat().st_size
    if size < 1_000_000:
        raise RuntimeError(
            f"Downloaded file looks too small ({size} bytes). Check the URL."
        )

    logger.info("Download complete (%s bytes).", size)


if __name__ == "__main__":
    main()