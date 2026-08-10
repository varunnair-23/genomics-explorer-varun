from pathlib import Path
import numpy as np
import pandas as pd
from anndata import AnnData

def main():
    output_dir = Path("data/processed") 
    output_dir.mkdir(parents=True, exist_ok=True)

    X = np.array([
        [5, 0, 2],
        [3, 1, 0],
        [0, 4, 6],
        [8, 0, 1]
    ])

    obs = pd.DataFrame({
        "cell_id": ["cell1", "cell2", "cell3", "cell4"],
        "cell_type": ["T cell", "T cell", "B cell", "Monocyte"],
        "tissue": ["lung", "lung", "blood", "blood"]
    }).set_index("cell_id")

    var = pd.DataFrame({
        "gene_symbol": ["CD3D", "CD3E", "CD19"],
        "gene_id": ["ENSG00000123456", "ENSG00000123457", "ENSG00000123458"]
    }).set_index("gene_symbol")

    adata = AnnData(X=X, obs=obs, var=var)

    print(adata)
    print()
    print("adata.X expression matirx:")
    print(adata.X)
    print()
    print("adata.obs cell metadata:")
    print(adata.obs)
    print()
    print("adata.var gene metadata:")
    print(adata.var)

    adata.write(output_dir / "anndata_practice.h5ad")

    print()
    print("Anndata object saved to data/processed/anndata_practice.h5ad")

if __name__ == "__main__":
    main()


