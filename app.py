from pathlib import Path

import pandas as pd
import streamlit as st


BASE_DIR = Path(__file__).parent

UMAP_IMAGE = BASE_DIR / "results" / "midway_showcase_umap.png"
QC_IMAGE = BASE_DIR / "results" / "midway_qc_overlay.png"
PATHWAY_FILE = BASE_DIR / "data" / "processed" / "pathway_enrichment.csv"
COMPOUND_FILE = BASE_DIR / "data" / "processed" / "compound_summary.csv"
STRUCTURE_DIR = BASE_DIR / "results" / "chemical_structures"


st.set_page_config(
    page_title="Chemical-Genomics Explorer",
    layout="wide"
)


@st.cache_data
def load_csv(path):
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


pathways = load_csv(PATHWAY_FILE)
compounds = load_csv(COMPOUND_FILE)


st.title("Chemical-Genomics Explorer")
st.caption("Single-cell clustering, pathway enrichment, and small-molecule structure visualization")

st.markdown(
    """
    This dashboard connects single-cell RNA-seq analysis with chemical biology.
    The pipeline streams public CELLxGENE data, performs QC, PCA, Leiden clustering,
    UMAP visualization, marker-gene ranking, KEGG pathway enrichment, and RDKit-based
    molecular structure rendering.
    """
)

tab1, tab2, tab3 = st.tabs([
    "Single-Cell Manifold",
    "Pathway Enrichment",
    "Chemical Compounds"
])

with tab1:
    st.header("Single-Cell UMAP Showcase")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Leiden Cluster UMAP")
        if UMAP_IMAGE.exists():
            st.image(str(UMAP_IMAGE), width="stretch")
        else:
            st.warning("UMAP image not found. Run src/plot_clusters.R first.")

    with col2:
        st.subheader("Library Complexity QC Overlay")
        if QC_IMAGE.exists():
            st.image(str(QC_IMAGE), width="stretch")
        else:
            st.warning("QC overlay image not found. Run src/plot_clusters.R first.")

    st.info(
        "UMAP is used here as a visualization layer. Clustering was performed on the PCA-based "
        "nearest-neighbor graph, not directly on the 2D UMAP coordinates."
    )

with tab2:
    st.header("KEGG Pathway Enrichment")

    if pathways.empty:
        st.warning("Pathway enrichment file not found. Run src/pathway_analysis.py first.")
    else:
        cluster_options = sorted(
        pathways["cluster"].astype(str).unique(),
        key=lambda x: int(x)
        )
        selected_cluster = st.selectbox(
            "Select Leiden cluster",
            cluster_options
        )

        cluster_pathways = pathways[
            pathways["cluster"].astype(str) == selected_cluster
        ].copy()

        cluster_pathways = cluster_pathways.sort_values(
            "adjusted_p_value",
            ascending=True
        )

        st.subheader(f"Top pathways for cluster {selected_cluster}")

        display_columns = [
            "pathway",
            "adjusted_p_value",
            "combined_score",
            "overlapping_genes"
        ]

        available_columns = [
            column for column in display_columns
            if column in cluster_pathways.columns
        ]

        st.dataframe(
            cluster_pathways[available_columns].head(15),
            width="stretch",
            hide_index=True
        )

with tab3:
    st.header("Small-Molecule Chemical Structures")

    if compounds.empty:
        st.warning("Compound summary file not found. Run src/chem_utils.py first.")
    else:
        compound_names = compounds["name"].tolist()

        selected_compound = st.selectbox(
            "Select compound",
            compound_names
        )

        compound_row = compounds[
            compounds["name"] == selected_compound
        ].iloc[0]

        col1, col2 = st.columns([1, 2])

        with col1:
            image_path = BASE_DIR / compound_row["structure_image"]

            st.subheader(selected_compound)

            if image_path.exists():
                st.image(str(image_path), width="stretch")
            else:
                st.warning(f"Structure image not found: {image_path}")

        with col2:
            st.subheader("Compound Metadata")

            st.write("Class:", compound_row["class"])
            st.write("Target:", compound_row["target"])
            st.write("Pathway:", compound_row["pathway"])

            descriptor_columns = [
                "molecular_weight",
                "logp",
                "h_donors",
                "h_acceptors",
                "rotatable_bonds"
            ]

            descriptors = compound_row[
                [col for col in descriptor_columns if col in compounds.columns]
            ]

            st.dataframe(
                descriptors.to_frame(name="value"),
                width = "stretch"
            )

            st.subheader("SMILES")
            st.code(compound_row["smiles"])