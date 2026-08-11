"""Streamlit dashboard linking compounds, perturbation UMAP, DEGs, and KEGG pathways."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

BASE_DIR = Path(__file__).parent

UMAP_TREATMENT = BASE_DIR / "results" / "umap_by_treatment.png"
UMAP_QC = BASE_DIR / "results" / "umap_qc_mito.png"
# fallbacks if you still have older plot names
UMAP_TREATMENT_FALLBACK = BASE_DIR / "results" / "midway_showcase_umap.png"
UMAP_QC_FALLBACK = BASE_DIR / "results" / "midway_qc_overlay.png"

PATHWAY_FILE = BASE_DIR / "data" / "processed" / "pathway_enrichment.csv"
DEG_FILE = BASE_DIR / "data" / "processed" / "drug_response_degs_significant.csv"
COMPOUND_FILE = BASE_DIR / "data" / "processed" / "compound_summary.csv"
UMAP_COORDS = BASE_DIR / "data" / "processed" / "umap_coordinates.csv"


st.set_page_config(page_title="Chemical-Genomics Explorer", layout="wide")


@st.cache_data
def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


pathways = load_csv(PATHWAY_FILE)
deg_sig = load_csv(DEG_FILE)
compounds = load_csv(COMPOUND_FILE)
umap_coords = load_csv(UMAP_COORDS)

st.title("Chemical-Genomics Explorer")
st.caption(
    "A549 lung adenocarcinoma (sci-Plex2): single-cell responses to Nutlin-3A, "
    "Vorinostat (SAHA), and BMS-345541 versus vehicle control."
)

st.markdown(
    """
    This dashboard connects **chemical structure**, **dose-aware single-cell manifolds**,
    **drug-vs-control differential expression**, and **KEGG pathway enrichment**.
    Clustering/UMAP are for visualization; pathway claims come from treated-vs-control DEGs.
    """
)

if compounds.empty:
    st.error("Compound summary missing. Run `python src/chem_utils.py`.")
    st.stop()

compound_names = compounds["name"].tolist()
selected_compound = st.sidebar.selectbox("Focus compound", compound_names)
compound_row = compounds[compounds["name"] == selected_compound].iloc[0]
perturbation_label = str(compound_row["perturbation_label"])

st.sidebar.markdown("---")
st.sidebar.markdown(
    f"**Target:** {compound_row['target']}  \n"
    f"**Class:** {compound_row['class']}  \n"
    f"**Expected pathway:** {compound_row['pathway']}"
)

tab1, tab2, tab3, tab4 = st.tabs(
    [
        "Compound",
        "Single-Cell Manifold",
        "Drug Response DEGs",
        "Pathway Enrichment",
    ]
)

with tab1:
    col1, col2 = st.columns([1, 2])
    with col1:
        image_path = BASE_DIR / compound_row["structure_image"]
        st.subheader(selected_compound)
        if image_path.exists():
            st.image(str(image_path), width="stretch")
        else:
            st.warning(f"Structure image not found: {image_path}")
    with col2:
        st.subheader("Chemical metadata")
        st.write(f"Perturbation label in dataset: `{perturbation_label}`")
        st.write(f"Class: {compound_row['class']}")
        st.write(f"Target: {compound_row['target']}")
        st.write(f"Pathway hypothesis: {compound_row['pathway']}")

        descriptor_columns = [
            "formula",
            "molecular_weight",
            "logp",
            "h_donors",
            "h_acceptors",
            "rotatable_bonds",
        ]
        descriptors = compound_row[
            [c for c in descriptor_columns if c in compounds.columns]
        ]
        st.dataframe(descriptors.to_frame(name="value"), width="stretch")

        st.subheader("SMILES")
        st.code(compound_row["smiles"])

        if not deg_sig.empty and "marker_genes" in compounds.columns:
            markers = [
                g.strip()
                for g in str(compound_row.get("marker_genes", "")).split(";")
                if g.strip()
            ]
            if markers:
                hit = deg_sig[
                    (deg_sig["compound"] == selected_compound)
                    & (deg_sig["names"].isin(markers))
                ].sort_values("pvals_adj")
                st.subheader("On-target marker genes in significant DEGs")
                if hit.empty:
                    st.info(
                        "Configured marker genes were not significant at current thresholds."
                    )
                else:
                    st.dataframe(
                        hit[["names", "logfoldchanges", "pvals_adj", "direction"]],
                        width="stretch",
                        hide_index=True,
                    )

with tab2:
    st.header("Perturbation manifold")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("UMAP by treatment")
        treatment_img = UMAP_TREATMENT if UMAP_TREATMENT.exists() else UMAP_TREATMENT_FALLBACK
        if treatment_img.exists():
            st.image(str(treatment_img), width="stretch")
        else:
            st.warning("Run plotting later to generate treatment UMAP images.")
    with col2:
        st.subheader("Mitochondrial QC overlay")
        qc_img = UMAP_QC if UMAP_QC.exists() else UMAP_QC_FALLBACK
        if qc_img.exists():
            st.image(str(qc_img), width="stretch")
        else:
            st.warning("Run plotting later to generate QC UMAP images.")

    st.info(
        "UMAP is a visualization layer. Neighborhood graph / Leiden clustering were "
        "computed in PCA space, not in 2D UMAP coordinates."
    )

if not umap_coords.empty and "treatment_group" in umap_coords.columns:
    st.subheader("Cells by treatment group")
    st.write(
        umap_coords["treatment_group"]
        .value_counts()
        .rename_axis("group")
        .reset_index(name="n_cells")
    )

with tab3:
    st.header(f"Differential expression: {selected_compound} vs control")
    if deg_sig.empty:
        st.warning(
            "Significant DEG table missing. Run `python src/differential_expression.py`."
        )
    else:
        compound_degs = deg_sig[deg_sig["compound"] == selected_compound].copy()
        if compound_degs.empty:
            st.warning(f"No significant DEGs found for {selected_compound}.")
        else:
            direction = st.radio("Direction", ["all", "up", "down"], horizontal=True)
            if direction != "all":
                compound_degs = compound_degs[compound_degs["direction"] == direction]
            compound_degs = compound_degs.sort_values("pvals_adj", ascending=True)
            st.metric("Significant DEGs", len(compound_degs))
            show_cols = [
                c
                for c in ["names", "logfoldchanges", "pvals_adj", "scores", "direction"]
                if c in compound_degs.columns
            ]
            st.dataframe(
                compound_degs[show_cols].head(50),
                width="stretch",
                hide_index=True,
            )

with tab4:
    st.header(f"KEGG pathways enriched by {selected_compound} response genes")
    if pathways.empty:
        st.warning("Pathway file missing. Run `python src/pathway_analysis.py`.")
    else:
        compound_pathways = pathways[pathways["compound"] == selected_compound].copy()
        if compound_pathways.empty:
            st.warning(f"No pathway results for {selected_compound}.")
        else:
            direction = st.radio(
                "DEG direction used for enrichment",
                sorted(compound_pathways["direction"].astype(str).unique()),
                horizontal=True,
                key="pathway_direction",
            )
            view = compound_pathways[
                compound_pathways["direction"].astype(str) == direction
            ].copy()
            view = view[view["adjusted_p_value"] < 0.05].sort_values(
                "adjusted_p_value", ascending=True
            )
            if view.empty:
                st.info("No pathways at adjusted p < 0.05 for this direction.")
            else:
                show_cols = [
                    c
                    for c in [
                        "pathway",
                        "adjusted_p_value",
                        "combined_score",
                        "overlapping_genes",
                    ]
                    if c in view.columns
                ]
                st.dataframe(
                    view[show_cols].head(20),
                    width="stretch",
                    hide_index=True,
                )
                st.caption(
                    "Enrichment is computed on significant drug-vs-control DEGs, "
                    "not on unsupervised cluster markers."
                )