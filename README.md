# Chemical-Genomics Explorer

The app is available here:
https://genomics-explorer-varun-amec2r67g5tdrqnhrmcamj.streamlit.app/

Single-cell chemical genomics for A549 lung adenocarcinoma cells treated with Nutlin-3A, Vorinostat (SAHA), and BMS-345541 versus vehicle control (sci-Plex2).

## Run the dashboard

python3 -m pip install -r requirements.txt
python3 -m streamlit run app.py

The app reads the committed CSV and PNG files. It does not need the raw .h5ad.

## Main finding

Nutlin-3A, an MDM2 inhibitor, increases canonical p53-response genes including CDKN1A, MDM2, and BAX. Pathway enrichment of those drug-versus-control genes recovers p53 signaling. UMAP is only a visualization; the drug claims come from the differential-expression tables.
