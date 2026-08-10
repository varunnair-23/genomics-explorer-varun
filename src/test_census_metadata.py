import cellxgene_census


def main():
    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        obs = (
            census["census_data"]["homo_sapiens"]
            .obs
            .read(
                value_filter=(
                    "tissue_general == 'liver' "
                    "and  disease == 'normal' "
                    "and assay == \"10x 3' v3\" "
                    "and is_primary_data == True"
                ),
                column_names=[
                    "soma_joinid",
                    "cell_type",
                    "tissue_general",
                    "tissue",
                    "disease",
                    "assay",
                    "is_primary_data"
                ],
            )
            .concat()
            .to_pandas()
        )

    print(obs.head())
    print("Matched cells:", obs.shape[0])
    print()
    print("Assays found:")
    print(obs["assay"].value_counts().head(20))

if __name__ == "__main__":
    main()