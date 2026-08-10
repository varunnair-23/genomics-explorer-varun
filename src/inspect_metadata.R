library(readr)
library(dplyr)
library(ggplot2)

metadata_path <- "data/processed/cell_metadata_qc.csv"
output_path <- "results/library_distributions.png"

cell_metadata <- read_csv(metadata_path)

print(head(cell_metadata))
print(colnames(cell_metadata))

cell_type_counts <- cell_metadata %>%
  count(cell_type, sort = TRUE)

print(cell_type_counts)

top_cell_types <- cell_type_counts %>%
  slice_head(n = 10) %>%
  pull(cell_type)

plot_data <- cell_metadata %>%
  filter(cell_type %in% top_cell_types)

png(
  filename = output_path,
  width = 1200,
  height = 800
)

ggplot(plot_data, aes(x = reorder(cell_type, total_counts, median), y = total_counts)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    title = "Library Size Distribution Across Liver Cell Types",
    subtitle = "CELLxGENE Census subset after Scanpy QC processing",
    x = "Cell type",
    y = "Total RNA counts per cell"
  ) +
  theme_minimal(base_size = 14)

dev.off()

print(paste("Saved plot to", output_path))