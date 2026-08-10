library(readr)
library(dplyr)
library(ggplot2)

umap_path <- "data/processed/umap_coordinates.csv"

umap_data <- read_csv(umap_path, show_col_types = FALSE)

umap_data <- umap_data %>%
  mutate(leiden = as.factor(leiden))

umap_plot <- ggplot(
  umap_data,
  aes(x = UMAP_1, y = UMAP_2, color = leiden)
) +
  geom_point(size = 0.6, alpha = 0.8) +
  labs(
    title = "Midway Showcase UMAP",
    subtitle = "Leiden clustering of CELLxGENE human liver subset",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Leiden Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "results/midway_showcase_umap.png",
  plot = umap_plot,
  width = 10,
  height = 7,
  dpi = 300
)

qc_overlay <- ggplot(
  umap_data,
  aes(x = UMAP_1, y = UMAP_2, color = total_counts)
) +
  geom_point(size = 0.6, alpha = 0.8) +
  labs(
    title = "UMAP Library Complexity Overlay",
    subtitle = "Total RNA counts per cell across the UMAP manifold",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Total counts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "results/midway_qc_overlay.png",
  plot = qc_overlay,
  width = 10,
  height = 7,
  dpi = 300
)

print("Saved results/midway_showcase_umap.png")
print("Saved results/midway_qc_overlay.png")