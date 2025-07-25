# Script: 03_celltype_composition_analysis.R
# Purpose: Visualize absolute and relative cell type distributions across tumor, PDX, and other sample types.

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# Step 1: Set identity class to cell type
Idents(harmony_obj) <- "celltype"

# Step 2: Extract metadata
meta_df <- harmony_obj@meta.data

# Step 3: Count number of cells per sample_type × celltype
cell_summary <- meta_df %>%
  group_by(sample_type, celltype) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(Proportion = Cell_Count / sum(Cell_Count)) %>%
  ungroup()

# Step 4A: Bar plot - Absolute cell counts per cell type per group
plot_counts <- ggplot(cell_summary, aes(x = sample_type, y = Cell_Count, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(
    title = "Absolute Cell Counts by Sample Type",
    x = "Sample Type",
    y = "Cell Count",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 4B: Bar plot - Proportional distribution
plot_proportions <- ggplot(cell_summary, aes(x = sample_type, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  labs(
    title = "Proportion of Cell Types by Sample Type",
    x = "Sample Type",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 5: Show plots side by side
plot_grid(plot_counts, plot_proportions, labels = c("A", "B"), ncol = 2)

# Step 6 (optional): Save combined figure
# ggsave("results/celltype_barplots_combined.png", width = 12, height = 6, dpi = 300)
