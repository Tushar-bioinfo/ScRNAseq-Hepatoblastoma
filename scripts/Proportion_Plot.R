library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)  # optional for combined plots

# Set identity to cell type
Idents(harmony_obj) <- "celltype"

# Extract metadata
meta_df <- harmony_obj@meta.data

# Summarize cell counts and proportions
cell_summary <- meta_df %>%
  group_by(sample_type, celltype) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(Proportion = Cell_Count / sum(Cell_Count)) %>%
  ungroup()

# Absolute cell count plot
plot_counts <- ggplot(cell_summary, aes(x = sample_type, y = Cell_Count, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(title = "Absolute Cell Counts by Sample Type",
       x = "Sample Type", y = "Cell Count", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Proportion plot
plot_proportions <- ggplot(cell_summary, aes(x = sample_type, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  labs(title = "Proportion of Cell Types by Sample Type",
       x = "Sample Type", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show plots side by side
plot_grid(plot_counts, plot_proportions, labels = c("A", "B"), ncol = 2)
