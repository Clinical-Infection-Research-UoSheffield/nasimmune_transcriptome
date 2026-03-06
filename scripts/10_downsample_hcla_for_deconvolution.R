library(dplyr)
library(cellxgene.census)
library(tiledbsoma)
library(Matrix)

set.seed(42)

#Speed up download speed
options(tiledbsoma.init_buffer_bytes = 4 * 1024^3)

# get Human Cell Atlas reference
census <- open_soma(census_version = "2025-11-08")

## Look at metadata from HCA
obs_df <- as.data.frame(
  census$get("census_data")$get("homo_sapiens")$obs$read(
    value_filter = "tissue == 'blood' & is_primary_data == TRUE & disease == 'normal' & suspension_type == 'cell'",
    column_names = c("cell_type", "soma_joinid", "suspension_type")
  )$concat()
)

n_cells <- nrow(obs_df)
print(paste("Total cells found:", n_cells))

cell_freq <- as.data.frame(table(obs_df$cell_type)) %>%
    filter(Freq != 0) %>%
    arrange(desc(Freq))

## Filter for specific cell types
cell_types <- c(
    "erythrocyte",
    "platelet",
    "neutrophil",
    "classical monocyte",
    "intermediate monocyte",
    "non-classical monocyte",
    "naive thymus-derived CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell",
    "central memory CD4-positive, alpha-beta T cell",
    "central memory CD8-positive, alpha-beta T cell",
    "effector memory CD4-positive, alpha-beta T cell",
    "effector memory CD8-positive, alpha-beta T cell",
    "T follicular helper cell",
    "regulatory T cell",
    "mucosal invariant T cell",
    "gamma-delta T cell",
    "natural killer cell",
    "naive B cell",
    "memory B cell",
    "transitional stage B cell",
    "plasmablast",
    "plasmacytoid dendritic cell",
    "myeloid dendritic cell",
    "mast cell"
)

obs_df <- filter(obs_df, cell_type %in% cell_types)

## Take a random sample of 10,000 cells per type. If less than 10,000, include all cells 
sampled_metadata <- obs_df %>%
  group_by(cell_type) %>%
  sample_n(size = min(n(), 1000)) %>% 
  ungroup()

sampled_cell_freq <- as.data.frame(table(sampled_metadata$cell_type)) %>%
  filter(Freq != 0) %>%
  arrange(desc(Freq))

sampled_ids <- sampled_ids <- sampled_metadata$soma_joinid

# Pull data from HCA
sig_mat_seurat <- get_seurat(
  census,
  organism = "Homo sapiens",
  obs_coords = sampled_ids,
  obs_column_names = c("cell_type", "assay", "suspension_type")
)

saveRDS(sig_mat_seurat, "./output/processed_data/blood_hca_for_sig_mat.rds")
