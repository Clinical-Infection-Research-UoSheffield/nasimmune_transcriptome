library(BayesPrism)

blood_deconvole <- readRDS("./output/processed_data/blood_bayesprism_object.rds")

# 1. Filter out cell types with low fraction and high cv

# 2. Describe average composition at day 0 and day 2.

# 3. Look at cell type specific gene expression. Cluster vst() normalised by day or immune response.