library(BayesPrism)
library(dplyr)
library(Seurat)
library(Matrix)
library(org.Hs.eg.db)

#Read and clean bulk data (convert genesymbols to ENSEMBL IDs)
Bld17 <- read.delim("./data/2017_blood_counts.tsv")
Bld18 <- read.delim("./data/2018_blood_counts.tsv")
bulk <- bind_cols(Bld17, Bld18)
bulk <- t(bulk)

current_symbols <- colnames(bulk)

# Map Symbols to Ensembl IDs
mapped_ids <- mapIds(org.Hs.eg.db,
                     keys = current_symbols,
                     column = "ENSEMBL",
                     keytype = "SYMBOL",
                     multiVals = "first")

valid_genes <- !is.na(mapped_ids)
bulk_mapped <- bulk[, valid_genes]
new_ensembl_ids <- mapped_ids[valid_genes]
colnames(bulk_mapped) <- new_ensembl_ids
bulk <- bulk_mapped
rm(bulk_mapped)
bulk <- bulk[, !duplicated(colnames(bulk))]

# Read and clean single cell data
sig_mat_seurat <- readRDS("./output/processed_data/blood_hca_for_sig_mat.rds")
raw_counts <- GetAssayData(sig_mat_seurat, assay = "RNA", slot = "counts")

sc <- t(raw_counts)

cell_type_labels <- sig_mat_seurat$cell_type
cell_state_labels <- sig_mat_seurat$cell_type

#QC cell type labels
pdf("./output/figs/bayesprism_cell_label_qc.pdf", width = 12, height = 12)
plot.cor.phi(sc, 
             input.labels = cell_type_labels,
             margins=c(20,20),
            )
dev.off()

pdf("./output/figs/bayesprism_outlier_genes.pdf", width = 8, height = 8)
plot.scRNA.outlier(sc,
                   cell.type.labels = cell_type_labels,
                   species = "hs")
dev.off()


pdf("./output/figs/bayesprism_bulk_outlier_genes.pdf", width = 8, height = 8)
plot.bulk.outlier(bulk.input = bulk,
                  sc.input = sc,
                  cell.type.labels = cell_type_labels,
                  species = "hs")
dev.off()

# Cleanup genes from scRNA-seq
sc_filtered <- cleanup.genes (input=sc,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

pdf("./output/figs/bayesprism_compare_sc_bulk.pdf", width = 8, height = 8)
plot.bulk.vs.sc(sc.input = sc_filtered,
                bulk.input = bulk)
dev.off()

# Filter to just include protein coding genes
sc_protein_coding <- select.gene.type (sc_filtered,
                                        gene.type = "protein_coding")

#Run deconvolution
bloodPrism <- new.prism(
  reference=sc_protein_coding,
  mixture=bulk,
  input.type="count.matrix",
  cell.type.labels = cell_type_labels,
  cell.state.labels = cell_state_labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

blood_deconvole <- run.prism(prism = bloodPrism, n.cores = 64)
saveRDS(blood_deconvole, "./output/processed_data/blood_bayesprism_object.rds")

