###############################################
# ONE-SHOT RCTD PIPELINE (H5 + H5AD VERSION)
# Author: Joseph's final single-file pipeline
# Requirements: merged_snRNA_counts.h5, processed_visium.h5ad
###############################################

library(hdf5r)
library(Matrix)
library(spacexr)
library(zellkonverter)

cat("\n=== 1) Loading snRNA HDF5 Reference ===\n")

# ---- LOAD H5 REFERENCE ----
h5 <- H5File$new("merged_snRNA_counts.h5", mode = "r")

genes    <- h5[["genes"]]$read()
barcodes <- h5[["barcodes"]]$read()

# LOAD ALL COUNT ROWS
count_keys <- names(h5[["counts"]])
count_list <- vector("list", length(count_keys))

for (i in seq_along(count_keys)) {
  count_list[[i]] <- h5[["counts"]][[count_keys[i]]]$read()
}

h5$close()

# BUILD MATRIX
mat <- do.call(rbind, count_list)
rownames(mat) <- genes
colnames(mat) <- barcodes

ref_mat <- Matrix(mat, sparse = TRUE)

cat("snRNA reference matrix:", nrow(ref_mat), "genes x", ncol(ref_mat), "cells\n")

# DUMMY CELL TYPES (you never had annotations)
cell_types <- rep("Unknown", ncol(ref_mat))

reference <- Reference(ref_mat, cell_types)

cat("Reference object created.\n")


###############################################
# 2) LOAD VISIUM ANNDATA
###############################################

cat("\n=== 2) Loading Visium processed H5AD ===\n")

adata <- readH5AD("processed_visium.h5ad")

sp_counts <- t(as.matrix(adata$X))
coords <- as.data.frame(adata$obsm$spatial)
rownames(coords) <- rownames(sp_counts)

cat("Visium matrix:", nrow(sp_counts), "spots x", ncol(sp_counts), "genes\n")

spatialRNA <- SpatialRNA(coords, sp_counts)
cat("SpatialRNA object prepared.\n")


###############################################
# 3) RUN RCTD
###############################################

cat("\n=== 3) Running RCTD ===\n")
rctd <- create.RCTD(spatialRNA, reference, max_cores = 4)
rctd <- run.RCTD(rctd)

cat("\nRCTD finished successfully.\n")


###############################################
# 4) SAVE OUTPUTS
###############################################

cat("\n=== 4) Saving outputs ===\n")

results <- get_deconvolved_results(rctd)

write.csv(results, "rctd_results.csv", row.names = TRUE)
saveRDS(rctd, "rctd_output.rds")

cat("\n=== DONE ===\n")
cat("Saved:\n")
cat(" - rctd_results.csv\n")
cat(" - rctd_output.rds\n")
###############################################