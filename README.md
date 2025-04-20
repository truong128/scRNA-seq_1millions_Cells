👏 Using R (BPCells with Seurat Objects) for single-cell RNA-seq (sc-RNA-seq)
https://lnkd.in/gZhR9YWD

💯 💥 The entire test (1M neurons) was completed in just 2 minutes on my MacBook Pro with 64GB of RAM!!!

A case study on using BPCells with Seurat Objects, as presented in the Seurat v5 vignette (compiled October 20, 2023):

💡 Background & Motivation
Challenge before: Processing large scRNA-seq datasets (e.g., 1M+ cells) was memory-intensive and often infeasible on standard machines.

👉 Now with BPCells: Efficient storage + on-disk operations enable scalable analysis, even on large datasets.

🤜 💥 Key Tools
BPCells: Enables fast, low-memory single-cell analysis using bit-packing and C++-based caching.
Seurat: Widely used R toolkit for single-cell genomics.
SeuratObject, SeuratDisk, Azimuth: Support seamless interaction between Seurat objects, on-disk storage, and reference mapping.
------------------------------snapshot of script------------------------------
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)

# Load 10X HDF5 file
brain.data <- open_matrix_10x_hdf5(
 path = "/Users/truong/Downloads/1M_neurons_filtered_gene_bc_matrices_h5.h5"
)

output_dir <- "/Users/truong/Documents/6_Single_Cells_v2"

if (!dir.exists(output_dir)) {
 write_matrix_dir(
 mat = brain.data,
 dir = output_dir
 )
} else {
 message("Directory already exists. Skipping write_matrix_dir().")
}

brain.mat <- open_matrix_dir(dir = output_dir)

brain.mat <- Azimuth:::ConvertEnsembleToSymbol(
 mat = brain.mat,
 species = "mouse"
)

brain <- CreateSeuratObject(counts = brain.mat)

print(brain)

VlnPlot(brain, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "counts", alpha = 0.1)

brain <- NormalizeData(brain, normalization.method = "LogNormalize")
VlnPlot(brain, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "data", alpha = 0.1)
-------------------------------------------------------------------------------
