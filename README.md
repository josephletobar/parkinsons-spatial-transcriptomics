# Analysis of Parkinson’s Disease Spatial Transcriptomics Dataset

## Data

### `GSE253975_data/`
Raw data directory containing the downloaded files from  
[GEO: GSE253975](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253975).

### `ckpt/`
Folder containing pretrained CellPLM checkpoint weights:  
[CellPLM Checkpoints (Dropbox)](https://www.dropbox.com/scl/fo/i5rmxgtqzg7iykt2e9uqm/h?rlkey=o8hi0xads9ol07o48jdityzv1&e=1&dl=0).

### `adata_preprocessed.h5ad`
Preprocessed AnnData object saved after running the preprocessing notebook (loading/merging counts, filtering, and CP10K + log-normalization).

## Code

### `preprocessing.ipynb`
Loads the raw Visium .txt.gz count matrices in GSE253975_data, inspects them, builds a combined AnnData with sample_id/batch metadata, filters zero-count genes, normalizes to CP10K + log1p, and saves the reusable adata_preprocessed.h5ad (with a small matrix preview at the end).

### `parkinsons_analysis.ipynb`
Starts from adata_preprocessed.h5ad, selects 2k HVGs, scales, runs PCA, corrects batch effects with Harmony, builds neighbor graph + UMAP, and clusters with K-Means. It computes silhouette/compactness/batch-mixing metrics for PCA vs Harmony.

### `parkinsons_analysis_cellplm.ipynb`
Also loads adata_preprocessed.h5ad but replaces PCA with 512-dim CellPLM embeddings (CPU fallback). Uses those embeddings for neighbors, UMAP, K-Means clustering, and silhouette/compactness/batch-mixing metrics.

## How to Run
1. Install the required Python packages listed in `requirements.txt`.
2. Make sure the `GSE253975_data/` directory and the `ckpt/` checkpoint folder are in place.
3. Run `preprocessing.ipynb` to load the GEO files and generate the unified normalized `adata_preprocessed.h5ad`.
4. Run `parkinsons_analysis.ipynb` and `parkinsons_analysis_cellplm.ipynb` to perform their respective analyses and automatically generate figures in the `figures/` folder.
