import numpy as np
import pandas as pd
import h5py

csv_path = "merged_snRNA_counts.csv"
h5_path = "merged_snRNA_counts.h5"

print("Opening CSV header...")
# Read only header line
with open(csv_path, "r") as f:
    header = f.readline().strip().split(",")
genes = header[1:]   # first col = gene names

# Create H5 file
print("Creating H5…")
h5 = h5py.File(h5_path, "w")
h5.create_dataset("genes", data=np.array(genes, dtype="S"))

# Prepare count group
counts_grp = h5.create_group("counts")
barcodes = []

chunksize = 5000
reader = pd.read_csv(csv_path, chunksize=chunksize)

row_idx = 0
for chunk in reader:
    print("Processing chunk starting row:", row_idx)

    gene_names = chunk.iloc[:, 0].values.astype("S")
    barcodes.append(gene_names)

    count_vals = chunk.iloc[:, 1:].to_numpy(dtype=np.float32)

    for i in range(count_vals.shape[0]):
        counts_grp.create_dataset(str(row_idx), data=count_vals[i])
        row_idx += 1

# Save barcodes
print("Saving barcodes…")
all_barcodes = np.concatenate(barcodes)
h5.create_dataset("barcodes", data=all_barcodes)

h5.close()
print("DONE — saved merged_snRNA_counts.h5")