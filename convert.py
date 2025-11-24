import pandas as pd
import numpy as np
import h5py

csv_path = "merged_snRNA_counts.csv"
h5_path  = "merged_snRNA_counts.h5"

chunksize = 10000  # adjust if needed

first = True
genes = None

with h5py.File(h5_path, "w") as h5:

    counts_grp = h5.create_group("counts")

    row_idx = 0

    for chunk in pd.read_csv(csv_path, sep=",", chunksize=chunksize):
        print("Loaded chunk:", row_idx)

        if first:
            genes = chunk.columns[1:]
            h5.create_dataset("genes", data=np.array(genes, dtype="S"))
            barcodes_list = []
            first = False

        barcodes = chunk.iloc[:, 0].values.astype("S")
        barcodes_list.append(barcodes)

        # store chunk row-by-row
        for i in range(chunk.shape[0]):
            row = chunk.iloc[i, 1:].values.astype(np.float32)
            counts_grp.create_dataset(str(row_idx), data=row)
            row_idx += 1

    # store barcodes
    h5.create_dataset("barcodes", data=np.array(np.concatenate(barcodes_list), dtype="S"))

print("Done. Saved:", h5_path)