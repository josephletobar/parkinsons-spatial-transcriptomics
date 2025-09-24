import pandas as pd

# path to dataset
file_path = "GSE253975_series_matrix.txt"

# read as a dataframe
df = pd.read_csv(file_path, sep="\t", comment="!", index_col=0)

print(df.shape)      # rows x columns
print(df.head())     # peek at the first few rows