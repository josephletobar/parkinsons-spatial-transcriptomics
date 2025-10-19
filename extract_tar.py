import pandas as pd
import tarfile

# path to tar file
tar_path = "GSE253975_RAW.tar"

# open and inspect
with tarfile.open(tar_path, "r") as tar:
    print(tar.getnames())  # lists all contained files

# extract contained files
with tarfile.open(tar_path, "r") as tar:
    tar.extractall("GSE253975_data")  # creates a new folder with the contents