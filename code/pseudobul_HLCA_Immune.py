import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

file_in = "HLCA_lung_immune_normal_COVID_filtered.h5ad"

sample_col = "donor_id"
disease_col = "disease"

adata = sc.read_h5ad(file_in)

if adata.raw is not None:
    X = adata.raw.X
    genes = adata.raw.var_names.astype(str)
else:
    X = adata.X
    genes = adata.var_names.astype(str)

obs = adata.obs.copy()

keep = (
    obs[sample_col].notna() &
    obs[disease_col].isin(["normal", "COVID-19"])
)

obs = obs.loc[keep].copy()
X = X[keep.values, :]

if not sparse.issparse(X):
    X = sparse.csr_matrix(X)
else:
    X = X.tocsr()

samples = obs[sample_col].astype(str)
sample_levels = pd.Index(sorted(samples.unique()))

codes = sample_levels.get_indexer(samples)

G = sparse.csr_matrix(
    (
        np.ones(len(codes)),
        (np.arange(len(codes)), codes)
    ),
    shape=(len(codes), len(sample_levels))
)

pb = X.T @ G

pb_df = pd.DataFrame(
    pb.toarray(),
    index=genes,
    columns=sample_levels
)

pb_df.index.name = "gene_id"

meta = obs[[sample_col, disease_col]].drop_duplicates()

check = meta.groupby(sample_col)[disease_col].nunique()
bad = check[check > 1]

if len(bad) > 0:
    raise ValueError("Hay donors con más de una disease: " + ",".join(bad.index.astype(str)))

meta = meta.set_index(sample_col).loc[sample_levels]
meta.index.name = "sample_id_pb"

n_cells = obs.groupby(sample_col).size()
meta["n_cells"] = n_cells.loc[sample_levels].values

pb_df.to_csv("HLCA_immune_pseudobulk_counts_by_donor.csv.gz")
meta.to_csv("HLCA_immune_pseudobulk_metadata_by_donor.csv")

print(pb_df.shape)
print(meta[disease_col].value_counts())
print(meta.head())
