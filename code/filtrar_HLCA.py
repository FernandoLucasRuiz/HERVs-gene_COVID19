from pathlib import Path
import gc

import numpy as np
from scipy import sparse
import scanpy as sc


BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "dbb5ad81-1713-4aee-8257-396fbabe7c6e.h5ad"

MIN_CELLS = 20
CHUNK_SIZE = 1000  


def count_cells_per_gene(adata, obs_mask, chunk_size=CHUNK_SIZE):

    obs_mask = np.asarray(obs_mask, dtype=bool)
    selected_total = int(obs_mask.sum())

    cells_per_gene = np.zeros(adata.n_vars, dtype=np.int64)
    processed = 0

    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        local_mask = obs_mask[start:stop]

        if not local_mask.any():
            continue

        X_chunk = adata.X[start:stop, :]
        X_chunk = X_chunk[local_mask, :]

        if sparse.issparse(X_chunk):
            cells_per_gene += np.asarray(
                X_chunk.getnnz(axis=0)
            ).ravel()
        else:
            cells_per_gene += np.count_nonzero(
                X_chunk,
                axis=0
            )

        processed += int(local_mask.sum())

        del X_chunk
        gc.collect()

    return cells_per_gene


def export_filtered_subset(label, mask_builder, output_name):
    print(f"\n=== {label} ===", flush=True)

    adata = sc.read_h5ad(INPUT_FILE, backed="r")

    obs_mask = np.asarray(mask_builder(adata.obs), dtype=bool)
    n_selected = int(obs_mask.sum())

    cells_per_gene = count_cells_per_gene(adata, obs_mask)
    gene_mask = cells_per_gene >= MIN_CELLS

    output_file = BASE_DIR / output_name
    if output_file.exists():
        output_file.unlink()

    adata[obs_mask, gene_mask].copy(filename=output_file)

    adata.file.close()
    del adata
    gc.collect()


export_filtered_subset(
    label="Células inmunes: pulmón, normal/COVID-19",
    mask_builder=lambda obs: (
        obs["tissue"].isin(["lung", "lung parenchyma"])
        & obs["disease"].isin(["normal", "COVID-19"])
        & (obs["ann_level_1"] == "Immune")
    ),
    output_name="HLCA_lung_immune_normal_COVID_filtered.h5ad",
)

export_filtered_subset(
    label="Fibroblastos: pulmón, normal/COVID-19",
    mask_builder=lambda obs: (
        obs["tissue"].isin(["lung"])
        & obs["disease"].isin(["normal", "COVID-19"])
        & (obs["ann_level_2"] == "Fibroblast lineage")
    ),
    output_name="HLCA_lung_Fibroblast_normal_COVID_filtered.h5ad",
)
