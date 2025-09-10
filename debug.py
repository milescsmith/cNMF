from pathlib import Path

import numpy as np
import scanpy as sc
from numpy import random

from cnmf import cNMF

random.seed(318)

# adata = sc.read_h5ad("20241205_refined_bcells_with_nmf.h5ad")
adata = sc.read_h5ad("tcell_counts.h5ad")

adata.var_names_make_unique()
numiter=20 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
numhvgenes=2000 ## Number of over-dispersed genes to use for running the actual factorizations

project_dir = Path.home().joinpath("workspace", "cNMF")
output_directory = project_dir / ("nmf")
run_name = "tcells"
if not output_directory.joinpath(run_name).exists():
    output_directory.joinpath(run_name).mkdir(parents=True)

print('Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10"')
K = " ".join([str(i) for i in range(5,11)])

print("Path to the filtered counts dataset we output previously")
countfn = project_dir / "tcell_counts.h5ad"
cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
cnmf_obj.prepare(counts_fn=countfn, components=np.arange(5,11), n_iter=numiter, seed=318, num_highvar_genes=numhvgenes)

# cnmf_obj.factorize()

print("combining")
cnmf_obj.combine()
print("selection_plot")
cnmf_obj.k_selection_plot(close_fig=False)
selected_k = 8
density_threshold = 0.05
print("reading tmp")
tpm = sc.read(cnmf_obj.paths["tpm"])

print("consensus")
cnmf_obj.consensus(
    k=selected_k,
    density_threshold=density_threshold,
    show_clustering=True,
    close_clustergram_fig=False
)
print("reading")
adata = sc.read(countfn)
print("to gpu")

print("Obtain high variance genes that were used for cNMF as these were saved to a text file")
hvgs = open("./trial/bcells/bcells.overdispersed_genes.txt").read().split("\n")

print("TPT normalization")
sc.pp.normalize_total(adata, target_sum=10**4) ## TPT normalization

print(
    "Set log-normalized data to the raw attribute of the AnnData object to make it easy to plot expression levels of individual genes. "
    "This does not log normalize the actual AnnData data matrix"
)
adata.raw = sc.pp.log1p(adata.copy(), copy=True)
## Subset out only the high-variance genes

adata = adata[:,hvgs]

print("scaling")
sc.pp.scale(adata)

print("pca")
sc.pp.pca(adata)

sc.pl.pca_variance_ratio(adata, log=True)
## Construct the nearest neighbor graph for UMAP

print("neighbors")
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=15)
## Run UMAP

print("umap")
sc.tl.umap(adata)
aifi_bcell_genes = [
    "AIM2",
    "CD9",
    "CD19",
    "CD24",
    "CD27",
    "CD38",
    "PAX5",
    "FCER2",
    "IGHG1",
    "IGHG2",
    "IGHG3",
    "IGHG4",
    "IGHA1",
    "IGHA2",
    "IGHE",
    "IL4R",
    "ITGAX",
    "MME",
    "MS4A1",
    "MZB1",
    "FCRL4",
    "FCRL5",
    "TBX21",
    "PDCD1",
    "PAX5",
    "PRDM1",
    "LAMF7",
    "XBP1",
    "ZEB2",
]
## Plot the UMAP with some cannonical marker genes to see that the apparent grouping makes sense

sc.pl.umap(adata, color=adata.var_names.intersection(aifi_bcell_genes), size=15,
           use_raw=True, ncols=3, alpha=0.75)
usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_k, density_threshold=density_threshold)
usage_norm.columns = [f"Usage_{i}" for i in usage_norm.columns]
