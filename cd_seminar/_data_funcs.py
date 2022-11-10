
import numpy as np
import anndata as a
import web_files as wf
import vinplots
import os

from ._plot_umaps import overview_plot

# -- define local paths to data / tools: -------------------------------------------------

PathDict = {
    "h5ad": "./content/adata.h5ad",
    "umap": "./content/umap.pkl",
    "pca":  "./content/pca.pkl",
    "ckpt": "./content/ckpt.pt",
}



# -- download tutorial content: ----------------------------------------------------------
def download_tutorial_content(data_dir="./content"):

    base_address = "https://figshare.com/ndownloader/files/{}"
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    figshare_files = {
        "adata.h5ad": 38171943,
        "umap.pkl": 38175252,  # modified to protocol=4
        "pca.pkl" : 38175255,  # modified to protocol=4
        "ckpt.pt" : 38171649,
    }

    for filename, http_key in figshare_files.items():
        print("Downloading {:.<10}...".format(filename), end="")
        f = wf.WebFile(
            http_address=base_address.format(http_key),
            local_path="{}/{}".format(data_dir, filename),
        )
        f.download()
        print("Done.")
        


# -- define functions: -------------------------------------------------------------------
def downsampled_hematopoiesis(h5ad_path=None, downsample=1, plot=True, quiet=False, data_dir="./content"):
    
    download_tutorial_content(data_dir=data_dir)
    
    if not h5ad_path:
        h5ad_path = PathDict['h5ad']
    
    adata = a.read_h5ad(h5ad_path)[::downsample].copy()
    adata.uns["cmaps"] = {
        "Time point": {2.0: "#3f88c5", 4.0: "#ffba08", 6.0: "#d00000"},
        "Cell type annotation": vinplots.colors.LARRY_in_vitro,
    }
    if not quiet:
        print("\n{}".format(adata))
        if plot:
            overview_plot(adata, groupby=list(adata.uns["cmaps"].keys()))
        
    return adata
