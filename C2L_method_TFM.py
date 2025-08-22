import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import squidpy as sq
import cell2location as c2l

import torch

adata_vis = sc.read("visium_d30_porta5_2.rds.h5ad")
# Set sample name
adata_vis.obs['sample'] = 'Visium_d30_porta5'

# Set gene names
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('SYMBOL', drop=False, inplace=True)

adata_vis.var['SYMBOL'] = adata_vis.var.index

# Single cell reference data
adata_ref = sc.read("seurat.anotado_C2L.h5ad")
# Set gene names
adata_ref.var['SYMBOL'] = adata_ref.var.index

adata_ref.X = adata_ref.layers["counts_SCT"]
sc.pl.tsne(adata_ref, color="orig.ident")

# Filter genes from the reference object
selected = c2l.utils.filtering.filter_genes(
    adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
)

# Prepare reference for the regression model
c2l.models.RegressionModel.setup_anndata(adata=adata_ref,
                        batch_key='orig.ident',
                        labels_key='Final_Labels'
                       )

# Regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup
mod.view_anndata_setup()

mod.train(max_epochs=250)
mod.plot_history()

ref_run_name = "my_model_run_250_epoch_bueno"
# Estimated cell abundance
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 100, 'batch_size': 250}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = c2l.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={'num_samples': 100, 'batch_size': 250}
)

mod.plot_QC()

ref_run_name = "my_model_run_250_epoch_bueno"
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = c2l.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# Export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# Find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# Prepare query for cell2location model
c2l.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# Create and train the model
mod = c2l.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # Expected average cell abundance: tissue-dependent
    N_cells_per_location=14,
    detection_alpha=200
)
mod.view_anndata_setup()

mod.train(max_epochs=10000,
          # train using full data
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
         )

# plot ELBO loss history during training
mod.plot_history()
plt.legend(labels=['full data training']);

run_name = "my_model_run_2_10000_epoch"
mod = c2l.models.Cell2location.load(f"{run_name}", adata_vis)

adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 10, 'batch_size': mod.adata.n_obs}
)

run_name = "my_model_run_2_10000_epoch_14_cellabund"
# Export the estimated cell abundance
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 10, 'batch_size': mod.adata.n_obs}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

mod.plot_QC()

from PIL import Image
import numpy as np

# Tissue image
image_path = "detected_tissue_image.jpg" 
tissue_image = Image.open(image_path)
image_flipud = np.flipud(tissue_image)
image_rotated = np.rot90(image_flipud, k=3)
# Image to array
tissue_image_array = np.array(image_rotated)

# Insert the image and the scale factor into query object
adata_vis.uns['spatial'] = {
    'Visium_d30_p5': { 
        'images': {
            'hires': tissue_image_array 
        },
        'scalefactors': {
            'spot_diameter_fullres': 15.02406,  
            'tissue_hires_scalef': 1 
        },
        'metadata': {
            'chemistry_description': 'Custom',
            'software_version': 'Manual'
        }
    }
}

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
adata_vis.obs['library_id'] = 'Visium_d30_p5'
slide = adata_vis[adata_vis.obs['library_id'] == 'Visium_d30_p5'].copy()
slide.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obs.loc[slide.obs.index, adata_vis.uns['mod']['factor_names']]

import matplotlib as mpl
# Plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
#magma, viridis, plasma, inferno, cividis
    sc.pl.spatial(slide, cmap="magma",
                  color=adata_vis.uns['mod']['factor_names'],
                  ncols=4, size=0.8,
                  img_key='hires',
                  vmin=0, vmax=30,
                  save="spatial_plot.png"
                 )
