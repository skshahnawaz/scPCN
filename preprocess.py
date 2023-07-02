import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

class Preprocess:

    def __init__(self):
        # self.dataset_url = dataset_url
        print()

    def configure_scanpy(self):
        sc.settings.verbosity = 3
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor='white')

    def create_anndata(self):
        self.adata = sc.read_10x_mtx('dataset/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)

    def preprocess_data(self):
        self.adata.var_names_make_unique()
        sc.pl.highest_expr_genes(self.adata, n_top=20, save='.pdf')

        # Basic filtering
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)

        # annotate the group of mitochondrial genes as 'mt'
        self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-')

        # Calculate QC Metrics
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # Plot of computed quality measures
        sc.pl.violin(self.adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='.pdf')

        # Removing cells that have too many mitochondrial genes expressed or too many total counts
        # sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt')
        # sc.pl.scatter(self.adata, x='total_counts', y='n_genes_by_counts')

        # Actual filtering by slicing the AnnData object
        self.adata = self.adata[self.adata.obs.n_genes_by_counts < 2500, :]
        self.adata = self.adata[self.adata.obs.pct_counts_mt < 5, :]

        # Count normalize the data metrics
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        
        # Logarithmize the data
        sc.pp.log1p(self.adata)

        # Export the anndata to CSV
        self.adata.to_df().to_csv('./anndata.csv')

        # Extract the highly variable genes
        sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(self.adata, save='.pdf')

        self.adata.raw = self.adata

        # Regress out the data
        self.adata = self.adata[:, self.adata.var.highly_variable]
        sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(self.adata, max_value=10)

        return self.adata

        # print(self.adata)