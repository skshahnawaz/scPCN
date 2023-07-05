# scPCN: Revealing pathway connection among cell types from single-cell data

Single-cell analysis has become increasingly popular in recent years, as it allows researchers to explore gene expression and cell types at a high resolution. However, interpreting the vast amount of data generated from the single-cell analysis can be challenging, particularly when trying to identify the functional significance of gene expression changes within specific cell types. One approach to address this challenge is to incorporate pathway information into single-cell analysis. In this research, we present the effective utilization of scPCN (single-cell Pathway Correlation Network) in capturing cell-specific pathways within a specific tissue. Our method  facilitates grouping of cells based on their functionality and also enables the discovery of novel relationships among enriched pathways. scPCN addresses the critical requirement of elucidating the interplay between pathways across diverse cell types within a specific tissue under various conditions. These relationships establish a comprehensive pathway interaction model, which contributes to a biologically driven understanding of phenotypes, reduces noise, and enhances overall performance. Therefore, scPCN represents a powerful advancement in interpreting results obtained through gene set methods.

![flowchart (1)-1](https://github.com/skshahnawaz/scpcn/assets/52563824/52dd4a1f-e1d3-49ed-9cfb-cd52bf2ba191)

**Fig:** The schematic diagram of the proposed method

## Steps to run
Run `python main.py` after installing dependencies via `pip install requirements.txt`

## Clustering Output

`sce.tl.phenograph(self.pathway_activity_matrix2, clustering_algo="louvain", k=30)`

![umap_pbmc_3k-1](https://github.com/skshahnawaz/scpcn/assets/52563824/4493b1c2-49bb-44a0-8621-2466901dbc04)

## Biological Significance

<p float="left" align="middle">
  <img src="https://github.com/skshahnawaz/scpcn/assets/52563824/72aeb283-c63d-4387-9b8e-10b9566311ab" width="33%" />
  <img src="https://github.com/skshahnawaz/scpcn/assets/52563824/78636cd3-5e36-4de6-bbaa-5efc63022d46" width="33%" /> 
  <img src="https://github.com/skshahnawaz/scpcn/assets/52563824/9619bb67-0680-43c6-9ae4-7a0de2eb34ec" width="33%" />
</p>

