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

### Unveiling Pathway Crosstalk through scPCN

Cervical Caner             |  Endometrial Cancer             |  Ovarian Cancer
:-------------------------:|:-------------------------:|:-------------------------:
![cerivical-1](https://github.com/skshahnawaz/scpcn/assets/52563824/2c390653-b542-41b2-b3d3-0a42c47f3503)  |  ![endo-1](https://github.com/skshahnawaz/scpcn/assets/52563824/16ddf310-2e18-4bd3-853f-c901e7ecb8dc)  |  ![ovary-1](https://github.com/skshahnawaz/scpcn/assets/52563824/cf20aea3-f075-4bec-ad62-d3faab09ee22)



