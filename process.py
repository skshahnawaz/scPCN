import scanpy as sc
import json
import requests
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
# import phenograph
import scanpy.external as sce
import time

# adatacheck -> cell_gene_df

class Process:
    
    def __init__(self, cell_gene_expression_data):
    # def __init__(self):
        self.adata = cell_gene_expression_data
        print(self.adata)
        # print()

    def convert_to_dataframe(self):
        print("-----Converting to dataframe-----")
        self.cell_gene_df = self.adata.to_df()
        # Create a list of the genes
        self.genesyms = list(self.cell_gene_df.columns)
        print(f'{len(self.genesyms)} genes found')

    def fetch_enrichr_data(self):
        print("-----Fetching EnrichR data-----")
        ENRICHR_URL_POST = 'https://maayanlab.cloud/Enrichr/addList'
        genes_str = '\n'.join(self.genesyms)
        description = 'My Gene List'
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
            }
        enrichr_list_response = requests.post(ENRICHR_URL_POST, files=payload)
        if not enrichr_list_response.ok:
            raise Exception('Error analyzing gene list')
        
        enrichr_list_data = json.loads(enrichr_list_response.text)

        ENRICHR_URL_GET = 'https://maayanlab.cloud/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'

        user_list_id = enrichr_list_data['userListId']

        gene_set_library = 'KEGG_2016'

        enrichr_fetch_response = requests.get(ENRICHR_URL_GET + query_string % (user_list_id, gene_set_library))
    
        if not enrichr_fetch_response.ok:
            raise Exception('Error fetching enrichment results')
        
        self.enrichr_fetch_data = json.loads(enrichr_fetch_response.text)
        print("-----Fetched EnrichR data-----")

    def process_enrichr_data(self):
        print("-----Processing EnrichR data-----")
        self.pathway_geneset = pd.DataFrame(columns=['Pathway_Name','adjusted_p_value','p_value','z_score','genes'])
        for result in self.enrichr_fetch_data['KEGG_2016']:
            pathway_name = result[1].split("Homo sapiens")[0]
            adjusted_p_value = result[6]
            p_value = result[2]
            z_score = result[3]
            genelist = list(result[5])
            if (p_value < 0.05):
                exclusions = ["disease","infection","cancer"]
                if (not any(ele in pathway_name for ele in exclusions)):
                    self.pathway_geneset.loc[len(self.pathway_geneset)] = [pathway_name,adjusted_p_value,p_value,z_score,genelist]
        self.pathway_geneset.to_csv('pathway_geneset.csv')
        print(f'Shape of pathway_geneset data : {self.pathway_geneset.shape}')

    def find_activated_pathways(self):
        # for storing cell-gene expression values of the genes associated with the pathway
        pathway_cell_gene_matrix = pd.DataFrame()
        self.pathway_activity_matrix = pd.DataFrame()

        self.pathway_activated_genes = pd.DataFrame()

        for row in self.pathway_geneset.iterrows():
            # self.start = time.time()
            pathway_name = row[1]['Pathway_Name'] 
            genes = row[1]['genes']
            pathway_cell_gene_matrix = self.cell_gene_df[genes]

            # calculate the gene-gene correlation matrix
            gene_correlation_matrix = pathway_cell_gene_matrix.corr(method='pearson')

            # Perform the thresholding
            for gene in list(gene_correlation_matrix.columns):
                gene_correlation_matrix[gene] = np.where(gene_correlation_matrix[gene] > 0.2, 1, gene_correlation_matrix[gene])
                gene_correlation_matrix[gene] = np.where(gene_correlation_matrix[gene] <= 0.2, 0, gene_correlation_matrix[gene])

            gene_correlation_matrix.to_csv('./gene_correlation_matrix.csv') 

            # Find activated genes for a pathway
            activated_genes = {''}
            for row, col in gene_correlation_matrix.iterrows():
                for g in list(gene_correlation_matrix.columns):
                    if (col[g] == 1.0):
                        if (row != g):
                            activated_genes.add(row)
                            activated_genes.add(g)

            activated_genes.remove('')

            # Add to the pathway-activated_genes df
            self.pathway_activated_genes = self.pathway_activated_genes.append({'Pathway Name':pathway_name, 'Activated Genes':activated_genes},ignore_index=True)

            if (len(activated_genes) == 0):
                continue

            cell_gene_activated = pd.DataFrame()
            for row in self.cell_gene_df.iterrows():
                cell_gene_activated = self.cell_gene_df[list(activated_genes)]

            cell_gene_activated.to_csv('./cell_gene_activated.csv')

            for row,col in cell_gene_activated.iterrows():
                exprsum = 0
                for gene in list(cell_gene_activated.columns):
                    exprsum += col[gene]
                exprsum /= math.sqrt(len(activated_genes))
                self.pathway_activity_matrix.loc[row,pathway_name] = exprsum

            # self.end = time.time()

        self.pathway_activated_genes.to_csv('./pathway_activated_genes.csv')

        # self.end = time.time()
        # print(f'Time is : {self.end - self.start}') # remove it

        # print(f'Shape is : {self.pathway_activity_matrix.shape}') #remove it

        self.pathway_activity_matrix.to_csv('./pathway_activity_matrix.csv')

    def cluster_data(self):
        # implementing the PhenoGraph clustering

        self.pathway_activity_matrix2 = pd.read_csv('./pathway_activity_matrix.csv', index_col=0)

        print("Clustering...")
        sce.tl.phenograph(self.pathway_activity_matrix2, clustering_algo="louvain", k=30)

        # save the clustering results as csv (with ranking)
        # end