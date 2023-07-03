from download import Download
from preprocess import Preprocess
from process import Process
import time

# https://discourse.scverse.org/t/how-to-fix-this-error-when-install-scanpy/812/2

dataset_dict = {'url': 'http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'}

# Download the dataset
download = Download(dataset_dict)
download.download()

# preprocess = Preprocess()

# preprocess.configure_scanpy()
# preprocess.create_anndata()
# cell_gene_data = preprocess.preprocess_data()

# # print(cell_gene_data)

# # process = Process()
# process = Process(cell_gene_data) #only if data to be loaded from saved file
# process.convert_to_dataframe()
# process.fetch_enrichr_data()
# process.process_enrichr_data()
# process.find_activated_pathways()
# process.cluster_data()