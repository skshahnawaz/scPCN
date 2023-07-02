import os

os.environ['KMP_DUPLICATE_LIB_OK']='True'

import scvi
scvi.data.brainlarge_dataset(save_path='data/')