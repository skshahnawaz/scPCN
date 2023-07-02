import os
import requests
import tarfile

class Download:

    def __init__(self, dataset):
        self.url = dataset['url']
        print(f'URL received is : {self.url}')

    def download(self):
        path = "dataset"
        # Create a directory if it does not exist
        if not os.path.exists(path):
            os.makedirs(path)
            print("New directory is created")

        response = requests.get(self.url, stream=True)
        file = tarfile.open(fileobj=response.raw, mode="r|gz")
        file.extractall(path="./dataset")
