import os
import gdown

# Run this file to download data necessary for information_fetch

os.makedirs('data', exist_ok=True)

ifeature_url = "https://drive.google.com/uc?id=1p0jY2Q1HRO7GE9vcbGngzCE24hGyIs0d"
prosite_wo_url = "https://drive.google.com/uc?id=1C3KqmV8FewuxYQ31SqLQj_xzyfliVqYg"
proteins_url = "https://drive.google.com/uc?id=1dTuW5r8oKeNvHwukVqD7fMmf5lYVmx5M"
uniprot_url = "https://drive.google.com/uc?id=1vLWv6ft0vJ9_jUJe_LJINptAVBJZUwOe"

ifeature_output = "data\\iFeature.json"
prosite_wo_output = "data\\Prosite_wo.json"
proteins_output = "data\\proteins-2023.06.08-20.05.13.50.fasta"
uniprot_output = "data\\UniProt.json"

gdown.download(ifeature_url, ifeature_output, quiet=False)
gdown.download(prosite_wo_url, prosite_wo_output, quiet=False)
gdown.download(proteins_url, proteins_output, quiet=False)
gdown.download(uniprot_url, uniprot_output, quiet=False)