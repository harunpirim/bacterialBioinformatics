from Bio import SeqIO
import json, csv

uni_file = "data\\UniProt.json"
ncbi_file = "data\\sequence.fasta"
output_file = "data\\names.tsv"

ncbi_proteins = {}

list = []
for record in SeqIO.parse(ncbi_file, "fasta"):
    protein_id = record.id.split("_")[2]
    protein_desc = record.description.split("protein=")[1].split("]")[0]
    ncbi_proteins[protein_id] = (protein_desc)                             

with open(uni_file, 'r') as file:
    data = json.load(file)
    results = data["results"]

uni_proteins = {}

for result in results:
    id = result["primaryAccession"]
    embl_id = ""
    evidences = result["organism"]["evidences"]
    for evidence in evidences:
        if evidence["source"] == "EMBL":
            embl_id = evidence["id"]

    uni_proteins[embl_id] = id

matched_proteins = []

for embl, uni in uni_proteins.items():
    if embl in ncbi_proteins.keys():
        matched_proteins.append((uni, embl, ncbi_proteins[embl]))
    else:
        print(uni, embl)

#hardcode missed proteins
matched_proteins.append(("B2KJ98", "ABK41839.1", "thioredoxin"))
matched_proteins.append(("A7YHC3", "ABK41836.1", "gliding protein GldC, partial"))

print(len(matched_proteins))

header = ["UniProt ID", "EMBL ID", "NCBI Protein Description"]

with open(output_file, "w", newline="") as file:
    writer = csv.writer(file, delimiter="\t")
    writer.writerow(header)
    writer.writerows(matched_proteins)