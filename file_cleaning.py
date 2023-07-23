import csv
import json
from Bio import SeqIO


def nested_dict_to_json(nested_dict, file_path):
    with open(file_path, 'w') as json_file:
        json.dump(nested_dict, json_file, indent=4)

def json_to_nested_dict(file_path):
    with open(file_path, 'r') as json_file:
        nested_dict = json.load(json_file)
    return nested_dict

#--- iFeature cleaning ---
def create_dictionary_from_iFeature_tsv(ifeature_file):
    dictionary = {}
    with open(ifeature_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader) 
        for row in reader:
            key = row[0].split("|")[1].strip()  #grabs only protein id
            values = {header: value for header, value in zip(headers[1:], row[1:])}
            dictionary[key] = values
    return dictionary

ifeature_files = {
"dpc" : "raw_data\\dpc.tsv",
#"tpc" : "data\\tpc.tsv",
"paac" : "raw_data\\paac.tsv",
"ctdc" : "raw_data\\ctdc.tsv",
"ctdt" : "raw_data\\ctdt.tsv",
"ctdd" : "raw_data\\ctdd.tsv",
"ctriad" : "raw_data\\ctriad.tsv",
"gaac" : "raw_data\\gaac.tsv",
"moran" : "raw_data\\moran.tsv"
}

combined_dictionary = {}

for key, value in ifeature_files.items():
    combined_dictionary[key] = create_dictionary_from_iFeature_tsv(value)

nested_dict_to_json(combined_dictionary, "data\\train\\train_iFeature.json")

# #--- Prosite cleaning ---
# motif_file = "data\\motifs_wo_profiles.fasta"

# motifs = {}
# for record in SeqIO.parse(motif_file, "fasta"):
#     id = record.id.split("|")[1].strip()
#     seq = str(record.seq)
#     if id not in motifs:
#         motifs[id] = [seq]
#     else:
#         motifs[id].append(seq)

# # nested_dict_to_json(motifs, "Prosite.json")

# #--- UniProt cleaning ---
# uniprot_file = "raw_data\\uniprot-download_true_format_json_query__28_28proteome_3AUP000058566-2023.06.08-21.46.26.53.json"

# dict = json_to_nested_dict(uniprot_file)
# # nested_dict_to_json(dict, "data\\UniProt.json")



