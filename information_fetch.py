import json
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def create_dictionary_from_tsv(tsv_file):
    dictionary = {}

    with open(tsv_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')

        headers = next(reader) 

        for row in reader:
            key = row[0].split("|")[1].strip()  #grabs only protein id
            values = {header: value for header, value in zip(headers[1:], row[1:])}
            dictionary[key] = values

    return dictionary

def get_go_terms(protein_id, results):
    go_terms = []
    for result in results:
        if result["primaryAccession"] == protein_id:
            references = result["uniProtKBCrossReferences"]
            for reference in references:
                if reference["database"] == "GO":
                    go_id = reference["id"]
                    go_term = reference["properties"][0]["value"]
                    go_terms.append((go_id, go_term))
    return go_terms

fasta_file = "uniprot-download_true_format_fasta_query__28_28proteome_3AUP00005856-2023.06.08-20.05.13.50.fasta"
json_file = "uniprot-download_true_format_json_query__28_28proteome_3AUP000058566-2023.06.08-21.46.26.53.json"
motif_file = "data\\motifs_wo_profiles.fasta"
dpc_file = "data\\dpc.tsv"
tpc_file = "data\\tpc.tsv"
paac_file = "data\\paac.tsv" 

# Load the JSON data
with open(json_file) as file:
    data = json.load(file)
    results = data["results"]

#organize motifs
motifs = {}
for record in SeqIO.parse(motif_file, "fasta"):
    id = record.id.split("|")[1].strip()
    seq = str(record.seq)
    if id not in motifs:
        motifs[id] = [seq]
    else:
        motifs[id].append(seq)

        import csv

#load composition features
dpcs = create_dictionary_from_tsv(dpc_file)
tpcs = create_dictionary_from_tsv(tpc_file)
paacs = create_dictionary_from_tsv(paac_file)

proteins = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    protein_id = record.id.split("|")[1].strip()
    seq = record.seq
    a_seq = ProteinAnalysis(seq)

    protein = {
        "seq" : seq,
        "length" : len(record),
        "m_weight": a_seq.molecular_weight(),
        "instab_index": a_seq.instability_index(),
        "isoele_point": a_seq.isoelectric_point(),
        "gravy": a_seq.gravy(),
        "amino_count": a_seq.count_amino_acids(),
        "aromaticity": a_seq.aromaticity(),
        "flexibility": a_seq.flexibility(),
        "sec_sruct_frac": a_seq.secondary_structure_fraction(),
        "go_terms": get_go_terms(protein_id, results),
        "motifs": motifs[id],
        "dpc": dpcs[id],
        "tpc": tpcs[id],
        "paac": paacs[id],
    }

    proteins[protein_id] = protein

for protein_id, protein_info in proteins.items():
    print(f"Protein ID: {protein_id}")
    print(f"Protein Sequence: {protein_info['seq']}")
    print(f"Protein Length: {protein_info['length']}")
    print(f"Molecular Weight: {protein_info['m_weight']}")
    print(f"Instability Index: {protein_info['instab_index']}")
    print(f"Isoelectric Point: {protein_info['isoele_point']}")
    print(f"Gravy: {protein_info['gravy']}")
    print(f"Amino Acid Count: {protein_info['amino_count']}")
    print(f"Aromaticity: {protein_info['aromaticity']}")
    print(f"Flexibility: {protein_info['flexibility']}")
    print(f"Secondary Structure Fraction: {protein_info['sec_sruct_frac']}")
    print(f"Gene Ontology Terms: {protein_info['go_terms']}")
    print(f"Motifs: {protein_info['motifs']}")
    print(f"Dipeptide Composition: {protein_info['dpc']}")
    print(f"Tripeptide Composition: {protein_info['tpc']}")
    print(f"Pseudo Amino Acid Composition: {protein_info['paac']}")
    print()
