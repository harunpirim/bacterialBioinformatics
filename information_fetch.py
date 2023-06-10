import json
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def get_go_ids(protein_id, results):
    go_ids = []
    for result in results:
        if result["primaryAccession"] == protein_id:
            references = result["uniProtKBCrossReferences"]
            for reference in references:
                if reference["database"] == "GO":
                    go_id = reference["id"]
                    go_ids.append(go_id)
    return go_ids

fasta_file = "uniprot-download_true_format_fasta_query__28_28proteome_3AUP00005856-2023.06.08-20.05.13.50.fasta"
json_file = "uniprot-download_true_format_json_query__28_28proteome_3AUP000058566-2023.06.08-21.46.26.53.json"

# Load the JSON data
with open(json_file) as file:
    data = json.load(file)
    results = data["results"]

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
        "go_ids": get_go_ids(protein_id, results)
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
    for go_id in protein_info['go_ids']:
        print(f"Gene Ontology IDs: {go_id}")
    print()
