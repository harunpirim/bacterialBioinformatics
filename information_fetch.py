from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import requests

def get_gene_ontology_annotations(protein_id):
    url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId={protein_id}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        annotations = data["results"]
        return annotations
    else:
        print(f"Failed to retrieve gene ontology annotations for protein {protein_id}.")
        return []

fasta_file = "uniprot-download_true_format_fasta_query__28_28proteome_3AUP00005856-2023.06.08-20.05.13.50.fasta"

# count = 0
proteins = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    id = record.id.split("|")[1].strip()
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
        "go_annotations": get_gene_ontology_annotations(id)
    }

    proteins[id] = protein

    # count+=1
    # print(count)

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
    # can also access other GO information such as goEvidence
    for annotation in protein_info['go_annotations']:
        go_term = annotation["goId"]
        go_aspect = annotation["goAspect"]
        print(f"Gene Ontology term: {go_term} ({go_aspect})")
    print()