from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# count = 0

proteins = {}
for record in SeqIO.parse("GCF_001534645.1_ASM153464v1_protein.faa", "fasta"):
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
    }
    # count+=1
    # print(count)

    proteins[record.id] = protein

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
    print()