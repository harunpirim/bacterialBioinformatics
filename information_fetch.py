import json
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# -- Variables: Protein Dictionary ---
proteins = {}

#--- Variables: File Locations ---
protein_file = "data\\proteins-2023.06.08-20.05.13.50.fasta"
uniprot_file = "data\\UniProt.json"
motif_file = "data\\Prosite_wo.json"
ifeature_file = "data\\iFeature.json"

#--- Functions: Assign Annotations (UniProt) ---
def assign_annotations():
    # Load JSON data from UniProt file
    with open(uniprot_file, 'r') as file:
        data = json.load(file)
        results = data["results"]
    
    global proteins
    for result in results:
        id = result["primaryAccession"]
        if "uniProtKBCrossReferences" in result: proteins[id]["go_terms"] = get_go_terms(result)
        if "comments" in result: proteins[id]["subcell_locations"] = get_subcellular_locations(result)
        if "features" in result: proteins[id]["transmembrane"] = get_transmembrane(result)
        if "keywords" in result: proteins[id]["binding_preference"] = get_binding_preference(result)

def get_go_terms(protein_info):
    go_terms = []
    references = protein_info["uniProtKBCrossReferences"]
    for reference in references:
        if reference["database"] == "GO":
            go_id = reference["id"]
            go_term = reference["properties"][0]["value"]
            go_terms.append((go_id, go_term))
    return go_terms

def get_subcellular_locations(protein_info):
    subcell_locations = []
    comments = protein_info["comments"]
    for comment in comments:
        if comment["commentType"] == "SUBCELLULAR LOCATION":
            locations = comment["subcellularLocations"]
            for location in locations:
                subcell_locations.append(location["location"]["value"])
    return subcell_locations

def get_transmembrane(protein_info):
    transmembrane = 0
    features = protein_info["features"]
    for feature in features:
        if feature["type"] == "Transmembrane":
            transmembrane = 1
    return transmembrane

def get_binding_preference(protein_info):
    bp = []
    keywords = protein_info["keywords"]
    for keyword in keywords:
        kw = keyword["name"]
        if (kw == "DNA-binding" or kw == "RNA-binding"):
            bp.append("DNA/RNA-binding")
        elif (kw == "Nucleotide-binding" or kw == "Metal-binding"):
            bp.append(kw)
    return bp

#--- Functions: Assign Motifs (Prosite) ---
def assign_motifs():
    global proteins
    with open(motif_file, 'r') as file:
        data = json.load(file)
    for id in data:
        proteins[id]['motifs'] = data[id]

#--- Functions: Assign iFeature Data (iFeature) ---
def assign_iFeature_data():
    global proteins
    with open(ifeature_file, 'r') as file:
        features = json.load(file)
    for feature, ids in features.items():
        for id, value in ids.items():
            proteins[id][feature] = value
                

#--- Main: Generate Protein Dictionary ---
for record in SeqIO.parse(protein_file, "fasta"):
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
        "go_terms": [],
        "motifs": [],
        "dpc": [],
        "tpc": [],
        "paac": [],
        "cdtc": [],
        "cdtd": [],
        "cdtt": [],
        "ctriad": [],
        "gaac": [],
        "moran": [],
        "subcell_locations": [],
        "transmembrane": 0,
        "binding_preference": [],
    }
    proteins[protein_id] = protein

assign_annotations()
assign_motifs()
assign_iFeature_data()


#--- Test: Print Proteins ---
for protein_id, protein_info in proteins.items():
    print(f"Protein ID: {protein_id}")
    # print(f"Protein Sequence: {protein_info['seq']}")
    # print(f"Protein Length: {protein_info['length']}")
    # print(f"Molecular Weight: {protein_info['m_weight']}")
    # print(f"Instability Index: {protein_info['instab_index']}")
    # print(f"Isoelectric Point: {protein_info['isoele_point']}")
    # print(f"Gravy: {protein_info['gravy']}")
    # print(f"Amino Acid Count: {protein_info['amino_count']}")
    # print(f"Aromaticity: {protein_info['aromaticity']}")
    # print(f"Flexibility: {protein_info['flexibility']}")
    # print(f"Secondary Structure Fraction: {protein_info['sec_sruct_frac']}")
    # print(f"Gene Ontology Terms: {protein_info['go_terms']}")
    # print(f"Motifs: {protein_info['motifs']}")
    # print(f"Dipeptide Composition: {protein_info['dpc']}")
    # print(f"Tripeptide Composition: {protein_info['tpc']}")
    # print(f"Pseudo Amino Acid Composition: {protein_info['paac']}")
    # print(f"Composition: {protein_info['cdtc']}")
    # print(f"Distribution: {protein_info['cdtd']}")
    # print(f"Translation: {protein_info['cdtt']}")
    # print(f"Conjoint Triad: {protein_info['ctriad']}")
    # print(f"Grouped Amino Acid Composition: {protein_info['gaac']}")
    # print(f"Moran Autocorrelation: {protein_info['moran']}")
    # print(f"Subcellular Locations: {protein_info['subcell_locations']}")
    # print(f"Transmembrane?: {protein_info['transmembrane']}")
    # print(f"Binding Preference: {protein_info['binding_preference']}")
    print()
