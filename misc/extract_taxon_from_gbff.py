from Bio import SeqIO
import pandas as pd

# Input/output files
gbff_file = "plastid.3.genomic.gbff"
output_file = "organelle_taxonomy3.xlsx"

# Define standard taxonomy ranks (trim if you prefer fewer)
standard_ranks = [
    "Superkingdom", "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
]

records_data = []

for record in SeqIO.parse(gbff_file, "genbank"):
    locus = record.name
    organism = record.annotations.get("organism", "NA")
    taxonomy = record.annotations.get("taxonomy", [])
    
    # Organelle
    organelle = "NA"
    for feature in record.features:
        if feature.type == "source":
            organelle = feature.qualifiers.get("organelle", ["NA"])[0]
            break
    
    # Align taxonomy list with ranks
    row = {
        "Locus": locus,
        "Organelle": organelle,
        "Organism": organism,
    }
    
    for i, rank in enumerate(standard_ranks):
        row[rank] = taxonomy[i].strip() if i < len(taxonomy) else ""
    
    records_data.append(row)

# Save to Excel
df = pd.DataFrame(records_data)
df.to_excel(output_file, index=False)

print(f"Processed {len(records_data)} records. Results saved to {output_file}")
