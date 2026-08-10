import pandas as pd
import os

dna_sequence = "acgtacgtggccaattxyz"

try:
    dna_sequence = dna_sequence.upper()

    counts = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0,
    }

    invalid_chararacters = []

    for nucleotide in dna_sequence:
        if nucleotide in counts:
            counts[nucleotide] += 1
        else:
            invalid_chararacters.append(nucleotide)
    
    metrics = pd.DataFrame({
        "Nucleotide": ["A", "C", "G", "T"],
        "Count": [
            counts["A"],
            counts["C"],
            counts["G"],
            counts["T"],
        ]
    })

    os.makedirs("results", exist_ok=True)

    metrics.to_csv("results/dna_metrics.csv", index=False)

    print("Dna sequence processed successfully")
    print(metrics)

    if invalid_chararacters:
        print(f"Warning: Invalid characters: {invalid_chararacters}")

except Exception as error:
    print(f"An error occurred: {error}")
   