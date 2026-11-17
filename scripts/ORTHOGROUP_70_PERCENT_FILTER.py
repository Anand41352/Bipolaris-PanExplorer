#!/usr/bin/env python3
import pandas as pd
import math

# -----------------------------
# Input file
# -----------------------------
input_file = "Orthogroups.GeneCount.tsv"
output_file = "Orthogroups_filtered_70pct_singlecopy.tsv"

# -----------------------------
# Read the file
# -----------------------------
df = pd.read_csv(input_file, sep="\t")

# Drop 'Total' column if it exists
if 'Total' in df.columns:
    df = df.drop(columns=['Total'])

# Extract strain columns (all except Orthogroup)
strain_cols = [c for c in df.columns if c != 'Orthogroup']
n_species = len(strain_cols)
threshold = math.ceil(0.7 * n_species)

print(f"Total species (columns): {n_species}")
print(f"Single-copy threshold (70%): {threshold} species")

# -----------------------------
# Filtering Step 1: Present in all species (no zeros)
# -----------------------------
present_all = df[(df[strain_cols] >= 1).all(axis=1)]

# -----------------------------
# Filtering Step 2: Single-copy (count == 1) in ≥70% of species
# -----------------------------
def count_single_copies(row):
    return (row[strain_cols] == 1).sum()

present_all['single_copy_count'] = present_all.apply(count_single_copies, axis=1)
filtered = present_all[present_all['single_copy_count'] >= threshold]

# -----------------------------
# Save the filtered data
# -----------------------------
filtered.to_csv(output_file, sep="\t", index=False)

# -----------------------------
# Summary
# -----------------------------
print("✅ Filtering complete!")
print(f"Input orthogroups: {len(df)}")
print(f"Present in all species: {len(present_all)}")
print(f"After 70% single-copy filtering: {len(filtered)}")
print(f"Filtered file saved as: {output_file}")

