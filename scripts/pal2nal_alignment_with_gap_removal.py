import os
import subprocess

# === User Configuration ===
protein_dir = "aligned_PROTEIN_10093"
cds_dir = "10093_CDS"
pal2nal_path = "PAL2NAL/pal2nal.pl"
output_dir = "codon_aligned_output_with_gap_removal_option"

# Create output folder
os.makedirs(output_dir, exist_ok=True)

# Log file
log_file = os.path.join(output_dir, "pal2nal_run_log.txt")
log = open(log_file, "w")

# List protein alignment files
protein_files = [f for f in os.listdir(protein_dir) if f.endswith(".fa")]

for prot_file in protein_files:
    og_id = prot_file.split("_")[0]  # e.g., OG0000151
    cds_file = f"{og_id}.fa"
    prot_path = os.path.join(protein_dir, prot_file)
    cds_path = os.path.join(cds_dir, cds_file)

    if not os.path.exists(cds_path):
        log.write(f"[MISSING] CDS not found for {og_id}\n")
        continue

    output_file = os.path.join(output_dir, f"{og_id}_codon_aln.fasta")

    cmd = [
        "perl", pal2nal_path,
        prot_path, cds_path,
        "-output", "fasta",
        "-nogap",
        "-codontable", "1"  # Universal code for fungal nuclear genes
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(output_file, "w") as out:
            out.write(result.stdout)
        log.write(f"[OK] {og_id}\n")
    except subprocess.CalledProcessError as e:
        log.write(f"[ERROR] {og_id} - {e.stderr}\n")

log.close()
print(f"\n✅ Codon alignments completed! Check folder: {output_dir}")
print(f"📝 Log file saved at: {log_file}")

