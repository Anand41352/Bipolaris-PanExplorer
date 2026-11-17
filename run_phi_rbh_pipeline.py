#!/usr/bin/env python3

"""
PHI-base homology and reciprocal best-hit (RBH) pipeline.

Input:
    core.fasta
    softcore.fasta
    shell.fasta
    cloud.fasta
    phi-base_current.fas

Primary filters:
    Identity >= 40%
    E-value <= 1e-5
    Query coverage >= 40%
    Subject coverage >= 40%

Sensitivity analysis:
    Query coverage >= 50%
    Subject coverage >= 50%

Outputs:
    - Combined B. sorokiniana proteome
    - Strain-specific proteomes
    - DIAMOND databases
    - Forward alignments with coverage
    - Coverage-filtered forward matches
    - Best qualifying forward hit per query
    - Reverse best-hit results
    - RBH classifications at 40%/40% and 50%/50%
    - Summary statistics

Coverage is calculated from the aligned coordinate span divided
by the full protein length.

Important:
    RBH is supporting evidence for candidate orthology, not
    definitive proof of orthology or conserved phenotype.
"""

import csv
import shutil
import subprocess
import sys

from collections import defaultdict
from pathlib import Path


# ============================================================
# 1. CONFIGURATION
# ============================================================

BASE_DIR = Path(
    "/home/naren/Desktop/ANAND/PAN_GENOME/"
    "MANUSCRIPT/REVISION/Anand papers/Pangenome/phi"
)

OUTPUT_DIR = BASE_DIR / "phi_rbh_analysis"

QUERY_FASTAS = [
    BASE_DIR / "core.fasta",
    BASE_DIR / "softcore.fasta",
    BASE_DIR / "shell.fasta",
    BASE_DIR / "cloud.fasta",
]

PHI_FASTA = BASE_DIR / "phi-base_current.fas"

IDENTITY_MIN = 40.0
EVALUE_MAX = 1e-5

COVERAGE_THRESHOLDS = [
    (40.0, 40.0),
    (50.0, 50.0),
]

# Search up to this many PHI-base targets per query.
# Increase this if you want to inspect more candidate matches.
MAX_TARGETS = 50

THREADS = 8

# DIAMOND alignment output fields.
OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


# ============================================================
# 2. GENERAL UTILITIES
# ============================================================

def run_command(command):
    """Execute a command and stop if it fails."""

    print("\nRUNNING:")
    print(" ".join(str(x) for x in command))

    subprocess.run(
        [str(x) for x in command],
        check=True
    )


def read_fasta(path):
    """
    Read a FASTA file.

    Returns:
        Dictionary: sequence ID -> amino-acid sequence.

    The sequence ID is the first whitespace-delimited token
    after the FASTA '>' symbol.
    """

    sequences = {}
    current_id = None
    chunks = []

    with open(path, "r") as handle:

        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:
                    sequences[current_id] = "".join(chunks)

                current_id = line[1:].split()[0]
                chunks = []

                if current_id in sequences:
                    raise ValueError(
                        f"Duplicate FASTA identifier: {current_id}"
                    )

            else:

                if current_id is None:
                    raise ValueError(
                        f"Sequence encountered before FASTA header "
                        f"in {path}"
                    )

                chunks.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(chunks)

    return sequences


def write_fasta(records, path):
    """Write a dictionary of sequence IDs and sequences to FASTA."""

    with open(path, "w") as handle:

        for sequence_id, sequence in records.items():

            handle.write(f">{sequence_id}\n")

            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(rows, path, fieldnames):
    """Write dictionaries to a TSV file."""

    with open(path, "w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    """Read a TSV file into a list of dictionaries."""

    with open(path, "r", newline="") as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        return list(reader)


def strain_from_protein_id(protein_id):
    """
    Extract strain from IDs such as:

        OG0014125|D2_01272
        OG0014426|WAI3285_06940
        OG0000002|BRIP10943a_02345

    Assumes the strain identifier is the part of the gene ID
    before its first underscore.
    """

    if "|" not in protein_id:
        raise ValueError(
            f"Unexpected protein ID format: {protein_id}"
        )

    gene_id = protein_id.split("|", 1)[1]

    if "_" not in gene_id:
        raise ValueError(
            f"Cannot identify strain in protein ID: {protein_id}"
        )

    return gene_id.split("_", 1)[0]


def alignment_coverage(start, end, full_length):
    """
    Calculate coordinate-span coverage as a percentage.

    Coordinates in DIAMOND output are 1-based and inclusive.
    """

    if full_length <= 0:
        return 0.0

    aligned_span = abs(int(end) - int(start)) + 1

    return 100.0 * aligned_span / full_length


def passes_filters(row, query_cov_min, subject_cov_min):
    """Check all four primary sequence filters."""

    return (
        float(row["pident"]) >= IDENTITY_MIN
        and float(row["evalue"]) <= EVALUE_MAX
        and float(row["query_coverage"]) >= query_cov_min
        and float(row["subject_coverage"]) >= subject_cov_min
    )


# ============================================================
# 3. VALIDATE INPUTS AND SOFTWARE
# ============================================================

print("\n" + "=" * 65)
print("PHI-BASE HOMOLOGY AND RBH PIPELINE")
print("=" * 65)

if shutil.which("diamond") is None:
    sys.exit(
        "ERROR: DIAMOND was not found in PATH. "
        "Activate the environment containing DIAMOND."
    )

for path in QUERY_FASTAS + [PHI_FASTA]:

    if not path.exists():
        sys.exit(f"ERROR: Required input file not found: {path}")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

PROTEOME_DIR = OUTPUT_DIR / "strain_proteomes"
DB_DIR = OUTPUT_DIR / "diamond_databases"
FORWARD_DIR = OUTPUT_DIR / "forward_search"
REVERSE_DIR = OUTPUT_DIR / "reverse_search"

for directory in [
    PROTEOME_DIR,
    DB_DIR,
    FORWARD_DIR,
    REVERSE_DIR,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True
    )


# ============================================================
# 4. LOAD AND COMBINE B. SOROKINIANA PROTEINS
# ============================================================

print("\nLoading B. sorokiniana proteins...")

all_bs_proteins = {}
strain_proteins = defaultdict(dict)

for fasta_path in QUERY_FASTAS:

    print(f"Reading: {fasta_path.name}")

    sequences = read_fasta(fasta_path)

    for protein_id, sequence in sequences.items():

        if protein_id in all_bs_proteins:
            raise ValueError(
                f"Protein occurs more than once across FASTA files: "
                f"{protein_id}"
            )

        strain = strain_from_protein_id(protein_id)

        all_bs_proteins[protein_id] = sequence
        strain_proteins[strain][protein_id] = sequence


if not all_bs_proteins:
    sys.exit("ERROR: No B. sorokiniana proteins were found.")


print(f"\nTotal B. sorokiniana proteins: {len(all_bs_proteins)}")
print(f"Number of strains identified: {len(strain_proteins)}")

for strain in sorted(strain_proteins):

    print(
        f"  {strain}: "
        f"{len(strain_proteins[strain])} proteins"
    )


# Write combined proteome.

ALL_BS_FASTA = OUTPUT_DIR / "all_Bsorokiniana_proteins.fasta"

write_fasta(
    all_bs_proteins,
    ALL_BS_FASTA
)


# Write strain-specific proteomes.

for strain, proteins in strain_proteins.items():

    strain_fasta = PROTEOME_DIR / f"{strain}.fasta"

    write_fasta(
        proteins,
        strain_fasta
    )


# ============================================================
# 5. LOAD PHI-BASE PROTEINS
# ============================================================

print("\nLoading PHI-base reference proteins...")

phi_proteins = read_fasta(PHI_FASTA)

print(f"Total PHI-base proteins: {len(phi_proteins)}")

if not phi_proteins:
    sys.exit("ERROR: No PHI-base proteins were found.")


# Save a simple mapping of sequence lengths.

bs_lengths = {
    protein_id: len(sequence)
    for protein_id, sequence in all_bs_proteins.items()
}

phi_lengths = {
    protein_id: len(sequence)
    for protein_id, sequence in phi_proteins.items()
}


# ============================================================
# 6. BUILD DIAMOND DATABASES
# ============================================================

print("\nBuilding DIAMOND databases...")

PHI_DB = DB_DIR / "phi_base_db"

run_command([
    "diamond",
    "makedb",
    "--in", PHI_FASTA,
    "--db", PHI_DB,
])

STRAIN_DBS = {}

for strain in sorted(strain_proteins):

    strain_fasta = PROTEOME_DIR / f"{strain}.fasta"

    strain_db = DB_DIR / f"{strain}_db"

    run_command([
        "diamond",
        "makedb",
        "--in", strain_fasta,
        "--db", strain_db,
    ])

    STRAIN_DBS[strain] = strain_db


# ============================================================
# 7. FORWARD DIAMOND SEARCH
# ============================================================

print("\n" + "=" * 65)
print("FORWARD SEARCH: B. sorokiniana -> PHI-base")
print("=" * 65)

FORWARD_RAW = FORWARD_DIR / "forward_all_alignments.tsv"

run_command([
    "diamond",
    "blastp",
    "--query", ALL_BS_FASTA,
    "--db", PHI_DB,
    "--out", FORWARD_RAW,
    "--outfmt", "6",
    *OUTFMT_FIELDS,
    "--evalue", EVALUE_MAX,
    "--max-target-seqs", MAX_TARGETS,
    "--sensitive",
    "--threads", THREADS,
])


# ============================================================
# 8. CALCULATE FORWARD COVERAGE
# ============================================================
# ============================================================
# 8. CALCULATE FORWARD COVERAGE
# ============================================================

print("\nCalculating forward query and subject coverage...")

forward_rows = []

with open(FORWARD_RAW, "r") as handle:

    for line_number, line in enumerate(handle, start=1):

        line = line.strip()

        if not line:
            continue

        fields = line.split("\t")

        if len(fields) != len(OUTFMT_FIELDS):
            raise ValueError(
                f"Unexpected number of columns in DIAMOND output "
                f"at line {line_number}: expected "
                f"{len(OUTFMT_FIELDS)}, found {len(fields)}"
            )

        # DIAMOND outfmt 6 is headerless.
        row = dict(zip(OUTFMT_FIELDS, fields))

        query_id = row["qseqid"]
        subject_id = row["sseqid"]

        if query_id not in bs_lengths:
            raise ValueError(
                f"Forward query ID not found in "
                f"B. sorokiniana FASTA: {query_id}"
            )

        if subject_id not in phi_lengths:
            raise ValueError(
                f"Forward subject ID not found in "
                f"PHI-base FASTA: {subject_id}"
            )

        qlen = bs_lengths[query_id]
        slen = phi_lengths[subject_id]

        qcov = alignment_coverage(
            row["qstart"],
            row["qend"],
            qlen
        )

        scov = alignment_coverage(
            row["sstart"],
            row["send"],
            slen
        )

        output_row = {
            "query": query_id,
            "sseqid": subject_id,
            "pident": row["pident"],
            "alignment_length": row["length"],
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "qstart": row["qstart"],
            "qend": row["qend"],
            "sstart": row["sstart"],
            "send": row["send"],
            "qlen": qlen,
            "slen": slen,
            "query_coverage": round(qcov, 4),
            "subject_coverage": round(scov, 4),
            "strain": strain_from_protein_id(query_id),
        }

        forward_rows.append(output_row)


FORWARD_FIELDS = [
    "query",
    "sseqid",
    "pident",
    "alignment_length",
    "evalue",
    "bitscore",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qlen",
    "slen",
    "query_coverage",
    "subject_coverage",
    "strain",
]

FORWARD_COVERAGE_FILE = (
    FORWARD_DIR / "forward_all_hits_with_coverage.tsv"
)

write_tsv(
    forward_rows,
    FORWARD_COVERAGE_FILE,
    FORWARD_FIELDS
)

print(f"Forward alignments processed: {len(forward_rows)}")
print(f"Output: {FORWARD_COVERAGE_FILE}")


# ============================================================
# 9. FILTER FORWARD HITS
# ============================================================

def select_best_qualifying_hits(
    rows,
    query_cov_min,
    subject_cov_min
):
    """
    Filter alignments and select the highest-bitscore qualifying
    PHI-base match for each B. sorokiniana query.

    Multiple qualifying matches are retained in a separate output.
    """

    passing = [
        row for row in rows
        if passes_filters(
            row,
            query_cov_min,
            subject_cov_min
        )
    ]

    best_by_query = {}

    for row in passing:

        query_id = row["query"]

        if query_id not in best_by_query:

            best_by_query[query_id] = row

        else:

            old_score = float(
                best_by_query[query_id]["bitscore"]
            )

            new_score = float(row["bitscore"])

            if new_score > old_score:
                best_by_query[query_id] = row

    best_rows = list(best_by_query.values())

    best_rows.sort(
        key=lambda x: (
            x["query"],
            -float(x["bitscore"])
        )
    )

    return passing, best_rows


ANALYSIS_RESULTS = {}

for qcov_min, scov_min in COVERAGE_THRESHOLDS:

    label = f"{int(qcov_min)}_{int(scov_min)}"

    passing_rows, best_rows = select_best_qualifying_hits(
        forward_rows,
        qcov_min,
        scov_min
    )

    passing_file = (
        FORWARD_DIR /
        f"forward_passing_{label}_coverage.tsv"
    )

    best_file = (
        FORWARD_DIR /
        f"forward_best_{label}_coverage.tsv"
    )

    write_tsv(
        passing_rows,
        passing_file,
        FORWARD_FIELDS
    )

    write_tsv(
        best_rows,
        best_file,
        FORWARD_FIELDS
    )

    ANALYSIS_RESULTS[label] = {
        "passing": passing_rows,
        "best": best_rows,
        "query_coverage": qcov_min,
        "subject_coverage": scov_min,
    }

    print("\n" + "-" * 60)
    print(f"Coverage analysis: query >= {qcov_min}%, "
          f"subject >= {scov_min}%")
    print(f"Passing alignment records: {len(passing_rows)}")
    print(f"Unique query proteins: {len(best_rows)}")
    print(
        "Unique orthogroups:",
        len({
            row["query"].split("|", 1)[0]
            for row in best_rows
        })
    )
    print(f"All passing hits: {passing_file}")
    print(f"Best qualifying hit per query: {best_file}")


# ============================================================
# 10. PREPARE REVERSE SEARCH QUERIES
# ============================================================

print("\nPreparing strain-specific reverse searches...")

# Collect PHI-base reference proteins that occur in either
# coverage-filtered best-hit set.

reverse_queries_by_strain = defaultdict(set)

for result in ANALYSIS_RESULTS.values():

    for row in result["best"]:

        strain = row["strain"]
        phi_id = row["sseqid"]

        if phi_id not in phi_proteins:
            raise ValueError(
                f"PHI-base reference absent from FASTA: {phi_id}"
            )

        reverse_queries_by_strain[strain].add(phi_id)


# ============================================================
# 11. REVERSE DIAMOND SEARCH
# ============================================================

print("\n" + "=" * 65)
print("REVERSE SEARCH: PHI-base -> SAME-STRAIN PROTEOME")
print("=" * 65)

reverse_best_hits = {}

for strain in sorted(reverse_queries_by_strain):

    reference_ids = sorted(
        reverse_queries_by_strain[strain]
    )

    if not reference_ids:
        continue

    print(
        f"\nReverse search for {strain}: "
        f"{len(reference_ids)} PHI-base proteins"
    )

    subset_sequences = {
        phi_id: phi_proteins[phi_id]
        for phi_id in reference_ids
    }

    reverse_query_fasta = (
        REVERSE_DIR / f"{strain}_PHI_queries.fasta"
    )

    write_fasta(
        subset_sequences,
        reverse_query_fasta
    )

    reverse_output = (
        REVERSE_DIR / f"{strain}_reverse_hits.tsv"
    )

    run_command([
        "diamond",
        "blastp",
        "--query", reverse_query_fasta,
        "--db", STRAIN_DBS[strain],
        "--out", reverse_output,
        "--outfmt", "6",
        *OUTFMT_FIELDS,
        "--evalue", EVALUE_MAX,
        "--max-target-seqs", 1,
        "--sensitive",
        "--threads", THREADS,
    ])

    reverse_rows = read_tsv(reverse_output)

    for row in reverse_rows:

        phi_id = row["qseqid"]
        bs_id = row["sseqid"]

        if phi_id not in phi_lengths:
            raise ValueError(
                f"Reverse query absent from PHI-base: {phi_id}"
            )

        if bs_id not in bs_lengths:
            raise ValueError(
                f"Reverse target absent from B. sorokiniana FASTA: "
                f"{bs_id}"
            )

        # Coverage in reverse direction:
        # query = PHI-base protein
        # subject = B. sorokiniana protein

        reverse_qcov = alignment_coverage(
            row["qstart"],
            row["qend"],
            phi_lengths[phi_id]
        )

        reverse_scov = alignment_coverage(
            row["sstart"],
            row["send"],
            bs_lengths[bs_id]
        )

        reverse_best_hits[(strain, phi_id)] = {
            "reverse_best_hit": bs_id,
            "reverse_pident": float(row["pident"]),
            "reverse_evalue": float(row["evalue"]),
            "reverse_bitscore": float(row["bitscore"]),
            "reverse_query_coverage": reverse_qcov,
            "reverse_subject_coverage": reverse_scov,
            "reverse_alignment_length": row["length"],
        }


# ============================================================
# 12. CLASSIFY RBH PAIRS
# ============================================================

print("\nClassifying reciprocal best hits...")

RBH_FIELDS = [
    "query",
    "strain",
    "phi_base_reference",
    "forward_identity",
    "forward_evalue",
    "forward_bitscore",
    "forward_query_coverage",
    "forward_subject_coverage",
    "reverse_best_hit",
    "reverse_identity",
    "reverse_evalue",
    "reverse_bitscore",
    "reverse_query_coverage",
    "reverse_subject_coverage",
    "reverse_passes_filters",
    "is_RBH",
    "classification",
]


def classify_rbh(rows, qcov_min, scov_min):

    output_rows = []

    for row in rows:

        query_id = row["query"]
        strain = row["strain"]
        phi_id = row["sseqid"]

        reverse = reverse_best_hits.get(
            (strain, phi_id)
        )

        if reverse is None:

            output_rows.append({
                "query": query_id,
                "strain": strain,
                "phi_base_reference": phi_id,
                "forward_identity": row["pident"],
                "forward_evalue": row["evalue"],
                "forward_bitscore": row["bitscore"],
                "forward_query_coverage": row["query_coverage"],
                "forward_subject_coverage": row["subject_coverage"],
                "reverse_best_hit": "NO_REVERSE_HIT",
                "reverse_identity": "",
                "reverse_evalue": "",
                "reverse_bitscore": "",
                "reverse_query_coverage": "",
                "reverse_subject_coverage": "",
                "reverse_passes_filters": "NO",
                "is_RBH": "NO",
                "classification": "NO_REVERSE_HIT",
            })

            continue

        reverse_passes = (
            reverse["reverse_pident"] >= IDENTITY_MIN
            and reverse["reverse_evalue"] <= EVALUE_MAX
            and reverse["reverse_query_coverage"] >= scov_min
            and reverse["reverse_subject_coverage"] >= qcov_min
        )

        # A strict RBH requires the reverse best hit to be
        # the original query protein, and the reverse alignment
        # must satisfy the same identity, E-value and coverage
        # criteria.

        is_rbh = (
            reverse["reverse_best_hit"] == query_id
            and reverse_passes
        )

        if is_rbh:
            classification = "RBH_SUPPORTED_CANDIDATE_ORTHOLOGUE"
        else:
            classification = "PHI_BASE_HOMOLOGUE_NON_RBH"

        output_rows.append({
            "query": query_id,
            "strain": strain,
            "phi_base_reference": phi_id,
            "forward_identity": row["pident"],
            "forward_evalue": row["evalue"],
            "forward_bitscore": row["bitscore"],
            "forward_query_coverage": row["query_coverage"],
            "forward_subject_coverage": row["subject_coverage"],
            "reverse_best_hit": reverse["reverse_best_hit"],
            "reverse_identity": reverse["reverse_pident"],
            "reverse_evalue": reverse["reverse_evalue"],
            "reverse_bitscore": reverse["reverse_bitscore"],
            "reverse_query_coverage": round(
                reverse["reverse_query_coverage"], 4
            ),
            "reverse_subject_coverage": round(
                reverse["reverse_subject_coverage"], 4
            ),
            "reverse_passes_filters": (
                "YES" if reverse_passes else "NO"
            ),
            "is_RBH": "YES" if is_rbh else "NO",
            "classification": classification,
        })

    return output_rows


for label, result in ANALYSIS_RESULTS.items():

    qcov_min = result["query_coverage"]
    scov_min = result["subject_coverage"]

    rbh_rows = classify_rbh(
        result["best"],
        qcov_min,
        scov_min
    )

    rbh_file = (
        OUTPUT_DIR / f"PHI_RBH_results_{label}_coverage.tsv"
    )

    write_tsv(
        rbh_rows,
        rbh_file,
        RBH_FIELDS
    )

    rbh_count = sum(
        row["is_RBH"] == "YES"
        for row in rbh_rows
    )

    non_rbh_count = sum(
        row["is_RBH"] == "NO"
        for row in rbh_rows
    )

    print("\n" + "-" * 60)
    print(f"RBH results for {label}%/{label}% coverage")
    print(f"Best qualifying forward hits: {len(rbh_rows)}")
    print(f"RBH-supported pairs: {rbh_count}")
    print(f"Non-RBH or no-reverse-hit pairs: {non_rbh_count}")
    print(f"Output: {rbh_file}")


# ============================================================
# 13. WRITE SUMMARY
# ============================================================

summary_rows = []

for label, result in ANALYSIS_RESULTS.items():

    passing = result["passing"]
    best = result["best"]

    rbh_file = (
        OUTPUT_DIR / f"PHI_RBH_results_{label}_coverage.tsv"
    )

    rbh_rows = read_tsv(rbh_file)

    rbh_count = sum(
        row["is_RBH"] == "YES"
        for row in rbh_rows
    )

    unique_orthogroups = {
        row["query"].split("|", 1)[0]
        for row in best
    }

    summary_rows.append({
        "coverage_threshold": label,
        "identity_min_percent": IDENTITY_MIN,
        "evalue_max": EVALUE_MAX,
        "query_coverage_min_percent": result["query_coverage"],
        "subject_coverage_min_percent": result["subject_coverage"],
        "passing_alignment_records": len(passing),
        "unique_query_proteins": len(best),
        "unique_orthogroups": len(unique_orthogroups),
        "RBH_supported_pairs": rbh_count,
        "non_RBH_or_no_reverse_hit": len(rbh_rows) - rbh_count,
    })


SUMMARY_FIELDS = [
    "coverage_threshold",
    "identity_min_percent",
    "evalue_max",
    "query_coverage_min_percent",
    "subject_coverage_min_percent",
    "passing_alignment_records",
    "unique_query_proteins",
    "unique_orthogroups",
    "RBH_supported_pairs",
    "non_RBH_or_no_reverse_hit",
]

SUMMARY_FILE = OUTPUT_DIR / "PHI_RBH_summary.tsv"

write_tsv(
    summary_rows,
    SUMMARY_FILE,
    SUMMARY_FIELDS
)


# ============================================================
# 14. FINAL REPORT
# ============================================================

print("\n" + "=" * 65)
print("PIPELINE COMPLETED")
print("=" * 65)

print(f"\nOutput directory:\n{OUTPUT_DIR}")

print("\nImportant output files:")

print(
    "1. Forward alignments with coverage:\n   "
    f"{FORWARD_COVERAGE_FILE}"
)

print(
    "2. Primary-filtered matches (40%/40%):\n   "
    f"{FORWARD_DIR / 'forward_passing_40_40_coverage.tsv'}"
)

print(
    "3. Best qualifying hit per query (40%/40%):\n   "
    f"{FORWARD_DIR / 'forward_best_40_40_coverage.tsv'}"
)

print(
    "4. Primary RBH results:\n   "
    f"{OUTPUT_DIR / 'PHI_RBH_results_40_40_coverage.tsv'}"
)

print(
    "5. Stricter RBH results (50%/50%):\n   "
    f"{OUTPUT_DIR / 'PHI_RBH_results_50_50_coverage.tsv'}"
)

print(
    "6. Summary:\n   "
    f"{SUMMARY_FILE}"
)

print(
    "\nInterpretation reminder: RBH supports candidate orthology "
    "but does not prove functional equivalence or establish "
    "the associated phenotype in B. sorokiniana."
)
