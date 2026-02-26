#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=60-00:00:00
#SBATCH --partition=long
#SBATCH --output=slurm/%x.%j.out
# Script name:  04.Predict_eukaryome_cutoffs.sh
# Description:  Predict local and global ITS sequence similarity cutoffs for
#               each taxonomic rank using DNAbarcoder. Sequence subsets are
#               built on-the-fly from the ID files produced by
#               03.Prepare_eukaryome_subsets.sh, avoiding large intermediate
#               FASTA files on disk outside the tmp folder.
#
#               Workflow per (target x parent) combination:
#                 1. Subset global FASTA + classification to tmp/ using ID file
#                 2. Run dnabarcoder predict with -rank <target> -higherrank <parent>
#               Additionally:
#                 - A master similarity matrix is computed once from the most
#                   complete FASTA (kingdom_unique_id.txt) and reused across runs
#                 - A global (no higherrank) prediction is run for each target
#                   rank using its unique-sequence FASTA
# Note:         This script must be run from the source directory.
# =============================================================================
# PARAMETER SETUP
# =============================================================================

# --- Input files (outputs from 01 and 03 scripts) ----------------------------
readonly GLOBAL_FASTA="./data/eukaryome_ITS_v2.0.fasta"
readonly GLOBAL_CLASS="./data/eukaryome_ITS_v2.0.classification"
readonly ID_DIR="./data"                   # directory containing *_pred_id_*.txt and *_unique_id.txt

# --- Temporary directory for subset FASTA and classification files -----------
readonly TMP_DIR="./tmp"

# --- DNAbarcoder settings -----------------------------------------------------
readonly DNABARCODER="./dnabarcoder/dnabarcoder.py"
readonly OUT_DIR="./data"                  # dnabarcoder writes all outputs here (-o flag)
readonly PREFIX="eukaryome"
readonly MIN_LEN=400                       # minimum BLAST alignment length (-ml)
readonly STEP=0.001                        # similarity step size (-s)
readonly END_THRESH=1                      # upper similarity bound (-et)

# Start thresholds vary by rank: lower ranks need a higher similarity floor
declare -A START_THRESH
START_THRESH["species"]="0.8"
START_THRESH["genus"]="0.7"
START_THRESH["family"]="0.6"
START_THRESH["order"]="0.5"
START_THRESH["class"]="0.5"
START_THRESH["phylum"]="0.5"

# Parent rank combinations for each target rank (mirrors prepare script)
declare -A PARENT_RANKS
PARENT_RANKS["species"]="genus family order class phylum kingdom"
PARENT_RANKS["genus"]="family order class phylum kingdom"
PARENT_RANKS["family"]="order class phylum kingdom"
PARENT_RANKS["order"]="class phylum kingdom"
PARENT_RANKS["class"]="phylum kingdom"
PARENT_RANKS["phylum"]="kingdom"

# Filename abbreviations matching 03.Prepare_eukaryome_subsets.sh output names
declare -A RANK_ABBR
RANK_ABBR["kingdom"]="kng"
RANK_ABBR["phylum"]="phy"
RANK_ABBR["class"]="cls"
RANK_ABBR["order"]="ord"
RANK_ABBR["family"]="fam"
RANK_ABBR["genus"]="gen"
RANK_ABBR["species"]="spe"

# Target ranks to process (kingdom is excluded — no valid parent ranks)
TARGET_RANKS=("species" "genus" "family" "order" "class" "phylum")

# Create required directories
mkdir -p "$TMP_DIR" "$OUT_DIR" slurm

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# subset_fasta_and_classification <id_file> <out_fasta> <out_class>
# Filters the global FASTA and classification to the IDs listed in id_file.
# Uses Biopython for FASTA (preserves full headers) and awk for the TSV.
subset_fasta_and_classification() {
    local id_file="$1"
    local out_fasta="$2"
    local out_class="$3"

    if [[ ! -f "$id_file" ]]; then
        echo "  WARNING: ID file not found, skipping: $id_file" >&2
        return 1
    fi

    # Filter FASTA: read ID list into a set, write matching records
    python3 - <<PYEOF
from Bio import SeqIO
ids = set(open("$id_file").read().split())
records = [r for r in SeqIO.parse("$GLOBAL_FASTA", "fasta") if r.id in ids]
SeqIO.write(records, "$out_fasta", "fasta")
print(f"  Wrote {len(records)} sequences to $out_fasta")
PYEOF

    # Filter classification: keep header row + rows whose first field is in the ID list
    awk 'NR == FNR { ids[$1] = 1; next }
         FNR == 1  { print; next }
         $1 in ids { print }' "$id_file" "$GLOBAL_CLASS" > "$out_class"

    local n_class
    n_class=$(( $(wc -l < "$out_class") - 1 ))
    echo "  Wrote $n_class classification rows to $out_class"
}

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# =============================================================================
# DNABARCODER SETUP
# =============================================================================

echo ""
echo "=== SETTING UP DNABARCODER ==="
echo $(date)

if [[ -d "dnabarcoder" ]]; then
    echo "Removing existing dnabarcoder directory..."
    rm -rf dnabarcoder
fi

if ! git clone https://github.com/vuthuyduong/dnabarcoder dnabarcoder; then
    echo "ERROR: Failed to clone dnabarcoder repository" >&2
    exit 1
fi


# Clear the bundled example data shipped with the repository
echo "Clearing bundled dnabarcoder/data folder..."
mkdir -p dnabarcoder/data
rm -rf dnabarcoder/data/*

echo "DNAbarcoder cloned successfully."
echo $(date)

# =============================================================================
# MASTER SIMILARITY MATRIX
# =============================================================================
# Compute once from the most complete FASTA (all kingdom-level unique sequences)
# and reuse across all prediction runs via the -sim flag.
# =============================================================================

echo ""
echo "=== COMPUTING MASTER SIMILARITY MATRIX ==="
echo $(date)

readonly KINGDOM_FASTA="$TMP_DIR/kingdom_unique.fasta"
readonly KINGDOM_ID="$ID_DIR/kingdom_unique_id.txt"
readonly SIM_MATRIX="$OUT_DIR/${PREFIX}.sim"

echo "Subsetting global FASTA to kingdom unique sequences..."
subset_fasta_and_classification "$KINGDOM_ID" \
    "$KINGDOM_FASTA" \
    "$TMP_DIR/kingdom_unique.classification"

echo "Computing similarity matrix (this may take several hours)..."
python3 "$DNABARCODER" sim \
    -i "$KINGDOM_FASTA" \
    -o "$OUT_DIR" \
    -ml "$MIN_LEN"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Similarity matrix computation failed." >&2
    exit 1
fi

echo "Similarity matrix written to: $SIM_MATRIX"
echo $(date)

# =============================================================================
# PREDICT LOCAL SIMILARITY CUTOFFS (per target x parent combination)
# =============================================================================
# For each combination a dedicated FASTA + classification subset is built in
# tmp/, and the similarity matrix is passed via -sim so BLAST is not re-run
# for the sequences already present in the master matrix.
# =============================================================================

echo ""
echo "=== PREDICTING LOCAL SIMILARITY CUTOFFS ==="
echo $(date)

for target in "${TARGET_RANKS[@]}"; do

    st="${START_THRESH[$target]}"

    for parent in ${PARENT_RANKS[$target]}; do

        parent_abbr="${RANK_ABBR[$parent]}"
        id_file="$ID_DIR/${target}_pred_id_${parent_abbr}.txt"
        tmp_fasta="$TMP_DIR/${target}_pred_${parent_abbr}.fasta"
        tmp_class="$TMP_DIR/${target}_pred_${parent_abbr}.classification"

        # Skip if the ID file was not produced (no groups passed filters)
        if [[ ! -f "$id_file" ]]; then
            echo "  Skipping ${target} within ${parent}: ID file not found ($id_file)"
            continue
        fi

        echo ""
        echo "--- ${target} within ${parent} ---"
        echo "Subsetting sequences at: $(date)"
        subset_fasta_and_classification "$id_file" "$tmp_fasta" "$tmp_class"

        echo "Predicting cutoffs at: $(date)"
        python3 "$DNABARCODER" predict \
            -i "$tmp_fasta" \
            -c "$tmp_class" \
            -st "$st" -et "$END_THRESH" -s "$STEP" \
            -rank "$target" \
            -higherrank "$parent" \
            -prefix "$PREFIX" \
            -sim "$SIM_MATRIX" \
            -o "$OUT_DIR" \
            -ml "$MIN_LEN"

        if [[ $? -ne 0 ]]; then
            echo "  WARNING: predict failed for ${target} within ${parent}" >&2
        else
            echo "  Finished ${target} within ${parent} at: $(date)"
        fi

        # Remove subset files immediately to keep tmp/ lean
        rm -f "$tmp_fasta" "$tmp_class"

    done
done

# =============================================================================
# PREDICT GLOBAL SIMILARITY CUTOFFS (no higherrank, per target rank)
# =============================================================================
# Uses the rank-level unique-sequence FASTA (from STEP 1 of the prepare script)
# which is the most complete set for each target rank.
# =============================================================================

echo ""
echo "=== PREDICTING GLOBAL SIMILARITY CUTOFFS ==="
echo $(date)

for target in "${TARGET_RANKS[@]}"; do

    st="${START_THRESH[$target]}"
    target_abbr="${RANK_ABBR[$target]}"
    id_file="$ID_DIR/${target}_unique_id.txt"
    tmp_fasta="$TMP_DIR/${target}_unique.fasta"
    tmp_class="$TMP_DIR/${target}_unique.classification"

    if [[ ! -f "$id_file" ]]; then
        echo "  Skipping global prediction for ${target}: ID file not found ($id_file)"
        continue
    fi

    echo ""
    echo "--- global: ${target} ---"
    echo "Subsetting sequences at: $(date)"
    subset_fasta_and_classification "$id_file" "$tmp_fasta" "$tmp_class"

    echo "Predicting global cutoffs at: $(date)"
    python3 "$DNABARCODER" predict \
        -i "$tmp_fasta" \
        -c "$tmp_class" \
        -st "$st" -et "$END_THRESH" -s "$STEP" \
        -rank "$target" \
        -prefix "$PREFIX" \
        -sim "$SIM_MATRIX" \
        -o "$OUT_DIR" \
        -ml "$MIN_LEN"

    if [[ $? -ne 0 ]]; then
        echo "  WARNING: global predict failed for ${target}" >&2
    else
        echo "  Finished global prediction for ${target} at: $(date)"
    fi

    rm -f "$tmp_fasta" "$tmp_class"

done

# =============================================================================
# CLEANUP
# =============================================================================

echo ""
echo "=== CLEANUP ==="
echo "Removing tmp directory..."
rm -rf "$TMP_DIR"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo $(date)
echo ""
echo "Cutoff files written to: $OUT_DIR"
echo "  Local cutoffs  : ${PREFIX}.cutoffs.json"
echo "  Global cutoffs : ${PREFIX}.cutoffs.json (appended)"
echo "  Similarity mat : $SIM_MATRIX"

conda deactivate