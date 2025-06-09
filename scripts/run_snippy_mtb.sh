#!/bin/bash

# === DEFAULTS ===
REFERENCE=""
CORES=4
SAMPLES_FILE=""

# === USAGE FUNCTION ===
usage() {
  echo "Usage: $0 -i sample_list.txt [-r reference.fasta] [-t threads]"
  echo "  -i    Sample list file (format: sampleID R1.fastq.gz R2.fastq.gz)"
  echo "  -r    Reference genome (default: H37Rv.fasta)"
  echo "  -t    Number of threads (default: 4)"
  exit 1
}

# === PARSE ARGUMENTS ===
while getopts ":i:r:t:" opt; do
  case ${opt} in
    i ) SAMPLES_FILE="$OPTARG" ;;
    r ) REFERENCE="$OPTARG" ;;
    t ) CORES="$OPTARG" ;;
    * ) usage ;;
  esac
done

# === VALIDATE INPUT ===
if [[ -z "$SAMPLES_FILE" ]]; then
  echo "ERROR: Sample list file (-i) is required."
  usage
fi

if [[ ! -f "$REFERENCE" ]]; then
  echo "ERROR: Reference genome $REFERENCE not found."
  exit 1
fi

if [[ ! -f "$SAMPLES_FILE" ]]; then
  echo "ERROR: Sample list file $SAMPLES_FILE not found."
  exit 1
fi

# === RUN SNIPPY FOR EACH SAMPLE ===
echo "Starting Snippy runs using $SAMPLES_FILE..."

while read -r SAMPLE R1 R2; do
  OUTDIR="${SAMPLE}_snippy"
  echo "Processing $SAMPLE..."

  if [[ -d "$OUTDIR" ]]; then
    echo "Skipping $SAMPLE (output exists)"
    continue
  fi

  mkdir -p "$OUTDIR"

  snippy \
    --outdir "$OUTDIR" \
    --ref "$REFERENCE" \
    --R1 "$R1" \
    --R2 "$R2" \
    --cpus "$CORES" \
    --rgid "$SAMPLE" \
    --force > "${OUTDIR}/snippy.log" 2>&1

  if [[ $? -ne 0 ]]; then
    echo "Snippy failed for $SAMPLE. Showing last 10 lines of log:"
    tail -n 10 "${OUTDIR}/snippy.log"
    continue
  fi

done < "$SAMPLES_FILE"

# === RUN SNIPPY-CORE ===
echo "Running snippy-core..."
SNIPPY_DIRS=$(awk '{print $1"_snippy"}' "$SAMPLES_FILE" | tr '\n' ' ')

snippy-core --ref "$REFERENCE" --prefix "snippy_core" $SNIPPY_DIRS

# === OPTIONAL SNP-DISTS ===
if command -v snp-dists &> /dev/null; then
  echo "Generating SNP distance matrix..."
  snp-dists snippy_core/core.full.aln > snippy_core_snp_dists.tsv
fi

# === ORGANIZE OUTPUT FILES ===
echo "Organizing output files..."
mkdir -p results/snippy_core_output
mv snippy_core.* results/snippy_core_output/ 2>/dev/null
[[ -f snippy_core_snp_dists.tsv ]] && mv snippy_core_snp_dists.tsv results/snippy_core_output/

echo "All done. Outputs are saved in: results/snippy_core_output"
