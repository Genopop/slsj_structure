#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=your-account
#SBATCH --mem=16G

# --- Configuration ---
# Check if an argument is provided for IN_DIR
if [ -z "$1" ]; then
  echo "Usage: sbatch your_script_name.sh <INPUT_DIRECTORY>"
  exit 1
fi

IN_DIR="$1" # Set IN_DIR to the first command-line argument
OUT_DIR="../../data"
BASE_NAME="cartagene_all"
SUFFIX=".merged.ibd"
OUTPUT_FILE="${OUT_DIR}/${BASE_NAME}_chr1-22${SUFFIX}"

# --- Command ---
# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Concatenate the files
cat "${IN_DIR}/${BASE_NAME}"_chr{1..22}"${SUFFIX}" > "$OUTPUT_FILE"

echo "Concatenation complete. Output saved to $OUTPUT_FILE"