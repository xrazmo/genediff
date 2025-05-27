#!/bin/bash

# GeneDiff Example Analysis Script
# This script demonstrates how to use GeneDiff with sample data

set -e  # Exit on any error

echo "============================================"
echo "GeneDiff Example Analysis"
echo "============================================"

SCRIPT_DIR="scripts"
# Check if we're in the correct directory
if [ ! -f "$SCRIPT_DIR/gene_diff.py" ]; then
    echo "Error: gene_diff.py not found. Please run this script from the genediff repository root."
    exit 1
fi

# Check if example files exist
EXAMPLE_DIR="examples"
REF_FILE="$EXAMPLE_DIR/sample_reference.fna"
QUERY_FILE="$EXAMPLE_DIR/sample_query.ffn"
GENES_FILE="$EXAMPLE_DIR/target_genes.txt"

if [ ! -f "$REF_FILE" ]; then
    echo "Error: $REF_FILE not found. Please ensure example files are present."
    exit 1
fi

if [ ! -f "$QUERY_FILE" ]; then
    echo "Error: $QUERY_FILE not found. Please ensure example files are present."
    exit 1
fi

# Create output directory
OUTPUT_DIR="example_results"
mkdir -p "$OUTPUT_DIR"

echo "Input files:"
echo "  Reference: $REF_FILE"
echo "  Query: $QUERY_FILE"
if [ -f "$GENES_FILE" ]; then
    echo "  Target genes: $GENES_FILE"
fi
echo "  Output directory: $OUTPUT_DIR"
echo ""

# Check if conda environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "genediff_env" ]]; then
    echo "Warning: genediff_env conda environment not detected."
    echo "Please activate with: conda activate genediff_env"
    echo ""
fi

# Check required tools
echo "Checking required tools..."
MISSING_TOOLS=()

if ! command -v makeblastdb &> /dev/null; then
    MISSING_TOOLS+=("makeblastdb")
fi

if ! command -v blastn &> /dev/null; then
    MISSING_TOOLS+=("blastn")
fi

if ! command -v mafft &> /dev/null; then
    MISSING_TOOLS+=("mafft")
fi

if [ ${#MISSING_TOOLS[@]} -ne 0 ]; then
    echo "Error: Missing required tools: ${MISSING_TOOLS[*]}"
    echo "Please install BLAST+ and MAFFT tools."
    exit 1
fi

echo "✓ All required tools found"
echo ""

# Show file information
echo "File information:"
echo "Reference file:"
grep -c "^>" "$REF_FILE" | xargs echo "  Number of sequences:"
echo "Query file:"
grep -c "^>" "$QUERY_FILE" | xargs echo "  Number of sequences:"

if [ -f "$GENES_FILE" ]; then
    echo "Target genes file:"
    wc -l < "$GENES_FILE" | xargs echo "  Number of target genes:"
fi
echo ""

# Run GeneDiff analysis
echo "Running GeneDiff analysis..."
echo "Command: python gene_diff.py -r $REF_FILE -q $QUERY_FILE -o $OUTPUT_DIR"

if [ -f "$GENES_FILE" ]; then
    echo "         (with gene filtering: -g $GENES_FILE)"
    python $SCRIPT_DIR/gene_diff.py -r "$REF_FILE" -q "$QUERY_FILE" -o "$OUTPUT_DIR" -g "$GENES_FILE" -v
else
    echo "         (analyzing all genes)"
    python $SCRIPT_DIR/gene_diff.py -r "$REF_FILE" -q "$QUERY_FILE" -o "$OUTPUT_DIR" -v
fi

echo ""
echo "============================================"
echo "Analysis Complete!"
echo "============================================"

# Show results
echo "Output files generated:"
for file in "$OUTPUT_DIR"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        size=$(ls -lh "$file" | awk '{print $5}')
        echo "  $filename ($size)"
    fi
done

echo ""
echo "Results summary:"

# Show summary from main results file
SUMMARY_FILE=$(find "$OUTPUT_DIR" -name "*.pwBlastn.genediff.csv" | head -1)
if [ -f "$SUMMARY_FILE" ]; then
    echo "  Genes analyzed: $(tail -n +2 "$SUMMARY_FILE" | wc -l)"
    echo "  Genes with mutations: $(tail -n +2 "$SUMMARY_FILE" | awk -F',' '$8 > 0' | wc -l)"
    
    # Show top mutated genes
    echo "  Top 5 most mutated genes:"
    tail -n +2 "$SUMMARY_FILE" | sort -t',' -k8 -nr | head -5 | while IFS=',' read -r gene_id ref_id query_id identity qcov scov evalue mutations rest; do
        echo "    $gene_id: $mutations mutations (${identity}% identity)"
    done
fi

echo ""
echo "Example analysis completed successfully!"
echo "Explore the results in the '$OUTPUT_DIR' directory."
echo "Open the .csv files in Excel or any spreadsheet application for detailed analysis."