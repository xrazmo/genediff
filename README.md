# GeneDiff

A comprehensive bioinformatics tool for gene-level mutation comparison between reference and query multi-FASTA files using BLAST searches, MAFFT alignments, and detailed mutation analysis.

## Features

- **BLAST-based gene matching**: Identifies homologous genes between reference and query sequences
- **Multiple sequence alignment**: Uses MAFFT for accurate sequence alignment
- **Comprehensive mutation analysis**: Detects substitutions, insertions, and deletions
- **Protein effect prediction**: Analyzes missense, nonsense, frameshift, and synonymous mutations
- **Flexible gene filtering**: Optional CSV-based gene selection
- **Detailed reporting**: Multiple output formats including summary tables and detailed mutation lists

## Requirements

### Software Dependencies
- Python ≥ 3.8
- BLAST+ tools (makeblastdb, blastn)
- MAFFT multiple sequence alignment tool

### Python Dependencies
- biopython
- pandas
- All other dependencies are part of Python standard library

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/xrazmo/genediff.git
cd genediff
```

### 2. Create Conda Environment
```bash
conda env create -f environment.yml
conda activate genediff_env
```

### 3. Install External Tools

#### BLAST+ (if not installed via conda)
```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# MacOS with Homebrew
brew install blast

# Or download from NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
```

#### MAFFT (if not installed via conda)
```bash
# Ubuntu/Debian
sudo apt-get install mafft

# MacOS with Homebrew
brew install mafft

# Or download from: https://mafft.cbrc.jp/alignment/software/
```

## Usage

### Basic Usage
```bash
python ./scripts/gene_diff.py -r reference.fna -q query.fna -o output_directory
```

### Advanced Usage with Gene Filtering
```bash
python ./scripts/gene_diff.py \
    -r reference_genome.fna \
    -q query_genome.fna \
    -o results/ \
    -g target_genes.csv \
    -v
```

### Parameters

| Parameter | Description | Required |
|-----------|-------------|----------|
| `-r, --reference` | Reference multi-FASTA file (.fna) | Yes |
| `-q, --query` | Query multi-FASTA file (.fna) | Yes |
| `-o, --output` | Output directory for results | Yes |
| `-g, --genes` | CSV file with target gene names (optional) | No |
| `-v, --verbose` | Enable verbose logging | No |

## Input File Formats

### FASTA Files
- Multi-FASTA format with CDS sequences
- Headers should contain gene annotations in standard format:
  ```
  >gene_id [gene=geneName] [locus_tag=locusTag] [protein=proteinDescription] [protein_id=proteinID]
  ```

### Gene Filter CSV (Optional)
- Simple CSV file with one gene name per row
- First column should contain gene names matching those in reference FASTA headers
- Example:
  ```
  gyrA
  rpoB
  katG
  ```

## Output Files

The tool generates three main output files:

### 1. Summary Results (`{query_name}.pwBlastn.genediff.csv`)
Contains per-gene comparison summary with:
- Gene identification information
- BLAST alignment statistics
- Mutation counts
- Protein effect predictions

### 2. Detailed Mutations (`{query_name}.detailed_mutations.genediff.csv`)
Lists all individual mutations with:
- Mutation type (substitution, insertion, deletion)
- Genomic positions
- Nucleotide changes
- Standard mutation notation

### 3. Full Results (`{query_name}.full_results.genediff.json`)
Complete analysis results in JSON format for programmatic access

## Example Workflow

1. **Prepare input files**:
   ```bash
   # Reference genome CDS sequences
   reference_genome.fna
   
   # Query genome CDS sequences  
   query_genome.fna
   
   # Optional: genes of interest
   target_genes.csv
   ```

2. **Run analysis**:
   ```bash
   python gene_diff.py \
       -r reference_genome.fna \
       -q query_genome.fna \
       -o analysis_results/ \
       -g target_genes.csv \
       -v
   ```

3. **Review results**:
   ```bash
   ls analysis_results/
   # query_genome.pwBlastn.genediff.csv
   # query_genome.detailed_mutations.genediff.csv  
   # query_genome.full_results.genediff.json
   ```

## Algorithm Overview

1. **Gene Loading**: Parse reference and query FASTA files, extract gene annotations
2. **BLAST Database Creation**: Build nucleotide BLAST database from query sequences
3. **Homolog Identification**: BLAST each reference gene against query database
4. **Quality Filtering**: Filter hits by identity (>50%) and coverage (>10%)
5. **Sequence Alignment**: Perform pairwise MAFFT alignment of matched sequences
6. **Mutation Detection**: Identify all nucleotide-level differences
7. **Protein Analysis**: Translate sequences and predict protein-level effects
8. **Result Export**: Generate comprehensive reports in multiple formats

## Performance Considerations

- Processing time scales with number of genes and sequence lengths
- BLAST database creation is performed once per analysis
- Temporary files are automatically cleaned up
- Memory usage depends on largest individual gene sequences

## Troubleshooting

### Common Issues

1. **Missing external tools**:
   ```
   Error: Missing required tools: makeblastdb, blastn, mafft
   ```
   Solution: Install BLAST+ and MAFFT tools

2. **No BLAST hits found**:
   - Check sequence similarity between reference and query
   - Verify FASTA file formats and gene annotations
   - Consider adjusting BLAST parameters in code if needed

3. **Gene annotation parsing issues**:
   - Ensure FASTA headers follow standard format with bracketed annotations
   - Check gene names in filter CSV match those in reference file

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
GeneDiff
https://github.com/xrazmo/genediff
```

## Contact

- **Issues**: https://github.com/xrazmo/genediff/issues

## Changelog

### v1.0.0 (Initial Release)
- Core mutation comparison functionality
- BLAST-based gene matching
- MAFFT sequence alignment
- Comprehensive mutation analysis
- Multiple output formats
- Gene filtering support
