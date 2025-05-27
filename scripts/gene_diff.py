#!/usr/bin/env python3
"""
Gene-level mutation comparison between reference and query multi-FASTA files.

This script performs BLAST searches, MAFFT alignments, and mutation analysis
to compare genes between reference and query sequences.

Requirements:
- Biopython
- BLAST+ tools (makeblastdb, blastn)
- MAFFT
- pandas
"""

import os
import sys
import re
import csv
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio.Blast import NCBIXML
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneComparator:
    """Main class for gene-level mutation comparison."""
    
    def __init__(self, reference_fasta: str, query_fasta: str, output_dir: str, 
                 gene_csv: Optional[str] = None):
        """
        Initialize the GeneComparator.
        
        Args:
            reference_fasta: Path to reference multi-FASTA file
            query_fasta: Path to query multi-FASTA file
            output_dir: Directory to save results
            gene_csv: Optional CSV file with gene names to filter
        """
        self.reference_fasta = Path(reference_fasta)
        self.query_fasta = Path(query_fasta)
        self.output_dir = Path(output_dir)
        self.gene_csv = Path(gene_csv) if gene_csv else None
        self.query_name = os.path.basename(query_fasta).split('.')[0]
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize data structures
        self.reference_genes = {}
        self.query_genes = {}
        self.target_genes = set()
        self.results = []
        
    def parse_fasta_header(self, header: str) -> Dict[str, str]:
        """
        Parse FASTA header to extract gene information.
        
        Args:
            header: FASTA header string
            
        Returns:
            Dictionary containing parsed information
        """
        info = {}
        
        # Extract gene name
        gene_match = re.search(r'\[gene=([^\]]+)\]', header)
        if gene_match:
            info['gene'] = gene_match.group(1)
            
        # Extract locus tag
        locus_match = re.search(r'\[locus_tag=([^\]]+)\]', header)
        if locus_match:
            info['locus_tag'] = locus_match.group(1)
            
        # Extract protein info
        protein_match = re.search(r'\[protein=([^\]]+)\]', header)
        if protein_match:
            info['protein'] = protein_match.group(1)
            
        # Extract protein ID
        protein_id_match = re.search(r'\[protein_id=([^\]]+)\]', header)
        if protein_id_match:
            info['protein_id'] = protein_id_match.group(1)
            
        return info
        
    def load_target_genes(self) -> Set[str]:
        """
        Load target gene names from CSV file if provided.
        
        Returns:
            Set of target gene names
        """
        if not self.gene_csv or not self.gene_csv.exists():
            return set()
            
        target_genes = set()
        try:
            with open(self.gene_csv, 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    if row:  # Skip empty rows
                        target_genes.add(row[0].strip())
            logger.info(f"Loaded {len(target_genes)} target genes from CSV")
        except Exception as e:
            logger.error(f"Error loading gene CSV: {e}")
            
        return target_genes
        
    def load_reference_genes(self) -> Dict[str, SeqRecord]:
        """
        Load and parse reference genes from FASTA file.
        
        Returns:
            Dictionary of gene records keyed by gene name/locus_tag
        """
        reference_genes = {}
        target_genes = self.load_target_genes()
        
        try:
            for record in SeqIO.parse(self.reference_fasta, "fasta"):
                header_info = self.parse_fasta_header(record.description)
                
                # Determine gene identifier
                gene_name = header_info.get('gene')
                locus_tag = header_info.get('locus_tag')
                
                identifier = gene_name if gene_name else locus_tag
                
                if not identifier:
                    logger.warning(f"No gene name or locus tag found for: {record.id}")
                    continue
                    
                # Filter by target genes if CSV provided
                if target_genes and gene_name not in target_genes:
                    continue
                    
                # Store additional info in record
                record.gene_info = header_info
                record.identifier = identifier
                
                reference_genes[identifier] = record
                
            logger.info(f"Loaded {len(reference_genes)} reference genes")
            
        except Exception as e:
            logger.error(f"Error loading reference FASTA: {e}")
            sys.exit(1)
            
        return reference_genes
        
    def load_query_genes(self) -> Dict[str, SeqRecord]:
        """
        Load query genes from FASTA file.
        
        Returns:
            Dictionary of query gene records
        """
        query_genes = {}
        
        try:
            for record in SeqIO.parse(self.query_fasta, "fasta"):
                query_genes[record.id] = record
                
            logger.info(f"Loaded {len(query_genes)} query genes")
            
        except Exception as e:
            logger.error(f"Error loading query FASTA: {e}")
            sys.exit(1)
            
        return query_genes
        
    def create_blast_database(self) -> str:
        """
        Create BLAST database from query FASTA file.
        
        Returns:
            Path to BLAST database
        """
        db_dir = os.path.join(self.output_dir, ".tmp_blastdb")
        os.makedirs(db_dir, exist_ok=True)
        db_path = os.path.join(db_dir, "query_blast_db")
  
        try:
            makedb_cmd = NcbimakeblastdbCommandline(
                cmd='makeblastdb',
                input_file=str(self.query_fasta),
                dbtype="nucl",
                out=str(db_path)
            )
            
            stdout, stderr = makedb_cmd()
            
            if stderr:
                logger.warning(f"BLAST database creation warning: {stderr}")
                
            logger.info("BLAST database created successfully")
            return str(db_path)
            
        except Exception as e:
            logger.error(f"Error creating BLAST database: {e}")
            sys.exit(1)
            
    def run_blast_search(self, query_sequence: SeqRecord, blast_db: str) -> Optional[Dict]:
        """
        Run BLAST search for a single gene against the database.
        
        Args:
            query_sequence: Reference gene sequence to search
            blast_db: Path to BLAST database
            
        Returns:
            Best BLAST hit information or None
        """
        # Create temporary query file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_query:
            SeqIO.write(query_sequence, temp_query, "fasta")
            temp_query_path = temp_query.name
        print(temp_query_path)
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.b6', delete=False) as temp_output:
            temp_output_path = temp_output.name
            
        try:
            # Run BLAST
            blast_cmd = NcbiblastnCommandline(
                cmd='blastn',
                query=temp_query_path,
                db=blast_db,
                evalue=0.001,
                outfmt="6 qseqid sseqid pident qlen slen length stitle qcovs evalue bitscore",  
                out=temp_output_path,
                max_target_seqs=1
            )
            print(temp_output_path)
            stdout, stderr = blast_cmd()
            # Parse BLAST results
            columns = ['qseqid', 'sseqid', 'pident', 'qlen', 'slen', 'length', 'stitle', 'qcovs','evalue','bitscore']
            b6_df = pd.read_csv(temp_output_path, sep='\t', names=columns)
            b6_df['scovs'] = b6_df['length'] / b6_df['slen'] * 100
            # the subject gene at least should cover 10% of the true reference gene.
            filtered = b6_df[(b6_df['qcovs'] > 10) & (b6_df['pident'] > 50)]
            if not filtered.empty:
                hit = filtered.iloc[0]
                return {
                    'hit_id': hit['sseqid'],
                    'identity_pct': round(float(hit['pident']),2),
                    'qcoverage_pct': round(float(hit['scovs']),2), # here the query sequence is actually the reference gene, reversing the coverage names for final output
                    'scoverage_pct': round(float(hit['qcovs']),2), # same here
                    'e_value': hit['evalue'],
                    'bit_score': hit['bitscore'],
                    'alignment_length': int(hit['length'])
                }
            
        except Exception as e:
            logger.error(f"Error running BLAST for {query_sequence.id}: {e}")
            
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_query_path)
                os.unlink(temp_output_path)
            except:
                pass
                
        return None
        
    def run_mafft_alignment(self, seq1: SeqRecord, seq2: SeqRecord) -> Tuple[str, str]:
        """
        Perform MAFFT alignment between two sequences.
        
        Args:
            seq1: First sequence (reference)
            seq2: Second sequence (query match)
            
        Returns:
            Tuple of aligned sequences
        """
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
            SeqIO.write([seq1, seq2], temp_input, "fasta")
            temp_input_path = temp_input.name
            
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='r', suffix='.fasta', delete=False) as temp_output:
            temp_output_path = temp_output.name
            
        try:
            # Run MAFFT
            cmd = ['mafft', '--auto', '--quiet', temp_input_path]
            
            with open(temp_output_path, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, 
                                      text=True, check=True)
            
            # Parse aligned sequences
            aligned_seqs = list(SeqIO.parse(temp_output_path, "fasta"))
            
            if len(aligned_seqs) == 2:
                return str(aligned_seqs[0].seq), str(aligned_seqs[1].seq)
            else:
                raise ValueError("MAFFT did not return expected number of sequences")
                
        except Exception as e:
            logger.error(f"Error running MAFFT alignment: {e}")
            return str(seq1.seq), str(seq2.seq)  # Return unaligned as fallback
            
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_input_path)
                os.unlink(temp_output_path)
            except:
                pass
                
    def extract_mutations(self, aligned_ref: str, aligned_query: str) -> List[Dict]:
        """
        Extract nucleotide-level mutations from aligned sequences.
        
        Args:
            aligned_ref: Aligned reference sequence
            aligned_query: Aligned query sequence
            
        Returns:
            List of mutation dictionaries
        """
        mutations = []
        ref_pos = 0  # Position in original reference sequence (0-based)
        
        for i, (ref_nt, query_nt) in enumerate(zip(aligned_ref, aligned_query)):
            if ref_nt != '-':
                ref_pos += 1
                
            if ref_nt != query_nt:
                if ref_nt == '-':
                    # Insertion in query
                    mutations.append({
                        'type': 'insertion',
                        'position': ref_pos,
                        'ref_nt': '-',
                        'alt_nt': query_nt,
                        'change': f"ins{query_nt}"
                    })
                elif query_nt == '-':
                    # Deletion in query
                    mutations.append({
                        'type': 'deletion',
                        'position': ref_pos,
                        'ref_nt': ref_nt,
                        'alt_nt': '-',
                        'change': f"{ref_pos}{ref_nt}del"
                    })
                else:
                    # Substitution
                    mutations.append({
                        'type': 'substitution',
                        'position': ref_pos,
                        'ref_nt': ref_nt,
                        'alt_nt': query_nt,
                        'change': f"{ref_pos}{ref_nt}>{query_nt}"
                    })
                    
        return mutations
        
    def analyze_protein_effects(self, aligned_ref: str, aligned_query: str) -> Dict:
        """
        Analyze protein-level effects of mutations.
        
        Args:
            aligned_ref: Aligned reference sequence
            aligned_query: Aligned query sequence
            
        Returns:
            Dictionary with protein effect analysis
        """
        # Remove gaps and translate
        ref_clean = aligned_ref.replace('-', '')
        query_clean = aligned_query.replace('-', '')
        
        # Ensure sequences are divisible by 3
        ref_clean = ref_clean[:len(ref_clean)//3*3]
        query_clean = query_clean[:len(query_clean)//3*3]
        
        try:
            ref_protein = Seq(ref_clean).translate(table=11)
            query_protein = Seq(query_clean).translate(table=11)
            
            # Analyze protein changes
            protein_changes = []
            effect_types = set()
            
            # Check for frameshift (length difference not divisible by 3)
            if len(aligned_ref.replace('-', '')) != len(aligned_query.replace('-', '')):
                length_diff = len(aligned_query.replace('-', '')) - len(aligned_ref.replace('-', ''))
                if length_diff % 3 != 0:
                    effect_types.add('frameshift')
                    
            # Compare protein sequences
            min_len = min(len(ref_protein), len(query_protein))
            
            for i in range(min_len):
                if ref_protein[i] != query_protein[i]:
                    if query_protein[i] == '*':
                        effect_types.add('stop_gained')
                        protein_changes.append(f"{ref_protein[i]}{i+1}*")
                    elif ref_protein[i] == '*':
                        effect_types.add('stop_lost')
                        protein_changes.append(f"*{i+1}{query_protein[i]}")
                    else:
                        effect_types.add('missense')
                        protein_changes.append(f"{ref_protein[i]}{i+1}{query_protein[i]}")
                        
            # Check for synonymous mutations (same protein, different DNA)
            if ref_clean != query_clean and str(ref_protein) == str(query_protein):
                effect_types.add('synonymous')
                
            return {
                'ref_protein': str(ref_protein),
                'query_protein': str(query_protein),
                'protein_changes': protein_changes,
                'effect_types': list(effect_types)
            }
            
        except Exception as e:
            logger.warning(f"Error in protein analysis: {e}")
            return {
                'ref_protein': '',
                'query_protein': '',
                'protein_changes': [],
                'effect_types': ['translation_error']
            }
            
    def process_gene_comparison(self, gene_id: str, ref_gene: SeqRecord, 
                              blast_db: str) -> Optional[Dict]:
        """
        Process comparison for a single gene.
        
        Args:
            gene_id: Gene identifier
            ref_gene: Reference gene sequence
            blast_db: Path to BLAST database
            
        Returns:
            Comparison results dictionary
        """
        logger.info(f"Processing gene: {gene_id}")
        
        # Find best BLAST hit
        blast_hit = self.run_blast_search(ref_gene, blast_db)
        
        if not blast_hit:
            logger.warning(f"No suitable BLAST hit found for {gene_id}")
            return None
            
        # Get matching query sequence
        hit_id = blast_hit['hit_id']
        if hit_id not in self.query_genes:
            logger.warning(f"Query sequence {hit_id} not found for {gene_id}")
            return None
            
        query_gene = self.query_genes[hit_id]
        
        # Perform alignment
        aligned_ref, aligned_query = self.run_mafft_alignment(ref_gene, query_gene)
        
        # Extract mutations
        mutations = self.extract_mutations(aligned_ref, aligned_query)
        
        # Analyze protein effects
        protein_analysis = self.analyze_protein_effects(aligned_ref, aligned_query)
        
        return {
            'gene_id': gene_id,
            'ref_sequence_id': ref_gene.id,
            'query_sequence_id': hit_id,
            'blast_identity': blast_hit['identity_pct'],
            'blast_scoverage': blast_hit['scoverage_pct'],
            'blast_qcoverage': blast_hit['qcoverage_pct'],
            'blast_evalue': blast_hit['e_value'],
            'num_mutations': len(mutations),
            'mutations': mutations,
            'protein_analysis': protein_analysis,
            'gene_info': getattr(ref_gene, 'gene_info', {})
        }
        
    def run_analysis(self):
        """Run the complete gene comparison analysis."""
        logger.info("Starting gene comparison analysis")
        
        # Load sequences
        self.reference_genes = self.load_reference_genes()
        self.query_genes = self.load_query_genes()
        
        if not self.reference_genes:
            logger.error("No reference genes to process")
            return
            
        # Create BLAST database
        blast_db = self.create_blast_database()
        
        # Process each reference gene
        results = []
        
        for gene_id, ref_gene in self.reference_genes.items():
            result = self.process_gene_comparison(gene_id, ref_gene, blast_db)
            if result:
                results.append(result)
                
        logger.info(f"Completed analysis for {len(results)} genes")
        
        # Save results
        self.save_results(results)
        
    def save_results(self, results: List[Dict]):
        """
        Save analysis results to output files.
        
        Args:
            results: List of analysis results
        """
        if not results:
            logger.warning("No results to save")
            return
            
        # Create summary DataFrame
        summary_data = []
        detailed_mutations = []
        
        for result in results:
            # Summary row
            summary_row = {
                'gene_id': result['gene_id'],
                'ref_sequence_id': result['ref_sequence_id'],
                'query_sequence_id': result['query_sequence_id'],
                'blast_identity_pct': round(result['blast_identity'], 2),
                'blast_qcoverage_pct': round(result['blast_qcoverage'], 2),
                'blast_scoverage_pct': round(result['blast_scoverage'], 2),
                'blast_evalue': result['blast_evalue'],
                'num_mutations': result['num_mutations'],
                'protein_effects': ';'.join(result['protein_analysis']['effect_types']),
                'protein_changes': ';'.join(result['protein_analysis']['protein_changes']),
                'gene_name': result['gene_info'].get('gene', ''),
                'locus_tag': result['gene_info'].get('locus_tag', ''),
                'protein_description': result['gene_info'].get('protein', '')
            }
            summary_data.append(summary_row)
            
            # Detailed mutations
            for mut in result['mutations']:
                mut_row = {
                    'gene_id': result['gene_id'],
                    'mutation_type': mut['type'],
                    'position': mut['position'],
                    'ref_nucleotide': mut['ref_nt'],
                    'alt_nucleotide': mut['alt_nt'],
                    'change_notation': mut['change']
                }
                detailed_mutations.append(mut_row)
                
        # Save summary
        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(self.output_dir,f"{self.query_name}.pwBlastn.genediff.csv")
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Summary saved to: {summary_file}")
        
        # Save detailed mutations
        if detailed_mutations:
            mutations_df = pd.DataFrame(detailed_mutations)
            mutations_file = os.path.join(self.output_dir,f"{self.query_name}.detailed_mutations.genediff.csv")
            mutations_df.to_csv(mutations_file, index=False)
            logger.info(f"Detailed mutations saved to: {mutations_file}")
            
        # Save full results as JSON for reference
        import json
        results_file = os.path.join(self.output_dir,f"{self.query_name}.full_results.genediff.json")
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        logger.info(f"Full results saved to: {results_file}")

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Gene-level mutation comparison between reference and query FASTA files"
    )
    
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference multi-FASTA file (.fna)"
    )
    
    parser.add_argument(
        "-q", "--query",
        required=True,
        help="Query multi-FASTA file (.fna)"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for results"
    )
    
    parser.add_argument(
        "-g", "--genes",
        help="Optional CSV file with gene names to filter"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        
    # Check required tools
    required_tools = ['makeblastdb', 'blastn', 'mafft']
    missing_tools = []
    
    for tool in required_tools:
        try:
            help_flag = '--help' if tool == 'mafft' else '-help'
            result = subprocess.run([tool,help_flag], capture_output=True, text=True)
            # For MAFFT, non-zero exit code (e.g., 1) is acceptable if help output is produced
            if tool == 'mafft' and result.returncode != 0 and result.stderr:
                continue  # MAFFT is present, ignore non-zero exit code
            elif result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, tool)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(e)
            missing_tools.append(tool)

    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        logger.error("Please ensure BLAST+ and MAFFT are installed and in PATH")
        sys.exit(1)
        
    # Run analysis
    try:
        comparator = GeneComparator(
            reference_fasta=args.reference,
            query_fasta=args.query,
            output_dir=args.output,
            gene_csv=args.genes
        )
        
        comparator.run_analysis()
        logger.info("Analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
