#!/usr/bin/env python3
"""
extract_af3_sequences.py

Extract WT and mutant protein sequences from the project reference genome and
GFF3 annotation, then generate the FASTA and JSON files needed for AlphaFold 3.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import warnings

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Suppress BioPython warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning, module="Bio")

# Standard amino acid conversion tables
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O",  # Extended amino acids
}

AA1_TO_3 = {v: k for k, v in AA3_TO_1.items()}

# Common gene ID patterns across different annotation sources
GENE_ID_PATTERNS = [
    "ID=", "gene_id=", "locus_tag=", "Name=", "gene_name=", 
    "gene=", "Gene=", "product=", "label="
]

@dataclass
class GeneInfo:
    """Gene annotation information from GFF3/GTF"""
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: Optional[str] = None
    product: Optional[str] = None
    source: str = "unknown"
    feature_type: str = "gene"

@dataclass
class VariantInfo:
    """Amino acid variant information"""
    gene: str
    mutation: str
    ref_aa: str
    position: int
    alt_aa: str
    condition: str = "unknown"
    priority: str = "medium"
    source: str = "user_input"
    confidence: str = "unspecified"
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
    
    @property
    def mutation_id(self) -> str:
        """Generate a unique mutation identifier"""
        return f"{self.gene}_{self.ref_aa}{self.position}{self.alt_aa}_{self.condition}"
    
    @property
    def hgvs_notation(self) -> str:
        """Return HGVS-style notation"""
        return f"p.{AA1_TO_3[self.ref_aa]}{self.position}{AA1_TO_3[self.alt_aa]}"

class UniversalLogger:
    """Enhanced logging with multiple output formats"""
    
    def __init__(self, log_file: Path, log_level: str = "INFO"):
        self.log_file = log_file
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.stats: Dict[str, int] = {"processed": 0, "successful": 0, "failed": 0}
        
        # Setup logging
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, encoding='utf-8'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def log(self, message: str, level: str = "INFO") -> None:
        """Log message with appropriate level"""
        level = level.upper()
        if level == "ERROR":
            self.logger.error(message)
            self.errors.append(message)
        elif level == "WARNING":
            self.logger.warning(message)
            self.warnings.append(message)
        elif level == "DEBUG":
            self.logger.debug(message)
        else:
            self.logger.info(message)

    def update_stats(self, stat: str, increment: int = 1) -> None:
        """Update processing statistics"""
        self.stats[stat] = self.stats.get(stat, 0) + increment

    def generate_summary(self) -> Dict[str, Any]:
        """Generate processing summary"""
        return {
            "timestamp": datetime.now().isoformat(),
            "statistics": self.stats.copy(),
            "errors": len(self.errors),
            "warnings": len(self.warnings),
            "error_list": self.errors[-10:],  # Last 10 errors
            "warning_list": self.warnings[-10:]  # Last 10 warnings
        }

class UniversalSequenceExtractor:
    """Universal sequence extractor for any organism/annotation format"""
    
    def __init__(self, fasta_file: Path, gff_file: Path, output_dir: Path, 
                 organism: str = "unknown", log_level: str = "INFO"):
        self.fasta_file = fasta_file
        self.gff_file = gff_file
        self.output_dir = output_dir
        self.organism = organism
        
        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "sequences").mkdir(exist_ok=True)
        (self.output_dir / "configs").mkdir(exist_ok=True)
        (self.output_dir / "reports").mkdir(exist_ok=True)
        
        # Initialize logger
        self.logger = UniversalLogger(
            self.output_dir / "extraction.log", 
            log_level
        )
        
        # Data storage
        self.genome: Dict[str, Seq] = {}
        self.genes: Dict[str, GeneInfo] = {}
        self.gene_aliases: Dict[str, str] = {}  # For flexible gene name matching

    def validate_inputs(self) -> None:
        """Validate input files and log system information"""
        if not self.fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_file}")
        if not self.gff_file.exists():
            raise FileNotFoundError(f"GFF3/GTF file not found: {self.gff_file}")
        
        # Log system information
        self.logger.log(f"=== UNIVERSAL AF3 SEQUENCE EXTRACTOR ===")
        self.logger.log(f"Organism: {self.organism}")
        self.logger.log(f"FASTA: {self.fasta_file}")
        self.logger.log(f"GFF3/GTF: {self.gff_file}")
        self.logger.log(f"Output: {self.output_dir}")
        self.logger.log(f"Timestamp: {datetime.now().isoformat()}")

    def load_genome(self) -> None:
        """Load genome sequence(s) with comprehensive validation"""
        self.logger.log("Loading genome sequences...")
        
        try:
            sequences = SeqIO.parse(str(self.fasta_file), "fasta")
            self.genome = SeqIO.to_dict(sequences)
        except Exception as e:
            raise RuntimeError(f"Failed to parse FASTA file: {e}")
        
        if not self.genome:
            raise RuntimeError("No sequences found in FASTA file")
        
        # Log genome statistics
        total_length = sum(len(seq.seq) for seq in self.genome.values())
        longest_contig = max(len(seq.seq) for seq in self.genome.values())
        
        self.logger.log(f"Loaded {len(self.genome)} sequence(s)")
        self.logger.log(f"Total genome size: {total_length:,} bp")
        self.logger.log(f"Longest contig: {longest_contig:,} bp")
        
        # Log chromosome/contig names
        contig_names = list(self.genome.keys())
        if len(contig_names) <= 10:
            self.logger.log(f"Contigs: {', '.join(contig_names)}")
        else:
            self.logger.log(f"Contigs: {', '.join(contig_names[:5])} ... ({len(contig_names)} total)")

    def detect_gff_format(self) -> str:
        """Detect GFF3 vs GTF format"""
        with open(self.gff_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if '\t' in line and len(line.split('\t')) >= 9:
                    # Check for GTF-style attributes (space-separated key-value pairs)
                    attributes = line.split('\t')[8]
                    if '; ' in attributes and '"' in attributes:
                        return "GTF"
                    else:
                        return "GFF3"
        return "GFF3"  # Default

    def parse_gff_attributes(self, attributes: str, format_type: str) -> Dict[str, str]:
        """Parse attributes section of GFF3/GTF"""
        attr_dict = {}
        
        if format_type == "GTF":
            # GTF format: gene_id "value"; transcript_id "value";
            for item in attributes.split(';'):
                item = item.strip()
                if ' ' in item:
                    key, value = item.split(' ', 1)
                    # Remove quotes
                    value = value.strip('"').strip("'")
                    attr_dict[key] = value
        else:
            # GFF3 format: ID=value;Name=value;
            for item in attributes.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    attr_dict[key] = value
        
        return attr_dict

    def extract_gene_id(self, attributes: Dict[str, str]) -> Tuple[Optional[str], Optional[str]]:
        """Extract gene ID and name using multiple strategies"""
        gene_id = None
        gene_name = None
        
        # Priority order for gene ID
        for pattern in ["locus_tag", "gene_id", "ID", "Name", "gene", "gene_name"]:
            if pattern in attributes:
                raw_id = attributes[pattern]
                # Clean common prefixes
                gene_id = raw_id.replace("gene-", "").replace("cds-", "").replace("transcript-", "")
                break
        
        # Extract gene name (often more human-readable)
        for pattern in ["Name", "gene_name", "gene", "product", "label"]:
            if pattern in attributes:
                gene_name = attributes[pattern]
                break
        
        return gene_id, gene_name

    def parse_annotations(self) -> None:
        """Parse GFF3/GTF annotations with comprehensive gene identification"""
        format_type = self.detect_gff_format()
        self.logger.log(f"Detected annotation format: {format_type}")
        
        genes = {}
        aliases = {}
        line_count = 0
        feature_counts = {"gene": 0, "CDS": 0, "exon": 0, "mRNA": 0, "other": 0}
        
        try:
            with open(self.gff_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line_count += 1
                    
                    if not line.strip() or line.startswith('#'):
                        continue
                    
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 9:
                        self.logger.log(f"Line {line_num}: Invalid format (insufficient columns)", "WARNING")
                        continue
                    
                    try:
                        seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
                        start_pos = int(start)
                        end_pos = int(end)
                    except ValueError as e:
                        self.logger.log(f"Line {line_num}: Invalid coordinates - {e}", "WARNING")
                        continue
                    
                    # Track feature types
                    if feature_type in feature_counts:
                        feature_counts[feature_type] += 1
                    else:
                        feature_counts["other"] += 1
                    
                    # Only process protein-coding features
                    if feature_type not in {"gene", "CDS", "mRNA", "transcript"}:
                        continue
                    
                    # Parse attributes
                    attr_dict = self.parse_gff_attributes(attributes, format_type)
                    gene_id, gene_name = self.extract_gene_id(attr_dict)
                    
                    if not gene_id:
                        continue
                    
                    # Create gene info
                    gene_info = GeneInfo(
                        chrom=seqid,
                        start=start_pos,
                        end=end_pos,
                        strand=strand,
                        gene_id=gene_id,
                        gene_name=gene_name,
                        product=attr_dict.get("product"),
                        source=source,
                        feature_type=feature_type
                    )
                    
                    # Prioritize CDS over gene features
                    if feature_type == "CDS" or gene_id not in genes:
                        genes[gene_id] = gene_info
                    
                    # Build aliases for flexible matching
                    aliases[gene_id] = gene_id
                    if gene_name and gene_name != gene_id:
                        aliases[gene_name] = gene_id
                    
                    # Add cleaned versions
                    for variant in [gene_id.upper(), gene_id.lower()]:
                        aliases[variant] = gene_id
                    
        except Exception as e:
            raise RuntimeError(f"Failed to parse annotation file: {e}")
        
        self.genes = genes
        self.gene_aliases = aliases
        
        # Log parsing statistics
        self.logger.log(f"Processed {line_count:,} lines")
        self.logger.log(f"Found {len(genes)} gene/CDS entries")
        self.logger.log(f"Feature counts: {dict(feature_counts)}")
        self.logger.log(f"Created {len(aliases)} gene aliases for flexible matching")

    @staticmethod
    def parse_mutation_notation(mutation: str) -> Tuple[str, int, str]:
        """Parse various mutation notation formats"""
        # Remove common prefixes
        clean_mutation = mutation.replace("p.", "").replace("c.", "").strip()
        
        # Try different patterns
        patterns = [
            r"([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})",  # Standard: Ala123Val or A123V
            r"([A-Z])(\d+)([A-Z])",                   # Single letter: A123V
            r"([A-Za-z]{3})(\d+)([A-Za-z]{3})"       # Three letter: Ala123Val
        ]
        
        for pattern in patterns:
            match = re.match(pattern, clean_mutation)
            if match:
                ref, pos, alt = match.groups()
                
                # Convert to single letter if needed
                if len(ref) > 1:
                    ref = AA3_TO_1.get(ref, ref)
                if len(alt) > 1:
                    alt = AA3_TO_1.get(alt, alt)
                
                # Validate
                if (len(ref) == 1 and len(alt) == 1 and 
                    ref in AA3_TO_1.values() and alt in AA3_TO_1.values()):
                    return ref, int(pos), alt
        
        raise ValueError(f"Unable to parse mutation format: {mutation}")

    def load_variants_from_file(self, variants_file: Path) -> List[VariantInfo]:
        """Load variants from various file formats"""
        variants = []
        file_ext = variants_file.suffix.lower()
        
        try:
            if file_ext == '.json':
                with open(variants_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                # Handle different JSON structures
                if isinstance(data, list):
                    for item in data:
                        if isinstance(item, dict):
                            variants.append(self._dict_to_variant(item))
                elif isinstance(data, dict):
                    if "variants" in data:
                        for item in data["variants"]:
                            variants.append(self._dict_to_variant(item))
                    else:
                        variants.append(self._dict_to_variant(data))
            
            else:  # CSV/TSV
                separator = '\t' if file_ext == '.tsv' else ','
                df = pd.read_csv(variants_file, sep=separator)
                
                # Flexible column mapping
                column_map = self._detect_column_mapping(df.columns)
                
                for _, row in df.iterrows():
                    variant_dict = {}
                    for standard_col, actual_col in column_map.items():
                        if actual_col:
                            variant_dict[standard_col] = str(row[actual_col]).strip()
                    
                    if variant_dict.get("gene") and variant_dict.get("mutation"):
                        variants.append(self._dict_to_variant(variant_dict))
        
        except Exception as e:
            raise RuntimeError(f"Failed to load variants from {variants_file}: {e}")
        
        self.logger.log(f"Loaded {len(variants)} variants from {variants_file}")
        return variants

    def _detect_column_mapping(self, columns: List[str]) -> Dict[str, Optional[str]]:
        """Detect column mapping for flexible input formats"""
        column_map = {
            "gene": None,
            "mutation": None,
            "condition": None,
            "priority": None,
            "source": None,
            "confidence": None
        }
        
        # Flexible column name matching
        for col in columns:
            col_lower = col.lower().strip()
            
            if any(x in col_lower for x in ["gene", "locus", "orf"]):
                column_map["gene"] = col
            elif any(x in col_lower for x in ["mutation", "variant", "change", "substitution"]):
                column_map["mutation"] = col
            elif any(x in col_lower for x in ["condition", "treatment", "sample", "group"]):
                column_map["condition"] = col
            elif any(x in col_lower for x in ["priority", "importance", "rank"]):
                column_map["priority"] = col
            elif any(x in col_lower for x in ["source", "origin", "method"]):
                column_map["source"] = col
            elif any(x in col_lower for x in ["confidence", "quality", "score"]):
                column_map["confidence"] = col
        
        return column_map

    def _dict_to_variant(self, data: Dict[str, Any]) -> VariantInfo:
        """Convert dictionary to VariantInfo object"""
        gene = str(data.get("gene", data.get("Gene", "")))
        mutation = str(data.get("mutation", data.get("Mutation", "")))
        
        if not gene or not mutation:
            raise ValueError("Missing required fields: gene and mutation")
        
        try:
            ref_aa, position, alt_aa = self.parse_mutation_notation(mutation)
        except ValueError as e:
            raise ValueError(f"Invalid mutation format for {gene}: {mutation} - {e}")
        
        return VariantInfo(
            gene=gene,
            mutation=mutation,
            ref_aa=ref_aa,
            position=position,
            alt_aa=alt_aa,
            condition=str(data.get("condition", data.get("Condition", "unknown"))),
            priority=str(data.get("priority", data.get("Priority", "medium"))),
            source=str(data.get("source", data.get("Source", "file_input"))),
            confidence=str(data.get("confidence", data.get("Confidence", "unspecified")))
        )

    def parse_command_line_variant(self, variant_string: str) -> VariantInfo:
        """Parse variant from command line format: GENE:p.Ala123Val:condition:priority"""
        parts = variant_string.split(":")
        
        if len(parts) < 2:
            raise ValueError("Variant format should be: GENE:MUTATION[:condition[:priority]]")
        
        gene = parts[0].strip()
        mutation = parts[1].strip()
        condition = parts[2].strip() if len(parts) > 2 else "unknown"
        priority = parts[3].strip() if len(parts) > 3 else "medium"
        
        try:
            ref_aa, position, alt_aa = self.parse_mutation_notation(mutation)
        except ValueError as e:
            raise ValueError(f"Invalid mutation format: {mutation} - {e}")
        
        return VariantInfo(
            gene=gene,
            mutation=mutation,
            ref_aa=ref_aa,
            position=position,
            alt_aa=alt_aa,
            condition=condition,
            priority=priority,
            source="command_line"
        )

    def find_gene(self, gene_query: str) -> Optional[GeneInfo]:
        """Find gene using flexible matching strategies"""
        # Direct match
        if gene_query in self.genes:
            return self.genes[gene_query]
        
        # Alias match
        if gene_query in self.gene_aliases:
            canonical_id = self.gene_aliases[gene_query]
            return self.genes[canonical_id]
        
        # Fuzzy matching
        query_variants = [
            gene_query,
            gene_query.upper(),
            gene_query.lower(),
            gene_query.replace("_", "-"),
            gene_query.replace("-", "_"),
        ]
        
        for variant in query_variants:
            for gene_id, gene_info in self.genes.items():
                if (variant in gene_id or gene_id in variant or
                    (gene_info.gene_name and variant in gene_info.gene_name)):
                    self.logger.log(f"Fuzzy match: '{gene_query}' -> '{gene_id}'", "INFO")
                    return gene_info
        
        return None

    def extract_protein_sequence(self, gene_info: GeneInfo) -> str:
        """Extract and translate protein sequence with validation"""
        if gene_info.chrom not in self.genome:
            raise KeyError(f"Chromosome '{gene_info.chrom}' not found in genome for gene {gene_info.gene_id}")
        
        genome_seq = self.genome[gene_info.chrom].seq
        seq_len = len(genome_seq)
        
        # Validate coordinates
        if gene_info.start < 1 or gene_info.end > seq_len:
            raise ValueError(
                f"Gene {gene_info.gene_id} coordinates ({gene_info.start}-{gene_info.end}) "
                f"exceed chromosome length ({seq_len} bp)"
            )
        
        if gene_info.start >= gene_info.end:
            raise ValueError(f"Invalid coordinates for {gene_info.gene_id}: start >= end")
        
        # Extract DNA sequence (GFF coordinates are 1-based, inclusive)
        dna_seq = genome_seq[gene_info.start - 1:gene_info.end]
        
        # Handle strand
        if gene_info.strand == "-":
            dna_seq = dna_seq.reverse_complement()
        
        # Translate to protein
        try:
            protein_seq = dna_seq.translate(to_stop=True)
            protein_str = str(protein_seq)
        except Exception as e:
            raise ValueError(f"Translation failed for {gene_info.gene_id}: {e}")
        
        # Validation
        if len(protein_str) == 0:
            raise ValueError(f"Empty protein sequence for {gene_info.gene_id}")
        
        # Check for premature stop codons
        internal_stops = protein_str.count('*')
        if internal_stops > 0:
            self.logger.log(f"Warning: {gene_info.gene_id} has {internal_stops} internal stop codon(s)", "WARNING")
        
        return protein_str.rstrip('*')  # Remove terminal stop

    def apply_mutation(self, protein_seq: str, variant: VariantInfo) -> Optional[str]:
        """Apply amino acid mutation with comprehensive validation"""
        # Check position bounds
        if variant.position < 1 or variant.position > len(protein_seq):
            self.logger.log(
                f"{variant.gene}: Position {variant.position} out of range "
                f"(protein length: {len(protein_seq)} aa)", "ERROR"
            )
            return None
        
        # Check reference amino acid
        actual_aa = protein_seq[variant.position - 1]
        if actual_aa != variant.ref_aa:
            self.logger.log(
                f"{variant.gene}: Reference mismatch at position {variant.position} "
                f"(expected {variant.ref_aa}, found {actual_aa})", "WARNING"
            )
            # Continue with observed amino acid
        
        # Apply mutation
        mutant_seq = list(protein_seq)
        mutant_seq[variant.position - 1] = variant.alt_aa
        
        return ''.join(mutant_seq)

    def write_fasta(self, filepath: Path, header: str, sequence: str) -> None:
        """Write FASTA file with validation"""
        if not sequence:
            raise ValueError("Cannot write empty sequence")
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f">{header}\n")
            # Write in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")

    def write_alphafold_config(self, filepath: Path, job_name: str, sequence: str, 
                             additional_params: Optional[Dict[str, Any]] = None) -> None:
        """Write AlphaFold 3 configuration JSON"""
        if not sequence:
            raise ValueError("Cannot create AF3 config with empty sequence")
        
        config = {
            "name": job_name,
            "modelSeeds": [42],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": sequence,
                        "count": 1
                    }
                }
            ]
        }
        
        # Add additional parameters if provided
        if additional_params:
            config.update(additional_params)
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2)

    def validate_variants(self, variants: List[VariantInfo]) -> List[VariantInfo]:
        """Validate all variants before processing"""
        valid_variants = []
        
        self.logger.log("=== VARIANT VALIDATION ===")
        
        for variant in variants:
            try:
                # Find gene
                gene_info = self.find_gene(variant.gene)
                if not gene_info:
                    self.logger.log(f"❌ {variant.gene}: Gene not found", "ERROR")
                    continue
                
                # Extract protein
                protein_seq = self.extract_protein_sequence(gene_info)
                
                # Check mutation
                mutant_seq = self.apply_mutation(protein_seq, variant)
                if mutant_seq is None:
                    self.logger.log(f"❌ {variant.gene} {variant.mutation}: Invalid mutation", "ERROR")
                    continue
                
                self.logger.log(f"✅ {variant.gene} {variant.mutation}: Valid")
                valid_variants.append(variant)
                
            except Exception as e:
                self.logger.log(f"❌ {variant.gene}: Validation error - {e}", "ERROR")
        
        self.logger.log(f"Validation complete: {len(valid_variants)}/{len(variants)} variants valid")
        return valid_variants

    def process_variants(self, variants: List[VariantInfo], skip_json: bool = False) -> pd.DataFrame:
        """Process all variants and generate outputs"""
        self.logger.log("=== PROCESSING VARIANTS ===")
        
        # Group variants by gene
        gene_groups = {}
        for variant in variants:
            gene_groups.setdefault(variant.gene, []).append(variant)
        
        summary_data = []
        
        # Process each gene
        for gene_name, gene_variants in gene_groups.items():
            self.logger.log(f"Processing gene: {gene_name}")
            self.logger.update_stats("processed")
            
            try:
                # Find gene info
                gene_info = self.find_gene(gene_name)
                if not gene_info:
                    self.logger.log(f"Gene {gene_name} not found", "ERROR")
                    self.logger.update_stats("failed")
                    continue
                
                # Extract wild-type sequence
                wt_sequence = self.extract_protein_sequence(gene_info)
                
                # Create gene-specific directory
                gene_dir = self.output_dir / "sequences" / gene_name
                gene_dir.mkdir(exist_ok=True)
                
                # Write wild-type files
                wt_fasta = gene_dir / f"{gene_name}_WT.fasta"
                self.write_fasta(wt_fasta, f"{gene_name}_WT", wt_sequence)
                
                wt_config = None
                if not skip_json:
                    wt_config = self.output_dir / "configs" / f"{gene_name}_WT_config.json"
                    self.write_alphafold_config(wt_config, f"{gene_name}_WT", wt_sequence)
                
                # Add WT to summary
                summary_data.append({
                    "Gene": gene_name,
                    "Variant_ID": f"{gene_name}_WT",
                    "Mutation": "WT",
                    "Condition": "-",
                    "Priority": "-",
                    "Source": "-",
                    "Sequence_Length": len(wt_sequence),
                    "FASTA_File": str(wt_fasta.relative_to(self.output_dir)),
                    "AF3_Config": str(wt_config.relative_to(self.output_dir)) if wt_config else "",
                    "Status": "Success"
                })
                
                # Process each variant
                for variant in gene_variants:
                    try:
                        mutant_sequence = self.apply_mutation(wt_sequence, variant)
                        if not mutant_sequence:
                            continue
                        
                        variant_id = variant.mutation_id
                        
                        # Write mutant files
                        mut_fasta = gene_dir / f"{variant_id}.fasta"
                        self.write_fasta(mut_fasta, variant_id, mutant_sequence)
                        
                        mut_config = None
                        if not skip_json:
                            mut_config = self.output_dir / "configs" / f"{variant_id}_config.json"
                            self.write_alphafold_config(mut_config, variant_id, mutant_sequence)
                        
                        # Add to summary
                        summary_data.append({
                            "Gene": gene_name,
                            "Variant_ID": variant_id,
                            "Mutation": variant.hgvs_notation,
                            "Condition": variant.condition,
                            "Priority": variant.priority,
                            "Source": variant.source,
                            "Sequence_Length": len(mutant_sequence),
                            "FASTA_File": str(mut_fasta.relative_to(self.output_dir)),
                            "AF3_Config": str(mut_config.relative_to(self.output_dir)) if mut_config else "",
                            "Status": "Success"
                        })
                        
                    except Exception as e:
                        self.logger.log(f"Failed to process variant {variant.mutation}: {e}", "ERROR")
                        summary_data.append({
                            "Gene": gene_name,
                            "Variant_ID": variant.mutation_id,
                            "Mutation": variant.hgvs_notation,
                            "Condition": variant.condition,
                            "Priority": variant.priority,
                            "Source": variant.source,
                            "Sequence_Length": 0,
                            "FASTA_File": "",
                            "AF3_Config": "",
                            "Status": f"Error: {e}"
                        })
                
                self.logger.update_stats("successful")
                
            except Exception as e:
                self.logger.log(f"Failed to process gene {gene_name}: {e}", "ERROR")
                self.logger.update_stats("failed")
        
        # Create summary DataFrame
        summary_df = pd.DataFrame(summary_data)
        summary_path = self.output_dir / "processing_summary.tsv"
        summary_df.to_csv(summary_path, sep='\t', index=False)
        
        self.logger.log(f"Processing complete. Summary saved to: {summary_path}")
        return summary_df

    def generate_report(self, summary_df: pd.DataFrame) -> None:
        """Generate comprehensive analysis report"""
        report_path = self.output_dir / "reports" / "analysis_report.md"
        
        # Calculate statistics
        total_variants = len(summary_df)
        successful_variants = len(summary_df[summary_df["Status"] == "Success"])
        unique_genes = summary_df["Gene"].nunique()
        
        conditions = summary_df["Condition"].value_counts()
        priorities = summary_df["Priority"].value_counts()
        
        # Generate report
        report_content = f"""# Universal AlphaFold 3 Sequence Extraction Report

**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  
**Organism:** {self.organism}  
**Genome:** {self.fasta_file.name}  
**Annotations:** {self.gff_file.name}  

## Summary Statistics

- **Total variants processed:** {total_variants}
- **Successful extractions:** {successful_variants} ({successful_variants/total_variants*100:.1f}%)
- **Unique genes:** {unique_genes}
- **Output directory:** {self.output_dir}

## Variant Distribution

### By Condition
{conditions.to_string()}

### By Priority
{priorities.to_string()}

## File Structure

```
{self.output_dir.name}/
├── sequences/           # FASTA files organized by gene
├── configs/            # AlphaFold 3 configuration files
├── reports/            # Analysis reports and logs
├── processing_summary.tsv  # Detailed processing results
└── extraction.log      # Complete processing log
```

## Next Steps

1. **Review the processing log** for any warnings or errors
2. **Submit AF3 jobs** using the configuration files in `configs/`
3. **Monitor AF3 results** and organize them for downstream analysis
4. **Validate predictions** using experimental data when available

## Processing Log Summary

{self.logger.generate_summary()}

---
*Report generated by Universal AlphaFold 3 Sequence Extractor*
"""

        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        self.logger.log(f"Analysis report generated: {report_path}")

def build_argument_parser() -> argparse.ArgumentParser:
    """Build comprehensive argument parser"""
    parser = argparse.ArgumentParser(
        description="Universal AlphaFold 3 sequence extractor for structural variant analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:

  # Basic usage with TSV file
  python extract_af3_sequences_universal.py \\
    --gff annotations.gff3 \\
    --fasta genome.fasta \\
    --variants variants.tsv \\
    --output af3_analysis/

  # Single variant from command line
  python extract_af3_sequences_universal.py \\
    --gff annotations.gff3 \\
    --fasta genome.fasta \\
    --add-variant "BRCA1:p.Arg1450His:tumor:high" \\
    --output af3_analysis/

  # Validation only (no file generation)
  python extract_af3_sequences_universal.py \\
    --gff annotations.gff3 \\
    --fasta genome.fasta \\
    --variants variants.tsv \\
    --validate-only

  # Multiple input formats
  python extract_af3_sequences_universal.py \\
    --gff annotations.gff3 \\
    --fasta genome.fasta \\
    --variants variants.json \\
    --add-variant "TP53:p.Arg273His:control:high" \\
    --organism "Homo sapiens" \\
    --output af3_analysis/

INPUT FORMATS:

  TSV/CSV files should contain columns:
  - Gene (required): Gene identifier
  - Mutation (required): Mutation in format p.Ala123Val or A123V
  - Condition (optional): Experimental condition or sample type
  - Priority (optional): Variant priority (high/medium/low)
  - Source (optional): Data source or method

  Command-line variants format:
  GENE:MUTATION[:condition[:priority]]
  
  Examples:
  - "BRCA1:p.Arg1450His"
  - "TP53:A273V:tumor:high"
  - "EGFR:p.Leu858Arg:drug_resistant:medium"
        """
    )
    
    # Required arguments
    parser.add_argument("--gff", required=True, type=Path,
                       help="GFF3/GTF annotation file")
    parser.add_argument("--fasta", required=True, type=Path,
                       help="Genome FASTA file")
    
    # Input variants
    variant_group = parser.add_mutually_exclusive_group()
    variant_group.add_argument("--variants", type=Path,
                              help="Variants file (TSV, CSV, or JSON)")
    variant_group.add_argument("--add-variant", action="append",
                              help="Add single variant (format: GENE:MUTATION[:condition[:priority]])")
    
    # Output options
    parser.add_argument("--output", type=Path,
                       help="Output directory (required unless --validate-only)")
    parser.add_argument("--organism", default="Unknown",
                       help="Organism name for documentation")
    
    # Processing options
    parser.add_argument("--validate-only", action="store_true",
                       help="Only validate variants without generating files")
    parser.add_argument("--skip-json", action="store_true",
                       help="Skip AlphaFold 3 JSON configuration generation")
    parser.add_argument("--priority-filter", choices=["high", "medium", "low"],
                       help="Only process variants with specified priority")
    parser.add_argument("--condition-filter", 
                       help="Only process variants from specified condition")
    
    # Advanced options
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       default="INFO", help="Logging level")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing output directory")
    
    return parser

def main():
    """Main execution function"""
    parser = build_argument_parser()
    args = parser.parse_args()
    
    # Validation
    if not args.validate_only and not args.output:
        parser.error("--output is required unless using --validate-only")
    
    if not args.variants and not args.add_variant:
        parser.error("Either --variants file or --add-variant must be specified")
    
    # Check for output directory conflicts
    if args.output and args.output.exists() and not args.force:
        response = input(f"Output directory {args.output} exists. Overwrite? (y/N): ")
        if response.lower() != 'y':
            print("Operation cancelled.")
            sys.exit(1)
    
    try:
        # Initialize extractor
        output_dir = args.output or Path("temp_validation")
        extractor = UniversalSequenceExtractor(
            fasta_file=args.fasta,
            gff_file=args.gff,
            output_dir=output_dir,
            organism=args.organism,
            log_level=args.log_level
        )
        
        # Validate inputs and load data
        extractor.validate_inputs()
        extractor.load_genome()
        extractor.parse_annotations()
        
        # Load variants
        variants = []
        
        if args.variants:
            variants.extend(extractor.load_variants_from_file(args.variants))
        
        if args.add_variant:
            for variant_str in args.add_variant:
                variants.append(extractor.parse_command_line_variant(variant_str))
        
        # Apply filters
        if args.priority_filter:
            variants = [v for v in variants if v.priority == args.priority_filter]
            extractor.logger.log(f"Priority filter applied: {len(variants)} variants remain")
        
        if args.condition_filter:
            variants = [v for v in variants if v.condition == args.condition_filter]
            extractor.logger.log(f"Condition filter applied: {len(variants)} variants remain")
        
        if not variants:
            extractor.logger.log("No variants to process after filtering", "ERROR")
            sys.exit(1)
        
        # Validate variants
        valid_variants = extractor.validate_variants(variants)
        
        if args.validate_only:
            # Validation mode only
            if len(valid_variants) == len(variants):
                print(f"\n✅ All {len(variants)} variants are valid!")
                sys.exit(0)
            else:
                print(f"\n❌ {len(variants) - len(valid_variants)}/{len(variants)} variants failed validation")
                sys.exit(1)
        else:
            # Full processing mode
            if not valid_variants:
                extractor.logger.log("No valid variants to process", "ERROR")
                sys.exit(1)
            
            # Process variants and generate outputs
            summary_df = extractor.process_variants(valid_variants, skip_json=args.skip_json)
            extractor.generate_report(summary_df)
            
            # Final summary
            stats = extractor.logger.generate_summary()
            print(f"\n🎉 Processing complete!")
            print(f"📊 Processed: {stats['statistics']['processed']} genes")
            print(f"✅ Successful: {stats['statistics']['successful']} genes")
            print(f"❌ Failed: {stats['statistics']['failed']} genes")
            print(f"📁 Output directory: {output_dir}")
            
            if stats['errors'] > 0:
                print(f"⚠️  {stats['errors']} errors occurred (see log for details)")
                sys.exit(1)
            
    except Exception as e:
        print(f"❌ Fatal error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
