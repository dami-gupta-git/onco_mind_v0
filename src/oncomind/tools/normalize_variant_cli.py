#!/usr/bin/env python3
"""Command-line utility for variant normalization.

This tool normalizes variant notations across different formats:
- One-letter amino acid codes (V600E)
- Three-letter amino acid codes (Val600Glu)
- HGVS protein notation (p.V600E, p.Val600Glu)
- Structural variants (fusion, amplification)
- Indels (del, ins, dup, fs)

Usage:
    python -m oncomind.tools.normalize_variant_cli BRAF V600E
    python -m oncomind.tools.normalize_variant_cli BRAF Val600Glu
    python -m oncomind.tools.normalize_variant_cli BRAF p.V600E
    python -m oncomind.tools.normalize_variant_cli ALK fusion
    python -m oncomind.tools.normalize_variant_cli --batch variants.txt
    echo "BRAF,V600E" | python -m oncomind.tools.normalize_variant_cli --stdin

    # With genomic lookup (requires network)
    python -m oncomind.tools.normalize_variant_cli BRAF V600E --lookup

Output formats:
    --format json    JSON output (default)
    --format table   Human-readable table
    --format tsv     Tab-separated values
"""

import argparse
import asyncio
import json
import sys
from typing import TextIO

# Core normalization logic is in utils - CLI just provides interface
from oncomind.utils.variant_normalization import (
    normalize_variant_extended,
    normalize_variant_with_lookup,
)


def format_json(results: list[dict], pretty: bool = True) -> str:
    """Format results as JSON."""
    if len(results) == 1:
        return json.dumps(results[0], indent=2 if pretty else None)
    return json.dumps(results, indent=2 if pretty else None)


def format_table(results: list[dict]) -> str:
    """Format results as a human-readable table."""
    lines = []

    for result in results:
        lines.append("=" * 60)
        lines.append(f"Gene:             {result['gene']}")
        lines.append(f"Original:         {result['variant_original']}")
        lines.append(f"Normalized:       {result['variant_normalized']}")
        lines.append(f"Type:             {result['variant_type']}")
        lines.append(f"Chromosome:       {result.get('chromosome', 'N/A')}")
        lines.append(f"HGVS Protein:     {result.get('hgvs_protein', 'N/A')}")
        lines.append(f"HGVS Genomic:     {result.get('hgvs_genomic', 'N/A')}")
        lines.append(f"Position:         {result.get('position', 'N/A')}")
        lines.append(f"Allowed Type:     {result.get('is_allowed_type', 'N/A')}")

        # Genomic details (if lookup was performed)
        if result.get('gene_name') or result.get('gene_id'):
            lines.append("-" * 40)
            lines.append("Genomic Details:")
            if result.get('gene_name'):
                lines.append(f"  Gene Name:      {result['gene_name']}")
            if result.get('gene_id'):
                lines.append(f"  Gene ID:        {result['gene_id']}")
            if result.get('transcript_id'):
                lines.append(f"  Transcript:     {result['transcript_id']}")
            if result.get('exon'):
                lines.append(f"  Exon:           {result['exon']}")
            if result.get('genomic_position'):
                lines.append(f"  Genomic Pos:    {result['genomic_position']}")
            if result.get('ref_allele') and result.get('alt_allele'):
                lines.append(f"  Alleles:        {result['ref_allele']}>{result['alt_allele']}")

        if result.get('protein_change'):
            pc = result['protein_change']
            lines.append("-" * 40)
            lines.append("Protein Change:")
            lines.append(f"  Ref AA:         {pc.get('ref_aa', 'N/A')}")
            lines.append(f"  Alt AA:         {pc.get('alt_aa', 'N/A')}")
            lines.append(f"  Long Form:      {pc.get('long_form', 'N/A')}")

        if result.get('query_formats'):
            qf = result['query_formats']
            lines.append("-" * 40)
            lines.append("Query Formats:")
            lines.append(f"  MyVariant:      {qf.get('myvariant', 'N/A')}")
            lines.append(f"  VICC:           {qf.get('vicc', 'N/A')}")
            lines.append(f"  CIViC:          {qf.get('civic', 'N/A')}")

    lines.append("=" * 60)
    return "\n".join(lines)


def format_tsv(results: list[dict]) -> str:
    """Format results as TSV."""
    headers = [
        "gene", "chromosome", "variant_original", "variant_normalized", "variant_type",
        "hgvs_protein", "hgvs_genomic", "position", "genomic_position",
        "is_allowed_type", "ref_aa", "alt_aa", "gene_name", "gene_id", "transcript_id", "exon"
    ]

    lines = ["\t".join(headers)]

    for result in results:
        pc = result.get('protein_change') or {}
        row = [
            result['gene'],
            result.get('chromosome') or '',
            result['variant_original'],
            result['variant_normalized'],
            result['variant_type'],
            result.get('hgvs_protein') or '',
            result.get('hgvs_genomic') or '',
            str(result.get('position') or ''),
            str(result.get('genomic_position') or ''),
            str(result.get('is_allowed_type', '')),
            pc.get('ref_aa') or '',
            pc.get('alt_aa') or '',
            result.get('gene_name') or '',
            result.get('gene_id') or '',
            result.get('transcript_id') or '',
            result.get('exon') or '',
        ]
        lines.append("\t".join(row))

    return "\n".join(lines)


def parse_batch_line(line: str) -> tuple[str, str] | None:
    """Parse a line from batch input (gene,variant or gene\tvariant)."""
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    # Try comma separator
    if ',' in line:
        parts = line.split(',', 1)
        if len(parts) == 2:
            return parts[0].strip(), parts[1].strip()

    # Try tab separator
    if '\t' in line:
        parts = line.split('\t', 1)
        if len(parts) == 2:
            return parts[0].strip(), parts[1].strip()

    # Try space separator (for "BRAF V600E" format)
    parts = line.split(None, 1)
    if len(parts) == 2:
        return parts[0].strip(), parts[1].strip()

    return None


def process_batch(input_file: TextIO, lookup: bool = False) -> list[dict]:
    """Process batch input from a file or stdin."""
    results = []

    for line_num, line in enumerate(input_file, 1):
        parsed = parse_batch_line(line)
        if parsed:
            gene, variant = parsed
            try:
                if lookup:
                    result = asyncio.run(normalize_variant_with_lookup(gene, variant))
                else:
                    result = normalize_variant_extended(gene, variant)
                result['line_number'] = line_num
                results.append(result)
            except Exception as e:
                results.append({
                    'line_number': line_num,
                    'gene': gene,
                    'variant_original': variant,
                    'error': str(e)
                })

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Normalize variant notations to standard formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s BRAF V600E              # Normalize single variant
  %(prog)s BRAF Val600Glu          # Three-letter to one-letter
  %(prog)s BRAF p.V600E            # HGVS protein notation
  %(prog)s ALK fusion              # Structural variant
  %(prog)s --batch variants.txt    # Process batch file
  %(prog)s --format table EGFR L858R  # Table output
  %(prog)s --lookup BRAF V600E     # With genomic lookup (requires network)

Batch file format (one variant per line):
  BRAF,V600E
  EGFR,L858R
  ALK,fusion
  # Comments start with #
"""
    )

    parser.add_argument(
        'gene',
        nargs='?',
        help='Gene symbol (e.g., BRAF, EGFR)'
    )
    parser.add_argument(
        'variant',
        nargs='?',
        help='Variant notation (e.g., V600E, Val600Glu, p.V600E, fusion)'
    )
    parser.add_argument(
        '--batch', '-b',
        type=argparse.FileType('r'),
        help='Batch input file (one gene,variant per line)'
    )
    parser.add_argument(
        '--stdin', '-i',
        action='store_true',
        help='Read from stdin (batch mode)'
    )
    parser.add_argument(
        '--lookup', '-l',
        action='store_true',
        help='Look up genomic info from MyVariant.info (requires network)'
    )
    parser.add_argument(
        '--format', '-f',
        choices=['json', 'table', 'tsv'],
        default='json',
        help='Output format (default: json)'
    )
    parser.add_argument(
        '--compact', '-c',
        action='store_true',
        help='Compact JSON output (no indentation)'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Only output the normalized variant (no metadata)'
    )

    args = parser.parse_args()

    # Determine input mode
    if args.batch:
        results = process_batch(args.batch, lookup=args.lookup)
    elif args.stdin:
        results = process_batch(sys.stdin, lookup=args.lookup)
    elif args.gene and args.variant:
        if args.lookup:
            results = [asyncio.run(normalize_variant_with_lookup(args.gene, args.variant))]
        else:
            results = [normalize_variant_extended(args.gene, args.variant)]
    else:
        parser.print_help()
        sys.exit(1)

    # Handle quiet mode
    if args.quiet:
        for result in results:
            if 'error' in result:
                print(f"ERROR: {result.get('gene', '?')} {result.get('variant_original', '?')}: {result['error']}", file=sys.stderr)
            else:
                chrom = result.get('chromosome', '')
                genomic = result.get('hgvs_genomic', '')
                print(f"{result['gene']}\t{result['variant_normalized']}\t{chrom}\t{genomic}")
        return

    # Format output
    if args.format == 'json':
        print(format_json(results, pretty=not args.compact))
    elif args.format == 'table':
        print(format_table(results))
    elif args.format == 'tsv':
        print(format_tsv(results))


if __name__ == '__main__':
    main()
