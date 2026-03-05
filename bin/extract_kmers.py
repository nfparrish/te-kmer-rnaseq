#!/usr/bin/env python3
"""
extract_kmers.py – Decompose TE-mapped reads into k-mers (sliding window).

Reads a BAM file of reads mapped to TE consensus sequences, filters by
mismatch/clip thresholds (à la Kojima et al. 2021 hervk_kmers), and
emits per-k-mer counts.  Both mates of a read-pair are processed.

Output: gzipped TSV with columns  kmer  te_family  count
"""

import argparse
import gzip
import sys
from collections import Counter

import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True, help="Input BAM (mapped to TE consensus)")
    p.add_argument("--k", type=int, default=50, help="K-mer length [50]")
    p.add_argument("--max-mismatch", type=int, default=10,
                   help="Discard reads with > N mismatches [10]")
    p.add_argument("--max-clip", type=int, default=10,
                   help="Discard reads with > N soft-clipped bases [10]")
    p.add_argument("--output", required=True, help="Output TSV(.gz)")
    return p.parse_args()


def count_mismatches(read):
    """Return number of mismatches from NM tag, falling back to MD parsing."""
    try:
        return read.get_tag("NM")
    except KeyError:
        return 0


def total_soft_clip(read):
    """Sum of soft-clipped bases from CIGAR."""
    if read.cigartuples is None:
        return 0
    return sum(length for op, length in read.cigartuples if op == 4)


def extract_kmers_from_seq(seq, k):
    """Yield all k-mers of length k from a sequence string."""
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            yield kmer


def main():
    args = parse_args()
    k = args.k
    counts = Counter()        # (kmer, te_family) -> count
    n_reads = 0
    n_pass = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            n_reads += 1

            # ── Filter by mismatch / clip (Kojima-style) ───────────
            if count_mismatches(read) > args.max_mismatch:
                continue
            if total_soft_clip(read) > args.max_clip:
                continue

            n_pass += 1
            te_family = read.reference_name  # consensus name = TE family
            seq = read.query_sequence
            if seq is None:
                continue

            for kmer in extract_kmers_from_seq(seq, k):
                counts[(kmer, te_family)] += 1

    # ── Write output ────────────────────────────────────────────────
    opener = gzip.open if args.output.endswith(".gz") else open
    with opener(args.output, "wt") as fout:
        fout.write("kmer\tte_family\tcount\n")
        for (kmer, te_family), cnt in sorted(counts.items()):
            fout.write(f"{kmer}\t{te_family}\t{cnt}\n")

    print(f"[extract_kmers] reads={n_reads}  passed={n_pass}  "
          f"unique_kmers={len(counts)}", file=sys.stderr)


if __name__ == "__main__":
    main()
