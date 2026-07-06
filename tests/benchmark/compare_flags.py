#!/usr/bin/env python3
"""
Compare duplicate flags between two BAMs record-by-record.

Usage:
    python compare_flags.py a.bam b.bam

Returns exit code 0 if every record's duplicate flag agrees, 1 otherwise.
Used by the benchmark harness to confirm bam-dedup is concordant with Picard.
"""
import sys

import pysam


def load_dup_status(path):
    """Map a stable per-record key -> is_duplicate bool."""
    status = {}
    with pysam.AlignmentFile(path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            key = (r.query_name,
                   r.flag & ~0x400,          # flag with the duplicate bit cleared
                   r.reference_id,
                   r.reference_start,
                   r.is_read1,
                   r.is_secondary,
                   r.is_supplementary)
            status[key] = r.is_duplicate
    return status


def compare(bam_a, bam_b, label_a="A", label_b="B"):
    a = load_dup_status(bam_a)
    b = load_dup_status(bam_b)

    shared = set(a) & set(b)
    only_a = len(set(a) - set(b))
    only_b = len(set(b) - set(a))

    agree = disagree = 0
    dup_a = dup_b = 0
    for k in shared:
        if a[k]:
            dup_a += 1
        if b[k]:
            dup_b += 1
        if a[k] == b[k]:
            agree += 1
        else:
            disagree += 1

    total = agree + disagree
    print("records compared : {}".format(total))
    print("{:<16} dups: {}".format(label_a, dup_a))
    print("{:<16} dups: {}".format(label_b, dup_b))
    print("AGREE            : {} ({:.4f}%)".format(
        agree, 100.0 * agree / total if total else 0.0))
    print("DISAGREE         : {}".format(disagree))
    if only_a or only_b:
        print("records only in {} / only in {}: {} / {}".format(
            label_a, label_b, only_a, only_b))
    return disagree == 0 and only_a == 0 and only_b == 0


if __name__ == "__main__":
    ok = compare(sys.argv[1], sys.argv[2], "picard", "bam-dedup")
    sys.exit(0 if ok else 1)
