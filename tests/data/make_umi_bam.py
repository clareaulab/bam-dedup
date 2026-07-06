#!/usr/bin/env python3
"""
Generate a small deterministic single-end UMI BAM for testing dedup.fgbiolike.

Reads carry UB (UMI) and CB (cell barcode) tags, mimicking MAESTER-style data:
several molecules (cell x UMI x position), each with a few reads and a sprinkling
of sequencing errors, plus a sub-threshold singleton and a reverse-strand family.

Run once to (re)create umi_consensus.bam; the fgbio golden outputs are produced
separately (see tests/data/README or the repo history).
"""
import os
import random

import pysam

REF = "chrT"
REF_LEN = 300
HERE = os.path.dirname(__file__)

# a fixed 120bp reference-ish template the reads are drawn from
random.seed(7)
TEMPLATE = "".join(random.choice("ACGT") for _ in range(120))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def make_reads(header, mol_id_start=0):
    """Return a list of AlignedSegment for the whole file."""
    reads = []
    counter = [0]

    def emit(cell, umi, pos, strand, n, length=80, err_positions=()):
        """n reads for one molecule; err_positions: list of (read_idx, base_idx, newbase)."""
        base_seq = TEMPLATE[:length]
        for k in range(n):
            seq = list(base_seq)
            for (ri, bi, nb) in err_positions:
                if ri == k:
                    seq[bi] = nb
            seq = "".join(seq)
            quals = [37] * length
            # give the error bases a decent quality so they count
            a = pysam.AlignedSegment(header)
            a.query_name = "READ:{}:{}:{}".format(cell, umi, counter[0])
            counter[0] += 1
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigarstring = "{}M".format(length)
            if strand == "-":
                a.is_reverse = True
                a.query_sequence = revcomp(seq)
                a.query_qualities = pysam.qualitystring_to_array(
                    "".join(chr(q + 33) for q in reversed(quals)))
            else:
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array(
                    "".join(chr(q + 33) for q in quals))
            a.set_tag("UB", umi, value_type="Z")
            a.set_tag("CB", cell, value_type="Z")
            reads.append(a)

    # molecule B: cell C0, umi GGGGTTTT, +strand, 4 reads, one error
    emit("C0-1", "GGGGTTTT", 20, "+", 4, err_positions=[(1, 5, "A")])
    # molecule A: cell C1, umi AAAACCCC, +strand, 5 reads, two errors
    emit("C1-1", "AAAACCCC", 20, "+", 5, err_positions=[(0, 10, "G"), (3, 40, "T")])
    # molecule C: cell C2, SAME umi AAAACCCC at same pos. C1 sorts immediately
    # before C2 and AAAACCCC is the only UMI for both, so the two families are
    # *adjacent* when grouped by UB -> a cross-cell collision (fgbio would abort
    # CallMolecularConsensusReads -t UB). Grouping by MI keeps them separate.
    emit("C2-1", "AAAACCCC", 20, "+", 3, err_positions=[(2, 15, "C")])
    # molecule D: reverse strand family, cell C1, umi TTTTGGGG, 4 reads
    emit("C1-1", "TTTTGGGG", 50, "-", 4, err_positions=[(0, 7, "A")])
    # singleton-ish: only 2 reads -> below min-reads=3, no consensus
    emit("C4-1", "CCCCAAAA", 20, "+", 2)
    return reads


def main():
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    reads = make_reads(header)
    # coordinate sort (stable by pos then name to be deterministic)
    reads.sort(key=lambda r: (r.reference_start, r.query_name))
    out = os.path.join(HERE, "umi_consensus.bam")
    with pysam.AlignmentFile(out, "wb", header=header) as bam:
        for r in reads:
            bam.write(r)
    print("wrote", out, "with", len(reads), "reads")


if __name__ == "__main__":
    main()
