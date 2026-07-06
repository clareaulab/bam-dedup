#!/usr/bin/env python3
"""
Build a stress-test BAM by amplifying real read pairs into duplicate sets.

For each real pair we emit COPIES copies at the SAME genomic position/CIGAR/seq
(so they are true duplicates), giving each copy:
  * a synthetic read name "SYNTH:1:<tile>:<x>:<y>" (5 fields -> parseable coords)
  * an optical sub-cluster (first 3 copies within a few px, same tile) plus
    PCR-only copies (hundreds of px away) and a far-tile copy
  * distinct, non-tied base-quality sums so the representative is unambiguous

We do NOT decide who is a duplicate -- both tools classify independently; the
point is that they must AGREE.
"""
import sys
import pysam

COPIES_SPEC = [
    # (tile, dx, dy)   relative to a per-pair base (x0=1000, y0=2000)
    (1101, 0, 0),      # optical cluster A
    (1101, 5, 3),      # within 100px of A -> optical dup
    (1101, 10, 8),     # within 100px of A -> optical dup
    (1101, 800, 0),    # same tile but >100px -> PCR dup, not optical
    (1101, 1600, 0),   # same tile, far        -> PCR dup, not optical
    (1202, 0, 0),      # different tile        -> PCR dup, never optical
]


def clone_read(r, name, qual_value):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = r.flag
    a.reference_id = r.reference_id
    a.reference_start = r.reference_start
    a.mapping_quality = r.mapping_quality
    a.cigartuples = r.cigartuples
    a.next_reference_id = r.next_reference_id
    a.next_reference_start = r.next_reference_start
    a.template_length = r.template_length
    seq = r.query_sequence
    a.query_sequence = seq
    if seq is not None:
        a.query_qualities = pysam.qualitystring_to_array(
            chr(33 + qual_value) * len(seq))
    # carry over read group so library grouping is exercised
    try:
        a.set_tag("RG", r.get_tag("RG"))
    except KeyError:
        pass
    return a


def main(in_bam, out_bam):
    with pysam.AlignmentFile(in_bam, "rb") as bam:
        header = bam.header.to_dict()
        # collect mate pairs (primary, mapped, both mates mapped)
        pending = {}
        pairs = []
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if not r.is_paired or r.mate_is_unmapped:
                continue
            k = r.query_name
            if k in pending:
                r1 = pending.pop(k)
                pairs.append((r1, r))
            else:
                pending[k] = r

    out_records = []
    pair_idx = 0
    for (r1, r2) in pairs:
        pair_idx += 1
        x0, y0 = 1000, 2000
        for ci, (tile, dx, dy) in enumerate(COPIES_SPEC):
            x = x0 + dx
            y = y0 + dy
            # unique per (pair, copy); 5 colon fields so tile/x/y parse out
            name = "SYNTH{}:1:{}:{}:{}".format(pair_idx, tile, x, y)
            # distinct score per copy: quality value 20 + ci  (no ties)
            qv = 20 + ci
            out_records.append(clone_read(r1, name, qv))
            out_records.append(clone_read(r2, name, qv))

    # sort by coordinate for a valid coordinate-sorted BAM
    out_records.sort(key=lambda a: (a.reference_id if a.reference_id >= 0 else 1 << 30,
                                    a.reference_start))
    header.setdefault("HD", {})
    header["HD"]["SO"] = "coordinate"
    with pysam.AlignmentFile(out_bam, "wb", header=header) as out:
        for a in out_records:
            out.write(a)
    print("wrote {} records ({} pairs x {} copies)".format(
        len(out_records), len(pairs), len(COPIES_SPEC)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
