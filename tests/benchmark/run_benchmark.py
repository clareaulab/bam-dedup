#!/usr/bin/env python3
"""
Head-to-head benchmark: Picard MarkDuplicates (Java) vs bam-dedup (this package).

For a given coordinate-sorted BAM it:
  1. runs Picard MarkDuplicates (using the bundled picard.jar),
  2. runs bam-dedup (dedup.picardlike.mark_duplicates),
  3. times both (wall clock), and
  4. compares the duplicate flags record-by-record for concordance.

Examples
--------
# benchmark on a bundled test BAM
python run_benchmark.py --input ../data/synthetic_optical.bam

# benchmark on your own BAM
python run_benchmark.py --input /path/to/coord_sorted.bam

# amplify a small BAM into a big duplicate-heavy BAM first, then benchmark
python run_benchmark.py --input ../data/NA12891.bam --amplify 700

Requires: Java on PATH, pysam, and this package installed (pip install -e ..).
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_JAR = os.path.join(HERE, "picard.jar")


def _count(path, dup_only=False):
    n = 0
    with pysam.AlignmentFile(path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if not dup_only or r.is_duplicate:
                n += 1
    return n


def run_picard(jar, in_bam, out_bam, metrics):
    t0 = time.time()
    subprocess.run(
        ["java", "-jar", jar, "MarkDuplicates",
         "I=" + in_bam, "O=" + out_bam, "M=" + metrics,
         "VALIDATION_STRINGENCY=LENIENT"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.time() - t0


def run_bamdedup(in_bam, out_bam):
    from dedup import picardlike
    t0 = time.time()
    picardlike.mark_duplicates(in_bam, out_bam)
    return time.time() - t0


def amplify(in_bam, out_bam, copies):
    """Amplify each pair into `copies` duplicate copies with jittered names."""
    import random
    random.seed(0)
    with pysam.AlignmentFile(in_bam, "rb") as b:
        header = b.header.to_dict()
        pending, pairs = {}, []
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if not r.is_paired or r.mate_is_unmapped:
                continue
            if r.query_name in pending:
                pairs.append((pending.pop(r.query_name), r))
            else:
                pending[r.query_name] = r

    def clone(r, name, qv):
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
        s = r.query_sequence
        a.query_sequence = s
        if s is not None:
            a.query_qualities = pysam.qualitystring_to_array(chr(33 + qv) * len(s))
        try:
            a.set_tag("RG", r.get_tag("RG"))
        except KeyError:
            pass
        return a

    recs = []
    for pi, (r1, r2) in enumerate(pairs):
        for c in range(copies):
            tile = 1101 + (c % 4)
            x, y = 1000 + random.randint(0, 30000), 2000 + random.randint(0, 30000)
            name = "S{}c{}:1:{}:{}:{}".format(pi, c, tile, x, y)
            qv = 15 + (c % 20)
            recs.append(clone(r1, name, qv))
            recs.append(clone(r2, name, qv))
    recs.sort(key=lambda a: (a.reference_id if a.reference_id >= 0 else 1 << 30,
                             a.reference_start))
    header.setdefault("HD", {})
    header["HD"]["SO"] = "coordinate"
    with pysam.AlignmentFile(out_bam, "wb", header=header) as o:
        for a in recs:
            o.write(a)
    return out_bam


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="coordinate-sorted BAM")
    p.add_argument("--jar", default=DEFAULT_JAR,
                   help="path to picard.jar (default: bundled)")
    p.add_argument("--amplify", type=int, default=0,
                   help="amplify each pair into N copies before benchmarking")
    args = p.parse_args(argv)

    from compare_flags import compare  # local module

    tmp = tempfile.mkdtemp(prefix="bamdedup_bench_")
    in_bam = args.input
    if args.amplify > 0:
        in_bam = os.path.join(tmp, "amplified.bam")
        print("Amplifying {} x{} ...".format(args.input, args.amplify))
        amplify(args.input, in_bam, args.amplify)

    n = _count(in_bam)
    print("input: {} ({} reads)".format(in_bam, n))
    print("-" * 60)

    picard_out = os.path.join(tmp, "picard.bam")
    dedup_out = os.path.join(tmp, "bamdedup.bam")
    metrics = os.path.join(tmp, "picard.metrics")

    if not os.path.exists(args.jar):
        sys.exit("picard.jar not found at {}".format(args.jar))

    t_picard = run_picard(args.jar, in_bam, picard_out, metrics)
    t_dedup = run_bamdedup(in_bam, dedup_out)

    print("Picard    : {:7.2f} s   ({} dups)".format(
        t_picard, _count(picard_out, dup_only=True)))
    print("bam-dedup : {:7.2f} s   ({} dups)".format(
        t_dedup, _count(dedup_out, dup_only=True)))
    speedup = t_picard / t_dedup if t_dedup else float("nan")
    print("speedup   : {:.2f}x {}".format(
        speedup, "(bam-dedup faster)" if speedup > 1 else "(Picard faster)"))
    print("-" * 60)
    ok = compare(picard_out, dedup_out, "picard", "bam-dedup")
    print("-" * 60)
    print("CONCORDANT" if ok else "MISMATCH DETECTED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    sys.path.insert(0, HERE)
    main()
