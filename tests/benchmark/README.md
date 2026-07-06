# Benchmarking bam-dedup against Picard

This directory contains a reproducible head-to-head harness comparing
`bam-dedup` with the reference Java tool, Picard MarkDuplicates.

## Contents

| file | purpose |
|---|---|
| `picard.jar` | Picard 2.18.21 (the reference implementation) |
| `run_benchmark.py` | runs both tools, times them, checks concordance |
| `compare_flags.py` | record-by-record duplicate-flag comparison |
| `make_stress_bam.py` | builds a duplicate-heavy BAM from real read pairs |

> **Note:** `picard.jar` (~15 MB) is kept here for local reproducibility but is
> **excluded from the published sdist/wheel** (see `MANIFEST.in`) so PyPI
> downloads stay small. It ships only in the source repository.

## Requirements

- Java (for `picard.jar`) on your `PATH`
- `samtools` (optional, for inspecting outputs)
- `bam-dedup` installed: `pip install -e ..` from this directory's parent

## Quick start

```bash
# from tests/benchmark/
# 1) concordance + timing on a bundled test BAM
python run_benchmark.py --input ../data/synthetic_optical.bam

# 2) a realistic large-scale run: amplify a small BAM into a big one first
python run_benchmark.py --input ../data/NA12891.bam --amplify 700

# 3) benchmark on your own coordinate-sorted BAM
python run_benchmark.py --input /path/to/your.sorted.bam
```

Example output:

```
input: .../amplified.bam (212800 reads)
------------------------------------------------------------
Picard    :    4.53 s   (212502 dups)
bam-dedup :    2.98 s   (212502 dups)
speedup   : 1.52x (bam-dedup faster)
------------------------------------------------------------
records compared : 212800
picard           dups: 212502
bam-dedup        dups: 212502
AGREE            : 212800 (100.0000%)
DISAGREE         : 0
------------------------------------------------------------
CONCORDANT
```

## Building a custom stress BAM

`make_stress_bam.py` amplifies each real read pair into a controlled duplicate
set (an optical sub-cluster + PCR-only copies + a different-tile copy), which
exercises the representative-selection and optical-duplicate graph paths:

```bash
python make_stress_bam.py ../data/NA12891.bam stress.bam
samtools sort -o stress.sorted.bam stress.bam
python run_benchmark.py --input stress.sorted.bam
```

## What "concordant" means

`compare_flags.py` keys every record by
`(query_name, flag-without-dup-bit, ref, pos, read1, secondary, supplementary)`
and asserts the `0x400` duplicate flag matches between the two outputs. On the
bundled and amplified BAMs, `bam-dedup` is 100% concordant with Picard for both
the library-duplicate flag and the optical/sequencing classification.
