# bam-dedup

A fast, **JVM-free** reimplementation of [Picard
MarkDuplicates](https://broadinstitute.github.io/picard/). The core
duplicate-marking algorithm is a faithful, independent port of the
Picard/htsjdk logic (coordinate-sorted, Illumina paired-end, default options),
with the performance-critical inner loops accelerated in Cython. All BAM/CRAM
I/O goes through [pysam](https://github.com/pysam-developers/pysam) (htslib) —
no Java required.

On real data it is **bit-identical to Picard MarkDuplicates 2.18.21** for the
duplicate flag *and* the optical/sequencing-duplicate classification, while
running faster than the Java tool.

> The importable module is named `dedup.picardlike` — this project is an
> independent reimplementation and is **not** affiliated with or derived from
> the Picard source code; "Picard" is referenced only to describe the behavior
> it reproduces.

---

## Installation

### From PyPI (once published)

```bash
pip install bam-dedup
```

### From source

Requires a C compiler (clang/gcc) and Python ≥ 3.8. Cython is pulled in
automatically as a build dependency via `pyproject.toml`; the only *runtime*
dependency is `pysam`.

```bash
git clone https://github.com/caleblareau/bam-dedup.git
cd bam-dedup
pip install .
```

### For development (editable install + tests)

An editable install compiles the Cython extension (`dedup._fast`) in place, so
`pytest` can find it:

```bash
pip install -e ".[test]"
pytest
```

Verify the compiled acceleration is active:

```python
import dedup
print(dedup.HAVE_CYTHON)   # True when the Cython extension is built
```

If `HAVE_CYTHON` is `False`, the package still runs correctly using a pure-Python
fallback that produces bit-identical results (just slower).

---

## Usage

pysam-like: point at an input BAM, call one function, get an output BAM. Input
must be **coordinate-sorted** (e.g. `samtools sort`).

```python
from dedup import picardlike

# Picard-style: flag duplicates in place (records kept, 0x400 flag set)
picardlike.mark_duplicates("input.sorted.bam", "marked.bam")

# Produce a deduplicated BAM (duplicate records removed)
picardlike.deduplicate("input.sorted.bam", "dedup.bam")
```

`mark_duplicates(input_bam, output_bam, **options)`:

| option | default | meaning |
|---|---|---|
| `remove_duplicates` | `False` | drop duplicate records instead of only flagging them |
| `remove_sequencing_duplicates` | `False` | drop optical/sequencing duplicates |
| `read_name_regex_enabled` | `True` | parse tile/x/y from read names for optical detection |
| `optical_pixel_distance` | `100` | optical duplicate pixel distance |

`deduplicate(...)` is a convenience wrapper equal to
`mark_duplicates(..., remove_duplicates=True)`.

### Command line

Installing the package also provides a `bam-dedup` console command:

```bash
bam-dedup -i input.sorted.bam -o dedup.bam --remove-duplicates
```

```
bam-dedup -i INPUT -o OUTPUT
    --remove-duplicates              remove duplicates instead of flagging
    --remove-sequencing-duplicates   remove optical/sequencing duplicates
    --no-optical                     disable optical-duplicate detection
    --optical-pixel-distance N       optical pixel distance (default 100)
```

---

## Scope

Faithful to Picard for coordinate-sorted input, the default
`SUM_OF_BASE_QUALITIES` scoring strategy, standard Illumina read names, and
single/multiple libraries. **Not yet handled:** queryname-sorted input,
UMI/barcode-aware marking, flow-based mode, and the non-default scoring
strategies (`TOTAL_MAPPED_REFERENCE_LENGTH`, `RANDOM`).

---

## Tests

```bash
pip install -e ".[test]"
pytest
```

The suite (`tests/test_dedup.py`) validates duplicate and optical-duplicate
counts against Picard's reference output on the bundled BAMs
(`tests/data/`) and checks that the Cython and pure-Python paths agree.

## Benchmarking against Picard

Reproducible head-to-head scripts (and a bundled `picard.jar`) live in
[`tests/benchmark/`](tests/benchmark/). See
[`tests/benchmark/README.md`](tests/benchmark/README.md) for the concordance and
timing harness.

## License

MIT © Caleb Lareau
