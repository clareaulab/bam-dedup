# bam-dedup

A fast, **JVM-free** toolkit for removing duplicate reads from BAM/CRAM. It
handles the **two major kinds of duplication event** in sequencing data, each a
faithful, independent port of the standard reference tool, with the
performance-critical inner loops accelerated in Cython. All BAM/CRAM I/O goes
through [pysam](https://github.com/pysam-developers/pysam) (htslib) — no Java
required.

| Duplication kind | Module | Reproduces | What it does |
|---|---|---|---|
| **PCR / optical duplicates** | `dedup.picardlike` | [Picard MarkDuplicates](https://broadinstitute.github.io/picard/) | flags or removes duplicate *records* by 5′ position + orientation |
| **Molecular / UMI duplicates** | `dedup.fgbiolike` | [fgbio](https://github.com/fulcrumgenomics/fgbio) GroupReadsByUmi + CallMolecularConsensusReads | collapses reads sharing a UMI into one error-corrected *consensus* read |

Both are validated against their reference tool:

* `picardlike` is **bit-identical to Picard MarkDuplicates 2.18.21** for the
  duplicate flag *and* the optical/sequencing-duplicate classification.
* `fgbiolike` is **bit-identical to fgbio 4.1.0** (`GroupReadsByUmi -s Identity`
  and `CallMolecularConsensusReads`) for the grouped `MI` assignment and every
  consensus base, quality, and `cD`/`cM`/`cE`/`cd`/`ce` tag — verified
  record-for-record on both a bundled test BAM and a 176k-read MAESTER file.

Both run faster than the corresponding Java tool.

> The importable modules `dedup.picardlike` / `dedup.fgbiolike` are independent
> reimplementations, **not** affiliated with or derived from the Picard or fgbio
> source code; those names are referenced only to describe the behavior
> reproduced.

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

pysam-like: point at an input BAM, call one function, get an output BAM.

### PCR / optical duplicates (`picardlike`)

Input must be **coordinate-sorted** (e.g. `samtools sort`).

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

### Molecular / UMI duplicates (`fgbiolike`)

Reads must carry a UMI tag (default `UB`); a cell-barcode tag (default `CB`) is
used when present. The one-shot `consensus()` runs the two fgbio steps —
grouping reads into molecules, then collapsing each molecule into a consensus
read:

```python
from dedup import fgbiolike

# group by UMI (+ cell barcode) and emit error-corrected consensus reads
fgbiolike.consensus("input.bam", "consensus.bam", umi_tag="UB", min_reads=3)
```

Or run the two stages explicitly:

```python
fgbiolike.group_reads_by_umi("input.bam", "grouped.bam", raw_tag="UB")
fgbiolike.call_molecular_consensus_reads("grouped.bam", "consensus.bam",
                                         tag="MI", min_reads=3)
```

`consensus(...)` / `call_molecular_consensus_reads(...)` options:

| option | default | meaning |
|---|---|---|
| `min_reads` | `3` | minimum reads in a molecule to emit a consensus |
| `umi_tag` | `"UB"` | raw UMI tag used for grouping |
| `cell_tag` | `"CB"` | cell-barcode tag; consensus is called per cell |
| `consensus_tag` | `"MI"` | tag grouped on in the consensus step (see note) |
| `on_cell_collision` | `"split"` | if a molecule spans >1 cell: `split`, `merge`, or `error` |

> **Grouping key.** The consensus step groups by the assigned molecular id `MI`
> (per position × cell × UMI). This is the correct, crash-free equivalent of
> fgbio's `-t UB`: fgbio *aborts* if a raw-UMI group ever spans two cell
> barcodes, whereas `on_cell_collision="split"` (the default) simply calls one
> consensus per cell.

### Command line

Installing the package provides two console commands:

```bash
# PCR/optical duplicate marking
bam-dedup -i input.sorted.bam -o dedup.bam --remove-duplicates

# UMI molecular-consensus deduplication (group + call in one step)
bam-consensus consensus -i input.bam -o consensus.bam -M 3
# or the individual stages:
bam-consensus group -i input.bam -o grouped.bam -t UB
bam-consensus call  -i grouped.bam -o consensus.bam -t MI -M 3
```

---

## Scope

**`picardlike`** — faithful to Picard for coordinate-sorted input, the default
`SUM_OF_BASE_QUALITIES` scoring strategy, standard Illumina read names, and
single/multiple libraries. *Not yet handled:* queryname-sorted input, flow-based
mode, and the non-default scoring strategies (`TOTAL_MAPPED_REFERENCE_LENGTH`,
`RANDOM`).

**`fgbiolike`** — faithful to fgbio for the `Identity` grouping strategy
(`edits=0`) and the single-end / fragment consensus path (the MAESTER-style
input), with fgbio's default error model (pre=45, post=40),
`min-input-base-quality=10`, and per-base tags. *Not yet handled:* the
Adjacency/Edit/Paired strategies (`edits > 0`) and paired-end / duplex /
overlapping-base consensus.

---

## Tests

```bash
pip install -e ".[test]"
pytest
```

`tests/test_dedup.py` validates duplicate and optical-duplicate counts against
Picard's reference output on the bundled BAMs; `tests/test_fgbiolike.py`
validates UMI grouping and consensus reads against **fgbio golden outputs**
(`tests/data/umi_consensus.*.golden.bam`). Both suites also check that the Cython
and pure-Python paths agree.

## Benchmarking against Picard

Reproducible head-to-head scripts (and a bundled `picard.jar`) live in
[`tests/benchmark/`](tests/benchmark/). See
[`tests/benchmark/README.md`](tests/benchmark/README.md) for the concordance and
timing harness.

## License

MIT © Caleb Lareau
