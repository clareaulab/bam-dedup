"""
Tests for bam-dedup.

Expected duplicate counts were validated to be bit-identical to Picard
MarkDuplicates 2.18.21 (default options, coordinate-sorted) on these BAMs.
"""
import os

import pysam
import pytest

from dedup import picardlike

DATA = os.path.join(os.path.dirname(__file__), "data")

# (filename, total_reads, picard_duplicate_count, picard_optical_count)
CASES = [
    ("NA12891.bam", 751, 6, 0),
    ("NA12892.bam", 742, 2, 2),
    ("synthetic_optical.bam", 1824, 1526, 632),
]


def _count_dups(path):
    n = 0
    with pysam.AlignmentFile(path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_duplicate:
                n += 1
    return n


def _count_reads(path):
    with pysam.AlignmentFile(path, "rb") as bam:
        return sum(1 for _ in bam.fetch(until_eof=True))


@pytest.mark.parametrize("fname,total,exp_dup,exp_opt", CASES)
def test_mark_duplicates_counts(tmp_path, fname, total, exp_dup, exp_opt):
    """Duplicate flag counts match Picard's reference output exactly."""
    out = str(tmp_path / "marked.bam")
    picardlike.mark_duplicates(os.path.join(DATA, fname), out)
    assert _count_reads(out) == total          # no records dropped when flagging
    assert _count_dups(out) == exp_dup


@pytest.mark.parametrize("fname,total,exp_dup,exp_opt", CASES)
def test_deduplicate_removes(tmp_path, fname, total, exp_dup, exp_opt):
    """deduplicate() drops exactly the duplicate records."""
    out = str(tmp_path / "dedup.bam")
    picardlike.deduplicate(os.path.join(DATA, fname), out)
    assert _count_reads(out) == total - exp_dup
    assert _count_dups(out) == 0               # nothing flagged remains


@pytest.mark.parametrize("fname,total,exp_dup,exp_opt", CASES)
def test_optical_counts(tmp_path, fname, total, exp_dup, exp_opt):
    """Optical/sequencing-duplicate detection matches Picard's DT:Z:SQ set."""
    src = os.path.join(DATA, fname)
    _, opt_idx = picardlike.mark_duplicates(
        src, str(tmp_path / "m.bam"), metrics_file="x")
    assert len(opt_idx) == exp_opt


def test_pure_python_matches_cython(tmp_path):
    """The Cython path and the pure-Python fallback produce identical flags."""
    src = os.path.join(DATA, "synthetic_optical.bam")

    saved = picardlike._HAVE_FAST
    try:
        picardlike._HAVE_FAST = True
        dup_fast, opt_fast = picardlike.mark_duplicates(
            src, str(tmp_path / "fast.bam"), metrics_file="x")
        picardlike._HAVE_FAST = False
        dup_pure, opt_pure = picardlike.mark_duplicates(
            src, str(tmp_path / "pure.bam"), metrics_file="x")
    finally:
        picardlike._HAVE_FAST = saved

    assert dup_fast == dup_pure
    assert opt_fast == opt_pure


def test_cython_extension_built():
    """The compiled acceleration module should be importable after install."""
    from dedup import HAVE_CYTHON
    assert HAVE_CYTHON, "dedup._fast Cython extension was not built"
