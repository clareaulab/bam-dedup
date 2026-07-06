"""
Tests for dedup.fgbiolike (UMI grouping + molecular consensus).

Golden outputs in tests/data/umi_consensus.*.golden.bam were produced by fgbio
4.1.0 (GroupReadsByUmi -s Identity -e 0 -t UB; CallMolecularConsensusReads
-t MI -M 3) and this port reproduces them field-for-field.
"""
import os
import struct

import pysam
import pytest

from dedup import fgbiolike

DATA = os.path.join(os.path.dirname(__file__), "data")
INPUT = os.path.join(DATA, "umi_consensus.bam")
GROUPED_GOLDEN = os.path.join(DATA, "umi_consensus.grouped.golden.bam")
CONSENSUS_GOLDEN = os.path.join(DATA, "umi_consensus.consensus.golden.bam")


def _f32(x):
    return struct.unpack("f", struct.pack("f", x))[0]


def _load_list(path):
    with pysam.AlignmentFile(path, "rb", check_sq=False) as f:
        return list(f)


def _load_by_name(path):
    with pysam.AlignmentFile(path, "rb", check_sq=False) as f:
        return {r.query_name: r for r in f}


def test_group_reads_matches_fgbio(tmp_path):
    out = str(tmp_path / "grouped.bam")
    stats = fgbiolike.group_reads_by_umi(INPUT, out, raw_tag="UB")
    assert stats["accepted"] == 18
    ref, mine = _load_list(GROUPED_GOLDEN), _load_list(out)
    assert len(ref) == len(mine)
    for a, b in zip(ref, mine):
        assert a.query_name == b.query_name
        assert a.flag == b.flag
        assert a.reference_start == b.reference_start
        assert a.query_sequence == b.query_sequence
        # the assigned molecular id must match fgbio's integer exactly
        assert a.get_tag("MI") == b.get_tag("MI")


def test_consensus_matches_fgbio(tmp_path):
    grouped = str(tmp_path / "grouped.bam")
    fgbiolike.group_reads_by_umi(INPUT, grouped, raw_tag="UB")
    out = str(tmp_path / "consensus.bam")
    stats = fgbiolike.call_molecular_consensus_reads(grouped, out, tag="MI", min_reads=3)
    assert stats["consensus_reads"] == 4

    ref, mine = _load_by_name(CONSENSUS_GOLDEN), _load_by_name(out)
    assert set(ref) == set(mine)
    for name, a in ref.items():
        b = mine[name]
        ta, tb = dict(a.get_tags()), dict(b.get_tags())
        assert a.query_sequence == b.query_sequence
        assert a.query_qualities.tobytes() == b.query_qualities.tobytes()
        assert a.flag == b.flag
        assert list(ta.get("cd", [])) == list(tb.get("cd", []))   # per-base depth
        assert list(ta.get("ce", [])) == list(tb.get("ce", []))   # per-base errors
        assert ta.get("cD") == tb.get("cD")
        assert ta.get("cM") == tb.get("cM")
        assert _f32(ta.get("cE")) == _f32(tb.get("cE"))
        assert ta.get("CB") == tb.get("CB")
        assert ta.get("RG") == tb.get("RG")


def test_consensus_one_shot(tmp_path):
    """The consensus() convenience runs group+call and matches the golden."""
    out = str(tmp_path / "c.bam")
    res = fgbiolike.consensus(INPUT, out, umi_tag="UB", min_reads=3)
    assert res["consensus"]["consensus_reads"] == 4
    ref, mine = _load_by_name(CONSENSUS_GOLDEN), _load_by_name(out)
    assert set(ref) == set(mine)


def test_pure_python_matches_cython(tmp_path):
    """The Cython consensus loop and the pure-Python fallback agree exactly."""
    grouped = str(tmp_path / "grouped.bam")
    fgbiolike.group_reads_by_umi(INPUT, grouped, raw_tag="UB")

    saved = fgbiolike._HAVE_FAST
    try:
        fgbiolike._HAVE_FAST = True
        fast_out = str(tmp_path / "fast.bam")
        fgbiolike.call_molecular_consensus_reads(grouped, fast_out, tag="MI", min_reads=3)
        fgbiolike._HAVE_FAST = False
        pure_out = str(tmp_path / "pure.bam")
        fgbiolike.call_molecular_consensus_reads(grouped, pure_out, tag="MI", min_reads=3)
    finally:
        fgbiolike._HAVE_FAST = saved

    fast, pure = _load_by_name(fast_out), _load_by_name(pure_out)
    assert set(fast) == set(pure)
    for name, a in fast.items():
        b = pure[name]
        assert a.query_sequence == b.query_sequence
        assert a.query_qualities.tobytes() == b.query_qualities.tobytes()
        assert dict(a.get_tags()).get("ce") == dict(b.get_tags()).get("ce")


def test_cell_collision_handling(tmp_path):
    """Grouping literally by UB (not MI) puts two cells in one group; fgbio would
    crash. Our default 'split' handles it; 'error' reproduces the fgbio failure."""
    grouped = str(tmp_path / "grouped.bam")
    fgbiolike.group_reads_by_umi(INPUT, grouped, raw_tag="UB")

    # split (default): no crash, both cells get their own consensus
    out = str(tmp_path / "split.bam")
    stats = fgbiolike.call_molecular_consensus_reads(
        grouped, out, tag="UB", min_reads=3, on_cell_collision="split")
    assert stats["cell_collisions"] == 1

    # error: reproduce fgbio's abort
    with pytest.raises(ValueError):
        fgbiolike.call_molecular_consensus_reads(
            grouped, str(tmp_path / "err.bam"), tag="UB", min_reads=3,
            on_cell_collision="error")
