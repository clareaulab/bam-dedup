"""
Tests for the barcode_tag option in dedup.picardlike.

barcode_tag is a deliberate design choice, not a port of Picard's own
BARCODE_TAG (see the module docstring in dedup/picardlike.py for why): it
groups duplicate candidates by (library, barcode_tag value, position,
orientation[, mate]), so reads that would otherwise look like duplicates are
only flagged as duplicates of one another if they also carry the same tag
value. This is what lets a single-cell pipeline dedup within a cell without
ever collapsing reads from different cells that happen to share a position.
"""
import pysam
import pytest

from dedup import picardlike

REF = "chrT"
REF_LEN = 1000


def _make_read(header, name, pos, barcode=None, barcode_tag="XB", qual=37, length=50):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigarstring = "{}M".format(length)
    a.query_sequence = "A" * length
    a.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * length)
    if barcode is not None:
        a.set_tag(barcode_tag, barcode, value_type="Z")
    return a


def _write_bam(path, reads, header=None):
    header = header or pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for r in reads:
            bam.write(r)
    return header


def _dup_flags_by_name(path):
    with pysam.AlignmentFile(path, "rb") as bam:
        return {r.query_name: r.is_duplicate for r in bam.fetch(until_eof=True)}


def test_without_barcode_tag_all_four_collapse(tmp_path):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    reads = [
        _make_read(header, "aaaa_best", 100, barcode="AAAA", qual=40),
        _make_read(header, "aaaa_dup", 100, barcode="AAAA", qual=20),
        _make_read(header, "cccc_best", 100, barcode="CCCC", qual=39),
        _make_read(header, "cccc_dup", 100, barcode="CCCC", qual=20),
    ]
    src = str(tmp_path / "in.bam")
    _write_bam(src, reads, header)

    out = str(tmp_path / "out.bam")
    picardlike.mark_duplicates(src, out, read_name_regex_enabled=False)

    flags = _dup_flags_by_name(out)
    assert sum(flags.values()) == 3
    assert flags["aaaa_best"] is False  # highest score overall wins as the single survivor


def test_barcode_tag_keeps_groups_separate(tmp_path):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    reads = [
        _make_read(header, "aaaa_best", 100, barcode="AAAA", qual=40),
        _make_read(header, "aaaa_dup", 100, barcode="AAAA", qual=20),
        _make_read(header, "cccc_best", 100, barcode="CCCC", qual=39),
        _make_read(header, "cccc_dup", 100, barcode="CCCC", qual=20),
    ]
    src = str(tmp_path / "in.bam")
    _write_bam(src, reads, header)

    out = str(tmp_path / "out.bam")
    picardlike.mark_duplicates(src, out, read_name_regex_enabled=False, barcode_tag="XB")

    flags = _dup_flags_by_name(out)
    # One duplicate per barcode group -- never a cross-barcode collapse.
    assert sum(flags.values()) == 2
    assert flags["aaaa_best"] is False
    assert flags["aaaa_dup"] is True
    assert flags["cccc_best"] is False
    assert flags["cccc_dup"] is True


def test_barcode_tag_deduplicate_removes_within_group_only(tmp_path):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    reads = [
        _make_read(header, "aaaa_best", 100, barcode="AAAA", qual=40),
        _make_read(header, "aaaa_dup", 100, barcode="AAAA", qual=20),
        _make_read(header, "cccc_best", 100, barcode="CCCC", qual=39),
        _make_read(header, "cccc_dup", 100, barcode="CCCC", qual=20),
    ]
    src = str(tmp_path / "in.bam")
    _write_bam(src, reads, header)

    out = str(tmp_path / "out.bam")
    picardlike.deduplicate(src, out, read_name_regex_enabled=False, barcode_tag="XB")

    remaining = set(_dup_flags_by_name(out))
    assert remaining == {"aaaa_best", "cccc_best"}


def test_barcode_tag_missing_value_grouped_as_no_barcode(tmp_path):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    # Neither read carries the XB tag at all; both fall into the same "" bucket.
    reads = [
        _make_read(header, "no_tag_best", 100, barcode=None, qual=40),
        _make_read(header, "no_tag_dup", 100, barcode=None, qual=20),
    ]
    src = str(tmp_path / "in.bam")
    _write_bam(src, reads, header)

    out = str(tmp_path / "out.bam")
    picardlike.mark_duplicates(src, out, read_name_regex_enabled=False, barcode_tag="XB")

    flags = _dup_flags_by_name(out)
    assert flags == {"no_tag_best": False, "no_tag_dup": True}


def test_barcode_tag_none_is_bit_identical_to_default(tmp_path):
    """barcode_tag defaults to None, so existing callers see no behavior change."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF, "LN": REF_LEN}],
    })
    reads = [
        _make_read(header, "r1", 100, qual=40),
        _make_read(header, "r2", 100, qual=20),
    ]
    src = str(tmp_path / "in.bam")
    _write_bam(src, reads, header)

    out_default = str(tmp_path / "default.bam")
    out_explicit_none = str(tmp_path / "explicit_none.bam")
    picardlike.mark_duplicates(src, out_default, read_name_regex_enabled=False)
    picardlike.mark_duplicates(src, out_explicit_none, read_name_regex_enabled=False, barcode_tag=None)

    assert _dup_flags_by_name(out_default) == _dup_flags_by_name(out_explicit_none)
