#!/usr/bin/env python3
"""
fgbiolike.py  --  a faithful, dependency-light re-implementation of the core of
two fgbio UMI tools, ``GroupReadsByUmi`` and ``CallMolecularConsensusReads``.

Where :mod:`dedup.picardlike` removes **PCR/optical** duplicates (Picard-style,
by flagging or dropping records), this module handles the other major class of
duplicate: **molecular (UMI) duplicates**. Reads that share a unique molecular
identifier are grouped into their source molecule and collapsed into a single,
error-corrected *consensus* read. This is the deduplication model used by, e.g.,
single-cell / duplex UMI protocols such as MAESTER.

It reproduces fgbio's decisions using the same algorithm:

    GroupReadsByUmi (strategy=Identity, edits=0):
        1. Filter secondary/supplementary, non-PF, unmapped, low-mapq, and
           N-containing-UMI reads.
        2. Sort survivors into fgbio "template-coordinate" order.
        3. Assign a monotonically increasing molecular id (MI) to each distinct
           (ref, 5' unclipped pos, strand, library, cell-barcode, UMI) group.

    CallMolecularConsensusReads (vanilla, fragment path):
        4. Group *consecutive* reads by the molecular tag.
        5. End-trim/mask each read by base quality, keep the most common
           alignment, and call a Bayesian consensus base + quality per position
           (pre/post-UMI error model), emitting cD/cM/cE/cd/ce tags.

The only runtime dependency is **pysam** (htslib -- no JVM). The per-base
consensus likelihood loop is accelerated in Cython (:mod:`dedup._fast`) with a
bit-identical pure-Python fallback.

SCOPE / faithful to fgbio for:
    * strategy=Identity, edits=0 grouping
    * single-end / fragment consensus (the MAESTER-style input)
    * default error model (pre=45, post=40), min-input-base-quality=10,
      per-base tags on, read-group "A"

NOT (yet) handled:
    * Adjacency/Edit/Paired UMI strategies (edits > 0)
    * paired-end / duplex consensus, overlapping-base consensus, quality trimming
"""

import argparse
import array
import math
import struct
import sys

import pysam

try:
    from dedup import _fast
    # guard against an older _fast built before the consensus loop was added
    _HAVE_FAST = hasattr(_fast, "consensus_call_molecule")
except ImportError:  # pragma: no cover - pure-Python fallback
    _fast = None
    _HAVE_FAST = False

INT_MAX = 2 ** 31 - 1
NOCALL = ord("N")
SHORT_MAX = 32767
DNA = b"ACGT"

# pysam cigar op codes
_OP_M, _OP_I, _OP_D, _OP_N, _OP_S, _OP_H, _OP_P, _OP_EQ, _OP_X = range(9)
_CONSUMES_QUERY = frozenset((_OP_M, _OP_I, _OP_S, _OP_EQ, _OP_X))
_CLIPS = frozenset((_OP_S, _OP_H))
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

# read-type
FRAGMENT, FIRST_OF_PAIR, SECOND_OF_PAIR = 0, 1, 2


# ============================================================================
# Log-probability arithmetic (fgbio NumericTypes.LogProbability / PhredScore).
# Natural-log space. libm is used for log/exp; results agree well within the
# 1e-3 phred rounding buffer that fixes the emitted integer quality.
# ============================================================================
_NEG_INF = float("-inf")
_LN10 = math.log(10.0)
_LN2 = math.log(2.0)
_LN3 = math.log(3.0)
_LN4 = math.log(4.0)
_LOG_FOUR_THIRDS = _LN4 - _LN3
_PHRED_MIN = 2
_PHRED_MAX = 93
_PRECISION = 0.001


def _from_phred(s):
    if s > 127:
        s = 127
    return _LN10 * s / -10.0


_MAX_VALUE_AS_LOG = _from_phred(_PHRED_MAX)


def _log1pexp(v):
    if v <= -37:
        return math.exp(v)
    elif v <= 18:
        return math.log1p(math.exp(v))
    elif v <= 33.3:
        return v + math.exp(-v)
    return v


def _log1mexp(v):
    if v <= _LN2:
        return math.log(-math.expm1(-v))
    return math.log1p(-math.exp(-v))


def _lor(a, b):
    if a == _NEG_INF:
        return b
    if b == _NEG_INF:
        return a
    if b < a:
        a, b = b, a
    return a + _log1pexp(b - a)


def _lor_array(vals):
    if all(v == _NEG_INF for v in vals):
        return _NEG_INF
    min_i = min(range(len(vals)), key=lambda i: vals[i])
    s = vals[min_i]
    for i in range(len(vals)):
        if i != min_i:
            s = _lor(s, vals[i])
    return s


def _a_or_not_b(a, b):
    if b == _NEG_INF:
        return a
    if a == b:
        return _NEG_INF
    if a < b:
        raise ValueError("Subtraction will be less than zero.")
    return a + _log1mexp(a - b)


def _lnot(a):
    if 0.0 < a:
        return _NEG_INF
    return _a_or_not_b(0.0, a)


def _prob_error_two_trials(a, b):
    if a < b:
        a, b = b, a
    if a - b >= 6:
        return a
    return _a_or_not_b(_lor(a, b), _LOG_FOUR_THIRDS + a + b)


def _phred_from_logprob(lp):
    if lp < _MAX_VALUE_AS_LOG:
        return _PHRED_MAX
    return int(math.floor(-10.0 * (lp / _LN10) + _PRECISION))


def _cap_phred(q):
    if q < _PHRED_MIN:
        return _PHRED_MIN
    if q > _PHRED_MAX:
        return _PHRED_MAX
    return q


def _float32(x):
    """Round a Python float to IEEE-754 single precision (fgbio stores cE as float)."""
    return struct.unpack("f", struct.pack("f", x))[0]


class ConsensusModel:
    """Precomputes the qual->logprob lookup tables for fixed pre/post error rates."""

    __slots__ = ("ln_pre", "ln_post", "p_err_third", "p_truth")

    def __init__(self, error_rate_pre_umi=45, error_rate_post_umi=40):
        self.ln_pre = _from_phred(error_rate_pre_umi)
        self.ln_post = _from_phred(error_rate_post_umi)
        self.p_err_third = []
        self.p_truth = []
        for q in range(0, 127):
            e = _prob_error_two_trials(self.ln_post, _from_phred(q))
            self.p_err_third.append(e - _LN3)   # normalizeByScalar(e, 3)
            self.p_truth.append(_lnot(e))

    def call_position(self, bases_at_pos, quals_at_pos, min_reads, min_consensus_q):
        ll = [0.0, 0.0, 0.0, 0.0]
        comp = [0.0, 0.0, 0.0, 0.0]
        obs = [0, 0, 0, 0]
        pet, ptt = self.p_err_third, self.p_truth
        for b, q in zip(bases_at_pos, quals_at_pos):
            if b == NOCALL:
                continue
            pe, pt = pet[q], ptt[q]
            for i in range(4):
                term = pt if b == DNA[i] else pe
                y = term - comp[i]
                t = ll[i] + y
                comp[i] = (t - ll[i]) - y
                ll[i] = t
            if b == 0x41:
                obs[0] += 1
            elif b == 0x43:
                obs[1] += 1
            elif b == 0x47:
                obs[2] += 1
            elif b == 0x54:
                obs[3] += 1

        depth = obs[0] + obs[1] + obs[2] + obs[3]
        raw_base, raw_qual = self._call(ll)
        if raw_base == NOCALL:
            errors = depth
        else:
            errors = depth - obs[DNA.index(raw_base)]
        if depth < min_reads:
            return NOCALL, 0, depth, errors
        if raw_qual < min_consensus_q:
            return NOCALL, 2, depth, errors
        return raw_base, raw_qual, depth, errors

    def _call(self, ll):
        ll_sum = _lor_array(ll)
        max_v = -1.7976931348623157e308
        max_i = -1
        assigned = False
        for i in range(4):
            v = ll[i]
            if (not assigned) or v > max_v:
                max_v = v
                max_i = i
                assigned = True
            elif abs(v - max_v) <= 2.0 ** -52:
                max_i = -1
        if max_i == -1:
            return NOCALL, _PHRED_MIN
        p_cons_err = _lnot(max_v - ll_sum)
        p = _prob_error_two_trials(self.ln_pre, p_cons_err)
        return DNA[max_i], _cap_phred(_phred_from_logprob(p))


# ============================================================================
# GroupReadsByUmi (strategy=Identity, edits=0)
# ============================================================================
def _clip_lengths(cigartuples):
    lead = trail = 0
    n = len(cigartuples)
    i = 0
    while i < n and cigartuples[i][0] in _CLIPS:
        lead += cigartuples[i][1]
        i += 1
    j = n - 1
    while j >= 0 and cigartuples[j][0] in _CLIPS:
        trail += cigartuples[j][1]
        j -= 1
    return lead, trail


def _unclipped_5prime(rec):
    """Strand-aware 5' unclipped position (1-based), matching fgbio/htsjdk."""
    lead, trail = _clip_lengths(rec.cigartuples)
    if rec.is_reverse:
        return rec.reference_end + trail
    return (rec.reference_start + 1) - lead


def _library_of(rec, rg_to_lib):
    try:
        rg = rec.get_tag("RG")
    except KeyError:
        return "unknown"
    return rg_to_lib.get(rg, "unknown")


def group_reads_by_umi(input_bam, output_bam,
                       raw_tag="UB", assign_tag="MI", cell_tag="CB",
                       min_map_q=1, include_non_pf=False, logger=None):
    """Port of ``fgbio GroupReadsByUmi -s Identity -e 0``.

    Writes ``output_bam`` in template-coordinate order with ``assign_tag`` (MI)
    set on every record. Returns a dict of counters.
    """
    log = logger or (lambda *a: None)
    src = pysam.AlignmentFile(input_bam, "rb")
    header = src.header.to_dict()
    rg_to_lib = {rg["ID"]: rg.get("LB", "unknown") for rg in header.get("RG", [])}

    filtered_poor = filtered_ns = filtered_non_pf = kept = 0
    records = []
    for rec in src:
        if rec.is_secondary or rec.is_supplementary:
            continue
        if not include_non_pf and rec.is_qcfail:
            filtered_non_pf += 1
            continue
        mapped = not rec.is_unmapped
        mate_mapped = rec.is_paired and not rec.mate_is_unmapped
        if not (mapped or mate_mapped):
            filtered_poor += 1
            continue
        if mapped and rec.mapping_quality < min_map_q:
            filtered_poor += 1
            continue
        try:
            umi = rec.get_tag(raw_tag)
        except KeyError:
            umi = None
        if umi is not None and "N" in umi:
            filtered_ns += 1
            continue

        ref1 = rec.reference_id if mapped else INT_MAX
        pos1 = _unclipped_5prime(rec) if mapped else INT_MAX
        neg1 = 1 if rec.is_reverse else 0
        cid = rec.get_tag(cell_tag) if (cell_tag and rec.has_tag(cell_tag)) else ""
        mid = (umi or "").upper()
        lib = _library_of(rec, rg_to_lib)
        # TemplateCoordinateKey.compare, fragment specialization
        key = (ref1, pos1, neg1, len(cid), cid, len(mid), mid, rec.query_name, lib)
        records.append((key, rec))
        kept += 1

    src.close()
    log("Accepted {:,} SAM records for grouping.".format(kept))

    records.sort(key=lambda t: t[0])

    out_header = _grouped_output_header(header)
    dst = pysam.AlignmentFile(output_bam, "wb", header=out_header)
    mi = 0
    prev_mol = None
    for key, rec in records:
        ref1, pos1, neg1, _cl, cid, _ml, mid, _name, lib = key
        mol = (ref1, pos1, neg1, lib, cid, mid)
        if prev_mol is None or mol != prev_mol:
            if prev_mol is not None:
                mi += 1
            prev_mol = mol
        rec.set_tag(assign_tag, str(mi), value_type="Z")
        dst.write(rec)
    dst.close()

    return {"accepted": kept, "discarded_non_pf": filtered_non_pf,
            "discarded_poor_alignment": filtered_poor,
            "discarded_ns_in_umi": filtered_ns,
            "molecules": (mi + 1) if prev_mol is not None else 0}


def _grouped_output_header(in_header):
    hd = dict(in_header)
    hd_line = dict(hd.get("HD", {"VN": "1.6"}))
    hd_line.update({"VN": "1.6", "SO": "unsorted", "GO": "query",
                    "SS": "unsorted:template-coordinate"})
    hd["HD"] = hd_line
    return pysam.AlignmentHeader.from_dict(hd)


# ============================================================================
# CallMolecularConsensusReads (vanilla, fragment path)
# ============================================================================
class SourceRead:
    __slots__ = ("id", "bases", "quals", "cigar", "sam")

    def __init__(self, id, bases, quals, cigar, sam):
        self.id = id
        self.bases = bases
        self.quals = quals
        self.cigar = cigar
        self.sam = sam

    @property
    def length(self):
        return len(self.bases)


def _reverse_cigar(cig):
    return list(reversed(cig))


def _truncate_to_query_length(cig, n):
    out = []
    pos = 1
    for op, ln in cig:
        if pos > n:
            break
        if op in _CONSUMES_QUERY:
            m = n - pos + 1
            out.append((op, ln if ln <= m else m))
            pos += ln
        else:
            out.append((op, ln))
    return out


def _coalesce(cig):
    out = []
    for op, ln in cig:
        if out and out[-1][0] == op:
            out[-1] = (op, out[-1][1] + ln)
        else:
            out.append((op, ln))
    return out


def _simplify_cigar(cig):
    if all(op in (_OP_M, _OP_I, _OP_D) for op, _ in cig):
        return _coalesce(cig)
    mapped = [((_OP_M, ln) if op in (_OP_S, _OP_EQ, _OP_X, _OP_H) else (op, ln))
              for op, ln in cig]
    return _coalesce(mapped)


def _is_prefix_of(this, that):
    if len(that) < len(this):
        return False
    last = len(this) - 1
    for i in range(len(this)):
        lop, llen = this[i]
        rop, rlen = that[i]
        if lop != rop:
            return False
        if i == last:
            if llen > rlen:
                return False
        elif llen != rlen:
            return False
    return True


def _cigar_order_key(cig):
    return ([(ln, op) for op, ln in cig], len(cig))


def to_source_read(rec, tag, min_base_q):
    seq = rec.query_sequence
    if seq is None:
        return None
    bases = bytearray(seq, "ascii")
    quals = list(rec.query_qualities)
    cig = rec.cigartuples or []
    if rec.is_reverse:
        bases = bytearray(bytes(bases).translate(_COMPLEMENT)[::-1])
        quals = quals[::-1]
        cig = _reverse_cigar(cig)

    for i in range(len(bases)):
        if quals[i] < min_base_q:
            bases[i] = NOCALL
            quals[i] = 2

    index = len(bases) - 1
    while index >= 0 and bases[index] == NOCALL:
        index -= 1
    length = index + 1

    mol_id = rec.get_tag(tag)
    if length == 0:
        return None
    if length == len(bases):
        return SourceRead(mol_id, bytes(bases), quals, cig, rec)
    return SourceRead(mol_id, bytes(bases[:length]), quals[:length],
                      _truncate_to_query_length(cig, length), rec)


def filter_to_most_common_alignment(recs):
    if len(recs) < 2:
        return recs
    n = len(recs)
    order = sorted(range(n), key=lambda i: -recs[i].length)
    groups = []  # [cigar, set(si), size]
    for si in range(len(order)):
        simple = _simplify_cigar(recs[order[si]].cigar)
        found = False
        for g in groups:
            if _is_prefix_of(simple, g[0]):
                g[1].add(si)
                g[2] += 1
                found = True
        if not found:
            groups.append([simple, {si}, 1])
    if not groups:
        return []
    best = None
    for g in groups:
        key = (g[2], _cigar_order_key(g[0]))
        if best is None or key[0] > best[1][0] or (key[0] == best[1][0] and key[1] < best[1][1]):
            best = (g, key)
    kept = {order[si] for si in best[0][1]}
    return [recs[i] for i in range(n) if i in kept]


def _consensus_read_length(reads, min_reads):
    return sorted((r.length for r in reads), reverse=True)[min_reads - 1]


def _call_molecule_py(model, capped, length, min_reads, min_consensus_q,
                      bases_out, quals_out, depths_out, errors_out):
    for pos in range(length):
        bs, qs = [], []
        for r in capped:
            if r.length > pos:
                bs.append(r.bases[pos])
                qs.append(r.quals[pos])
        base, qual, depth, errors = model.call_position(bs, qs, min_reads, min_consensus_q)
        bases_out[pos] = base
        quals_out[pos] = qual
        depths_out[pos] = depth if depth <= SHORT_MAX else SHORT_MAX
        errors_out[pos] = errors if errors <= SHORT_MAX else SHORT_MAX


def consensus_call(model, reads, min_reads, min_consensus_q):
    if len(reads) < min_reads:
        return None
    length = _consensus_read_length(reads, min_reads)
    bases = bytearray(length)
    quals = array.array("b", bytes(length))
    depths = array.array("h", bytes(2 * length))
    errors = array.array("h", bytes(2 * length))
    if _HAVE_FAST:
        _fast.consensus_call_molecule(model.p_err_third, model.p_truth, model.ln_pre,
                                      reads, length, min_reads, min_consensus_q,
                                      bases, quals, depths, errors)
    else:
        _call_molecule_py(model, reads, length, min_reads, min_consensus_q,
                          bases, quals, depths, errors)
    return bases, quals, depths, errors, reads[0].id


def _consensus_output_header(in_header, read_group_id):
    old_rgs = in_header.get("RG", [])

    def collapse(field):
        vals = []
        for rg in old_rgs:
            v = rg.get(field)
            if v is not None and v not in vals:
                vals.append(v)
        return ",".join(vals) if vals else None

    rg = {"ID": read_group_id}
    for f in ("DS", "LB", "SM", "PL", "PU", "CN"):
        v = collapse(f)
        if v is not None:
            rg[f] = v
    out = {"HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"}, "RG": [rg],
           "CO": ["Read group {} contains consensus reads generated from {} input "
                  "read groups.".format(read_group_id, len(old_rgs))]}
    return pysam.AlignmentHeader.from_dict(out)


def _emit_record(header, prefix, mol_id, read_type, cell_tag, cell_barcode,
                 bases, quals, depths, errors, per_base_tags):
    rec = pysam.AlignedSegment(header)
    rec.query_name = "{}:{}".format(prefix, mol_id)
    rec.flag = 0
    rec.is_unmapped = True
    if read_type == FIRST_OF_PAIR:
        rec.is_paired = True
        rec.is_read1 = True
        rec.mate_is_unmapped = True
    elif read_type == SECOND_OF_PAIR:
        rec.is_paired = True
        rec.is_read2 = True
        rec.mate_is_unmapped = True
    rec.reference_id = -1
    rec.reference_start = -1
    rec.query_sequence = bases.decode("ascii")
    rec.query_qualities = array.array("B", quals)
    rec.set_tag("RG", "A", value_type="Z")
    rec.set_tag("MI", str(mol_id), value_type="Z")
    if cell_tag and cell_barcode is not None:
        rec.set_tag(cell_tag, cell_barcode, value_type="Z")
    rec.set_tag("cD", int(max(depths)), value_type="i")
    rec.set_tag("cM", int(min(depths)), value_type="i")
    rec.set_tag("cE", _float32(sum(errors) / _float32(sum(depths))), value_type="f")
    if per_base_tags:
        rec.set_tag("cd", array.array("h", depths))
        rec.set_tag("ce", array.array("h", errors))
    return rec


def call_molecular_consensus_reads(input_bam, output_bam, tag="MI", min_reads=3,
                                   error_rate_pre_umi=45, error_rate_post_umi=40,
                                   min_input_base_quality=10, cell_tag="CB",
                                   read_group_id="A", read_name_prefix="",
                                   output_per_base_tags=True,
                                   on_cell_collision="split", logger=None):
    """Port of ``fgbio CallMolecularConsensusReads`` (fragment path).

    Groups consecutive records by ``tag`` and writes unmapped consensus reads.
    ``on_cell_collision`` controls what happens when a molecular group spans more
    than one cell barcode (fgbio *aborts*; we default to ``split``): one of
    ``"split"`` (per cell), ``"merge"`` (single consensus, no CB), ``"error"``.
    """
    log = logger or (lambda *a: None)
    model = ConsensusModel(error_rate_pre_umi, error_rate_post_umi)
    min_consensus_q = 2

    src = pysam.AlignmentFile(input_bam, "rb", check_sq=False)
    out_header = _consensus_output_header(src.header.to_dict(), read_group_id)
    dst = pysam.AlignmentFile(output_bam, "wb", header=out_header)

    stats = {"total_reads": 0, "consensus_reads": 0, "cell_collisions": 0}

    def emit(recs):
        cell_barcode = None
        if cell_tag:
            bcs = []
            for r in recs:
                if r.has_tag(cell_tag):
                    b = r.get_tag(cell_tag)
                    if b not in bcs:
                        bcs.append(b)
            cell_barcode = bcs[0] if len(bcs) == 1 else None

        fragments = [r for r in recs if not r.is_paired]
        first = [r for r in recs if r.is_paired and r.is_read1]
        second = [r for r in recs if r.is_paired and r.is_read2]

        out = _consensus_from_sam(fragments)
        if out is not None:
            dst.write(_emit_record(out_header, read_name_prefix, out[4], FRAGMENT,
                                   cell_tag, cell_barcode, out[0], out[1], out[2],
                                   out[3], output_per_base_tags))
            stats["consensus_reads"] += 1

        if first or second:
            o1, o2 = _consensus_from_sam(first), _consensus_from_sam(second)
            if o1 is not None and o2 is not None:
                for o, rtype in ((o1, FIRST_OF_PAIR), (o2, SECOND_OF_PAIR)):
                    dst.write(_emit_record(out_header, read_name_prefix, o[4], rtype,
                                           cell_tag, cell_barcode, o[0], o[1], o[2],
                                           o[3], output_per_base_tags))
                    stats["consensus_reads"] += 1

    def _consensus_from_sam(records):
        if len(records) < min_reads:
            return None
        source = [sr for sr in (to_source_read(r, tag, min_input_base_quality)
                                for r in records) if sr is not None]
        filtered = filter_to_most_common_alignment(source)
        if len(filtered) < min_reads:
            return None
        return consensus_call(model, filtered, min_reads, min_consensus_q)

    def process(group):
        if not group:
            return
        batches = [group]
        if cell_tag:
            bcs = []
            for r in group:
                if r.has_tag(cell_tag):
                    b = r.get_tag(cell_tag)
                    if b not in bcs:
                        bcs.append(b)
            if len(bcs) > 1:
                stats["cell_collisions"] += 1
                if on_cell_collision == "error":
                    raise ValueError(
                        "Multiple different cell barcodes found for tag {}: {}"
                        .format(cell_tag, bcs))
                if on_cell_collision == "split":
                    by_cell = {}
                    for r in group:
                        key = r.get_tag(cell_tag) if r.has_tag(cell_tag) else None
                        by_cell.setdefault(key, []).append(r)
                    batches = list(by_cell.values())
        for b in batches:
            emit(b)

    cur = object()
    group = []
    for rec in src:
        stats["total_reads"] += 1
        try:
            k = rec.get_tag(tag)
        except KeyError:
            k = None
        if k != cur:
            process(group)
            group = []
            cur = k
        group.append(rec)
    process(group)

    dst.close()
    src.close()
    log("Total raw reads considered: {:,}. Consensus reads emitted: {:,}."
        .format(stats["total_reads"], stats["consensus_reads"]))
    return stats


def consensus(input_bam, output_bam, grouped_bam=None, umi_tag="UB", min_reads=3,
              consensus_tag="MI", cell_tag="CB", on_cell_collision="split",
              logger=None, **kwargs):
    """Convenience one-shot: group by UMI then call molecular consensus.

    Runs :func:`group_reads_by_umi` followed by :func:`call_molecular_consensus_reads`,
    the fgbio two-step pipeline. ``grouped_bam`` is an optional path for the
    intermediate; a temporary file is used if omitted. By default the consensus
    step groups by the assigned molecular id (``MI``) -- the correct, crash-free
    equivalent of fgbio ``-t UB`` for cell-barcoded data.
    """
    import os
    import tempfile

    tmp = grouped_bam
    if tmp is None:
        fd, tmp = tempfile.mkstemp(suffix=".grouped.bam")
        os.close(fd)
    try:
        s1 = group_reads_by_umi(input_bam, tmp, raw_tag=umi_tag, cell_tag=cell_tag,
                                logger=logger)
        s2 = call_molecular_consensus_reads(tmp, output_bam, tag=consensus_tag,
                                            min_reads=min_reads, cell_tag=cell_tag,
                                            on_cell_collision=on_cell_collision,
                                            logger=logger, **kwargs)
    finally:
        if grouped_bam is None and os.path.exists(tmp):
            os.remove(tmp)
    return {"group": s1, "consensus": s2}


# ============================================================================
# CLI
# ============================================================================
def main(argv=None):
    p = argparse.ArgumentParser(
        prog="bam-consensus",
        description="fgbio-style UMI grouping + molecular consensus deduplication.")
    sub = p.add_subparsers(dest="command", required=True)

    g = sub.add_parser("group", help="GroupReadsByUmi (Identity, edits=0)")
    g.add_argument("-i", "--input", required=True)
    g.add_argument("-o", "--output", required=True)
    g.add_argument("-t", "--umi-tag", default="UB")
    g.add_argument("-c", "--cell-tag", default="CB")

    c = sub.add_parser("call", help="CallMolecularConsensusReads")
    c.add_argument("-i", "--input", required=True)
    c.add_argument("-o", "--output", required=True)
    c.add_argument("-t", "--tag", default="MI")
    c.add_argument("-M", "--min-reads", type=int, default=3)
    c.add_argument("-c", "--cell-tag", default="CB")
    c.add_argument("--on-cell-collision", choices=["split", "merge", "error"], default="split")

    a = sub.add_parser("consensus", help="group + call in one step")
    a.add_argument("-i", "--input", required=True)
    a.add_argument("-o", "--output", required=True)
    a.add_argument("-t", "--umi-tag", default="UB")
    a.add_argument("-M", "--min-reads", type=int, default=3)
    a.add_argument("-c", "--cell-tag", default="CB")
    a.add_argument("--consensus-tag", default="MI")
    a.add_argument("--on-cell-collision", choices=["split", "merge", "error"], default="split")

    args = p.parse_args(argv)
    log = lambda m: sys.stderr.write(m + "\n")

    if args.command == "group":
        print(group_reads_by_umi(args.input, args.output, raw_tag=args.umi_tag,
                                 cell_tag=(args.cell_tag or None), logger=log), file=sys.stderr)
    elif args.command == "call":
        print(call_molecular_consensus_reads(args.input, args.output, tag=args.tag,
              min_reads=args.min_reads, cell_tag=(args.cell_tag or None),
              on_cell_collision=args.on_cell_collision, logger=log), file=sys.stderr)
    else:
        print(consensus(args.input, args.output, umi_tag=args.umi_tag,
              min_reads=args.min_reads, consensus_tag=args.consensus_tag,
              cell_tag=(args.cell_tag or None), on_cell_collision=args.on_cell_collision,
              logger=log), file=sys.stderr)


if __name__ == "__main__":
    main()
