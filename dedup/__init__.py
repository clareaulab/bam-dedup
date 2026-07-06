"""
bam-dedup: a fast, JVM-free toolkit for removing duplicate reads from BAM/CRAM.

It handles the two major classes of duplicate with faithful, independent ports of
the reference tools, both accelerated in Cython (:mod:`dedup._fast`) and doing all
I/O through pysam (htslib) -- no Java dependency:

* **PCR / optical duplicates** -- :mod:`dedup.picardlike`, a port of Picard/htsjdk
  MarkDuplicates (coordinate-sorted, Illumina paired-end, default options). Marks
  or removes duplicate *records*.

* **Molecular (UMI) duplicates** -- :mod:`dedup.fgbiolike`, a port of fgbio
  GroupReadsByUmi + CallMolecularConsensusReads (Identity strategy, fragment
  path). Collapses reads sharing a UMI into one error-corrected *consensus* read.

Typical use::

    from dedup import picardlike, fgbiolike

    # PCR/optical: flag duplicates (records kept, 0x400 flag set)
    picardlike.mark_duplicates("input.bam", "marked.bam")

    # molecular/UMI: group by UMI and call consensus reads
    fgbiolike.consensus("input.bam", "consensus.bam", umi_tag="UB", min_reads=3)
"""

from dedup.picardlike import mark_duplicates, deduplicate
from dedup.fgbiolike import (
    group_reads_by_umi,
    call_molecular_consensus_reads,
    consensus,
)

try:
    from dedup import _fast  # noqa: F401
    HAVE_CYTHON = True
except ImportError:  # pragma: no cover
    HAVE_CYTHON = False

__version__ = "0.2.0"
__all__ = [
    "picardlike", "fgbiolike",
    "mark_duplicates", "deduplicate",
    "group_reads_by_umi", "call_molecular_consensus_reads", "consensus",
    "HAVE_CYTHON", "__version__",
]
