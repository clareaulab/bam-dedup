"""
bam-dedup: a fast, JVM-free reimplementation of Picard MarkDuplicates.

The core duplicate-marking algorithm is a faithful port of Picard/htsjdk
MarkDuplicates (coordinate-sorted, Illumina paired-end, default options), with
the profiled hot spots accelerated in Cython (:mod:`dedup._fast`). It reads and
writes BAM/CRAM through pysam (htslib) -- no Java dependency.

Typical use::

    from dedup import picardlike

    # flag duplicates (Picard-style: records kept, 0x400 flag set)
    picardlike.mark_duplicates("input.bam", "marked.bam")

    # produce a deduplicated BAM (duplicate records removed)
    picardlike.deduplicate("input.bam", "dedup.bam")
"""

from dedup.picardlike import mark_duplicates, deduplicate

try:
    from dedup import _fast  # noqa: F401
    HAVE_CYTHON = True
except ImportError:  # pragma: no cover
    HAVE_CYTHON = False

__version__ = "0.1.0"
__all__ = ["picardlike", "mark_duplicates", "deduplicate", "HAVE_CYTHON", "__version__"]
