"""Command-line entry point for bam-dedup."""
import argparse

from dedup import __version__
from dedup.picardlike import mark_duplicates, DEFAULT_OPTICAL_PIXEL_DISTANCE


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="bam-dedup",
        description="Fast, JVM-free Picard MarkDuplicates reimplementation.")
    p.add_argument("-i", "--input", required=True,
                   help="coordinate-sorted input BAM")
    p.add_argument("-o", "--output", required=True, help="output BAM")
    p.add_argument("--remove-duplicates", action="store_true",
                   help="remove duplicate records instead of only flagging them")
    p.add_argument("--remove-sequencing-duplicates", action="store_true",
                   help="remove optical/sequencing duplicates")
    p.add_argument("--no-optical", action="store_true",
                   help="disable optical-duplicate detection (READ_NAME_REGEX=null)")
    p.add_argument("--optical-pixel-distance", type=int,
                   default=DEFAULT_OPTICAL_PIXEL_DISTANCE,
                   help="optical duplicate pixel distance (default: %(default)s)")
    p.add_argument("--version", action="version",
                   version="bam-dedup {}".format(__version__))
    args = p.parse_args(argv)

    mark_duplicates(
        args.input, args.output,
        remove_duplicates=args.remove_duplicates,
        remove_sequencing_duplicates=args.remove_sequencing_duplicates,
        read_name_regex_enabled=not args.no_optical,
        optical_pixel_distance=args.optical_pixel_distance,
        metrics_file=None,
    )


if __name__ == "__main__":
    main()
