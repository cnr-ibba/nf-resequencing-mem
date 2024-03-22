#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script will write a list of region from a reference file relying on
coverage depth

Author: Paolo Cozzi
Date: 2024/03/21

Usage:
    python split_ref_by_depth.py --depth_file <depth_file>

Arguments:
    depth_file: The path of samtools depth output file

Example:
    python split_ref_by_depth.py --depth_file all-samples.depth.tsv
"""

import csv
import gzip
import logging
import argparse
from typing import List

DEFAULT_LOGGING_LEVEL = logging.INFO
MAX_LOGGING_LEVEL = logging.CRITICAL


def setup_logger(verbose_level):
    fmt=('%(levelname)s %(asctime)s [%(module)s:%(lineno)s %(funcName)s] :: '
            '%(message)s')
    logging.basicConfig(format=fmt, level=max((0, min((MAX_LOGGING_LEVEL,
                        DEFAULT_LOGGING_LEVEL-(verbose_level*10))))))


def __append(regions: List[List], region: List, overlap_size: int):
    """Resize the ending position of the last region and the starting position
    of the current region in order to create an overlap"""

    last_region = regions[-1]

    logging.debug(f"extending last region {last_region}")

    # this is the ending interval. If i extend both ending, I will get a
    # double overlap region
    last_region[2] += int(overlap_size / 2)

    logging.debug(f"last region is now {last_region}")

    # can this be greater than last region end?
    if last_region[2] > region[2]:
        last_region[2] = region[2]

    logging.debug(f"extending current region {region}")

    # this is the start position of the new interval
    region[1] -= int(overlap_size / 2)

    # check overlap sizes
    if region[1] < 1:
        region[1] = 1

    logging.debug(f"current region is now {region}")

    # this is a normal list extension
    regions.append(region)
    return regions


def append_or_extend(
        regions: List[List], region: List, min_length: int, overlap_size: int):
    """
    Appends or extends a region based on the minimum length.

    Args:
        regions (List[List]): A list of regions.
        region (List): The region to append or extend.
        min_length (int): The minimum length required to append a new region.
        overlap_size (int): overlap size between two regions

    Returns:
        List[List]: The updated list of regions.
    """
    # simply add a region if is the first one
    if len(regions) == 0:
        logging.debug("Simply adding the current region")
        regions.append(region)
        return regions

    logging.debug(f"Last region: {regions[-1]}")
    logging.debug(f"Received region: {region}")

    # get the previous region and the current position
    last_region = regions[-1]

    current_start = region[1]
    current_end = region[2]

    # determine the current size of the last region
    current_length = current_end - current_start

    if current_length >= min_length:
        logging.debug("append the new region to the regions list")
        regions = __append(regions, region, overlap_size)
    else:
        logging.debug("extend the last region with the new end position")
        last_region[2] = current_end
        logging.debug(f"last region is now {last_region}")

    # there are all memory references. However declare a value to return
    return regions


def split_ref_by_coverage(
        depthfile: str, max_coverage: int, min_length: int, overlap_size: int):
    with gzip.open(depthfile, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t", lineterminator="\n")
        header = next(reader)

        # header: ["#CHROM", "POS", "Sample_1.bam", "Sample_2.bam", ...]
        logging.debug(f"Got header: {header}")

        # Inizialize some variables
        line = next(reader)
        start_chrom = line[0]
        start_pos = int(line[1])
        old_pos = start_pos
        total_sum = 0

        regions = []

        logging.debug(f"Starting from chrom: '{start_chrom}'")

        for i, line in enumerate(reader):
            chrom, pos = line[0], int(line[1])

            if chrom != start_chrom:
                logging.debug(f"{i}: Got a new chromosome '{chrom}'")
                logging.debug(
                    f"{i}: Test for a new region with: "
                    f"{[start_chrom, start_pos, old_pos]}"
                    f" ({old_pos-start_pos} bp)"
                )

                # add and open a new region
                regions = append_or_extend(
                    regions,
                    [start_chrom, start_pos, old_pos],
                    min_length,
                    overlap_size
                )

                # reset variables
                start_chrom = chrom
                start_pos = pos
                total_sum = 0

            # determine region size
            length = pos - start_pos

            # determine cumulative coverage in this region
            coverage_at_position = sum(int(sample_cov) for sample_cov in line[2:])
            total_sum += coverage_at_position

            if total_sum > max_coverage and length >= min_length:
                logging.debug(
                    f"{i}: Test for a new region with: "
                    f"{[start_chrom, start_pos, old_pos]}"
                    f" ({old_pos-start_pos} bp)"
                )

                # add and open a new region
                regions = append_or_extend(
                    regions,
                    [start_chrom, start_pos, old_pos],
                    min_length,
                    overlap_size
                )

                # reset variables
                start_chrom = chrom
                start_pos = pos
                total_sum = 0

            # track last position (useful when chromosome changes)
            old_pos = pos

        # check for an open region
        if total_sum > 0:
            logging.debug(
                f"{i}: Test for last region with: "
                f"{[start_chrom, start_pos, old_pos]}"
                f" ({old_pos-start_pos} bp)"
            )

            # add and open a new region
            regions = append_or_extend(
                regions,
                [start_chrom, start_pos, old_pos],
                min_length,
                overlap_size
            )

    return regions


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-d', '--depth_file', required=True,
        help="The output of samtools depth file for all samples"
    )
    parser.add_argument(
        "--min_length", default=10_000, type=int,
        help="minimum fragment length in bp"
    )
    parser.add_argument(
        "--max_coverage", default=500_000_000, type=int,
        help=(
            "max cumulative coverage per region "
            "(the sum of all coverages for each positions)"
        )
    )
    parser.add_argument(
        "--overlap_size", default=1_000, type=int,
        help="Overlapping size between two adjacent regions in bp"
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help="increase logging verbosity"
    )
    parser.add_argument(
        '-q', '--quiet', action='count', default=0,
        help="decrease logging verbosity"
    )
    args = parser.parse_args()

    # setup logger
    setup_logger(verbose_level=args.verbose-args.quiet)
    logging.info(f"args: {args}")

    # split reference by coverage depth
    regions = split_ref_by_coverage(
        args.depth_file, args.max_coverage, args.min_length, args.overlap_size)

    logging.info(f"Number of regions: {len(regions)}")

    for chrom, start, end in regions:
        print(f"{chrom}:{start}-{end}")
