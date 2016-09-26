#! /usr/bin/env python
"""
Produce summary statistics of a FASTQ.

This script reads a FASTQ file from STDIN and outputs summary statistics
of the sequencing. Assumes the FASTQ is in 4-line per entry format.

Example Usage: zcat FASTQ.gz | ./02-fastq-stats.py
"""
from collections import OrderedDict

def parse_fastq(fh):
    """Read through a FASTQ and return stats."""
    stats = {'lengths': [], 'per_read': [], 'per_base': OrderedDict()}
    while 1:
        head = fh.readline()
        if not head:
            break
        fh.readline()  # Sequence
        fh.readline()  # Plus line
        qual = fh.readline().rstrip()
        vals = [ord(i) for i in qual]

        # Update everything
        length = float(len(qual))
        stats['lengths'].append(length)
        stats['per_read'].append(sum(vals) / length)

        # Per base quality
        for pos, val in enumerate(vals):
            if pos not in stats['per_base']:
                stats['per_base'][pos] = []
            stats['per_base'][pos].append(val)

    return stats


def determine_per_base_quality(bases):
    """Calculate the mean per base quality."""
    mean = OrderedDict()
    for pos, quals in bases.iteritems():
        mean[pos] = "{0:.3f}".format(sum(quals) / float(len(quals)) - 33)
    return mean


if __name__ == '__main__':
    import argparse
    import gzip
    import json
    from collections import Counter

    parser = argparse.ArgumentParser(
        prog='fastq_stats.py',
        conflict_handler='resolve',
        description=('Generate summary statistics for a given FASTQ file.')
    )
    parser.add_argument('fastq', metavar="FASTQ", type=str,
                        help='FASTQ file to parse.')
    parser.add_argument('--output', metavar="OUTPUT_NAME", type=str,
                        help='Write output to a file.')
    parser.add_argument('--compressed', action='store_true',
                        help='Input is compressed.')

    args = parser.parse_args()
    stats = None
    if args.compressed:
        with gzip.open(args.fastq, 'rb') as fh:
            stats = parse_fastq(fh)
    else:
        with open(args.fastq, 'r') as fh:
            stats = parse_fastq(fh)

    json_string = json.dumps(
        OrderedDict([
            ('total_reads', len(stats['lengths'])),
            ('total_bp', sum(stats['lengths'])),
            ('mean_length', "{0:.3f}".format(
                sum(stats['lengths']) / float(len(stats['lengths']))
            )),
            ('read_lengths', Counter(stats['lengths'])),
            ('mean_quality', "{0:.3f}".format(
                sum(stats['per_read']) / len(stats['per_read']) - 33
            )),
            ('per_base_quality', determine_per_base_quality(stats['per_base']))
        ]),
        indent=4
    )

    # Output the Stats
    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(json_string)
    else:
        print(json_string)
