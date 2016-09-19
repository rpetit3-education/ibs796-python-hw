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
    import sys
    import json
    from collections import Counter

    stats = parse_fastq(sys.stdin)
    per_base_stats = determine_per_base_quality(stats['per_base'])

    # Print the Stats
    len_counts = Counter(stats['lengths'])
    print(json.dumps(
        OrderedDict([
            ('total_reads', len(stats['lengths'])),
            ('total_bp', sum(stats['lengths'])),
            ('mean_length', "{0:.3f}".format(
             sum(stats['lengths']) / float(len(stats['lengths'])))),
            ('read_lengths', len_counts),
            ('mean_quality', "{0:.3f}".format(
                sum(stats['per_read']) / len(stats['per_read']) - 33
            )),
            ('per_base_quality', per_base_stats)
        ]),
        indent=4
    ))
