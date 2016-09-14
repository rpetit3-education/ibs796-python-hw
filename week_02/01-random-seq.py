#! /usr/bin/env python
"""
Generate random DNA sequences.

This script produces random DNA sequences with a bias towards the given
GC_CONTENT value.

Example Usage: ./01-random-seq.py TOTAL_READS GC_CONTENT READ_LENGTH
"""
import random
random.seed(123456)


def generate_sequence(expected_gc, length):
    """Create a random sequence of a given length based on an expected % GC."""
    seq = []
    gc = 0.0
    for i in xrange(length):
        if random.random() > gc_expected:
            if random.random() >= 0.50:
                seq.append('A')
            else:
                seq.append('T')
        else:
            gc += 1
            if random.random() >= 0.50:
                seq.append('G')
            else:
                seq.append('C')

    return [''.join(seq), gc / length]


if __name__ == '__main__':
    import sys
    try:
        total = int(sys.argv[1])
        gc_expected = float(sys.argv[2])
        read_length = int(sys.argv[3])
    except IndexError:
        print("Error, missing parameters. Please correct.")
        print("Usage: ./01-random-seq.py TOTAL_READS GC_CONTENT READ_LENGTH")
        sys.exit(1)

    for i in xrange(total):
        seq, gc_obs = generate_sequence(gc_expected, read_length)
        print(">{0} gc={1} length={2}\n{3}".format(
            i + 1, gc_obs, read_length, seq
        ))
