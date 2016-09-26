#! /usr/bin/env python
"""
Generate random DNA sequences.

This script produces random DNA sequences with a bias towards the given
GC_CONTENT value.

Example Usage: ./01-random-seq.py TOTAL_READS GC_CONTENT READ_LENGTH
"""
import argparse
import random


def restricted_float(x):
    """
    Shamelessly taken from Stackoverflow Thanks Chepner

    http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


def generate_sequence(expected_gc, length, seed=None):
    """Create a random sequence of a given length based on an expected % GC."""
    if seed:
        random.seed(seed)
    seq = []
    gc = 0.0
    for i in xrange(length):
        if random.random() > expected_gc:
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
    parser = argparse.ArgumentParser(
        prog='random_seq.py',
        conflict_handler='resolve',
        description=('Generate random FASTA sequences.')
    )
    parser.add_argument('total', metavar="INT", type=int,
                        help='Total number of reads to simulate.')
    parser.add_argument('gc', metavar="FLOAT [0.0-1.0]", type=restricted_float,
                        help='An expected GC content.')
    parser.add_argument('length', metavar="INT", type=int,
                        help='Total length of each read.')
    parser.add_argument('--seed', metavar="INT", type=int, default=0,
                        help='A random seed value.')

    args = parser.parse_args()
    for i in xrange(args.total):
        seq, gc_obs = generate_sequence(args.gc, args.length, seed=args.seed)
        print(">{0} gc={1} length={2}\n{3}".format(
            i + 1, gc_obs, args.length, seq
        ))
