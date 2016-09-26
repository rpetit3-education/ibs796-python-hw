#! /usr/bin/env python
"""
Extract protein sequences from GenBank file.
"""

if __name__ == '__main__':
    import argparse
    import gzip
    from Bio import SeqIO

    parser = argparse.ArgumentParser(
        prog='extract-cds.py',
        conflict_handler='resolve',
        description=('Generate summary statistics for a given FASTQ file.')
    )
    parser.add_argument('gb', metavar="GENBANK", type=str,
                        help='GENBANK file to parse.')
    parser.add_argument('dna', metavar="DNA_FASTA", type=str,
                        help='GENBANK file to parse.')
    parser.add_argument('aa', metavar="AA_FASTA", type=str,
                        help='GENBANK file to parse.')
    parser.add_argument('--gzip', action='store_true',
                        help='Input is compressed.')

    args = parser.parse_args()
    types = {}
    genbank_fh = gzip.open(args.gb, 'rb') if args.gzip else open(args.gb, 'r')
    with open(args.dna, 'w') as dna_fh, open(args.aa, 'w') as aa_fh:
        for record in SeqIO.parse(genbank_fh, "genbank"):
            for feature in record.features:
                if feature.type in ['tRNA', 'mRNA', 'rRNA', 'CDS']:
                    header = ">{0} product={1} type={2}".format(
                        feature.qualifiers['locus_tag'][0],
                        feature.qualifiers['product'][0],
                        feature.type
                    )
                    dna = feature.extract(record.seq)
                    dna_fh.write("{0}\n{1}\n".format(header, dna))
                    if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
                        aa_fh.write("{0}\n{1}\n".format(
                            header, dna.translate()
                        ))

