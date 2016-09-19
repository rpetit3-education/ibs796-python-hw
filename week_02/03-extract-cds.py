#! /usr/bin/env python
"""
Extract protein sequences from GenBank file.
"""

import gzip

CODON = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


def reset_cds():
    """Reset CDS variables."""
    return [False, False, None, None, None]


def parse_genbank(genbank):
    """Parse a genbank file."""
    is_feature = False
    is_sequence = False
    is_cds, is_pseudo, locus_tag, product, pos = reset_cds()

    seq = []
    cds = {}
    with gzip.open(genbank, 'r') as fh:
        for line in fh:
            if line.startswith('FEATURES'):
                is_feature = True
            elif line.startswith('ORIGIN'):
                is_sequence = True
                is_feature = False
            elif is_sequence and line.startswith('//'):
                is_sequence = False
                break
            elif is_feature:
                line = line.strip()
                if line.startswith('CDS'):
                    is_cds = True
                    pos = line.split()[1]
                elif is_cds and line.startswith('/locus_tag'):
                    locus_tag = line.replace('"', '').replace("/locus_tag=", '')
                elif is_cds and line.startswith('/pseudo'):
                    is_pseudo = True
                elif is_cds and line.startswith('/product'):
                    product = line.replace('"', '').replace('/product=', '')
                    if not is_pseudo:
                        cds[locus_tag] = {'pos': pos, 'product': product}
                    is_cds, is_pseudo, locus_tag, product, pos = reset_cds()
            elif is_sequence:
                line = line.strip()
                line = line.split(' ')[1:]
                seq.append(''.join(line).upper())
            else:
                pass

        return [''.join(seq), cds]


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([complement[b] for b in seq[::-1]])


def translate(seq):
    """Translate a DNA sequence into a Protein sequence."""
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    aa = []
    for codon in codons:
        aa.append(CODON[codon])
    return ''.join(aa)

if __name__ == '__main__':
    import sys

    sequence, cds_features = parse_genbank(sys.argv[1])
    for key, vals in cds_features.iteritems():
        start = None
        stop = None
        seq = None
        pos = vals['pos']
        if 'complement' in pos:
            pos = pos.replace('complement(', '').replace(')', '')
            start, stop = pos.split('..')
            seq = reverse_complement(sequence[int(start) - 1:int(stop)])
        else:
            start, stop = pos.split('..')
            seq = sequence[int(start) - 1:int(stop)]

        print('>{0} product={1}\n{2}'.format(
            key, vals['product'], translate(seq)
        ))
