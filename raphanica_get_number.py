#!/usr/bin/env python3

'''
Get the number of expanded gene families in L. raphanica genomes
  - Output: Gene list of Lenrap1
Byoungnam Min. Jan 8, 2021
'''

import os
import re
from collections import Counter
from argparse import ArgumentParser

RAPHANICA_GENOMES = set([
    'Lenra_T_1', 'Lenrap1', 'Lenrap1_155', 'LraINPA1820_1', 'LraJLM1587_1',
    'LraTFB9207_1'])


def main():
    '''Main funciton'''
    argparse_usage = 'parse_orthofinder.py -p <protein_dir> -o <orthogroups_txt>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-i', '--orthogroups_file', nargs=1, required=True,
        help='Orthogroups.txt file from the OrthoFinder')
    
    args = parser.parse_args()
    orthogroups_file = os.path.abspath(args.orthogroups_file[0])
    
    # Run functions :)
    parse_ortho(orthogroups_file)


def parse_ortho(ortho_file):
    '''Parse ortho file'''
    ortho_txt = import_file(ortho_file)
    num_gene_fams = 0
    c_count = Counter()
    for line in ortho_txt:
        line = re.sub(r'^OG\d+: ', '', line)
        genomes = [x.split('|')[1] for x in line.split(' ')]
        if set(genomes) != RAPHANICA_GENOMES:
            continue
        num_gene_fams += 1
        c_count += Counter(genomes)
    print('The number of raphanica-specific orthogroups: {}'.format(
        num_gene_fams))
    print(c_count)



def import_file(input_file):
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


if __name__ == '__main__':
    main()
