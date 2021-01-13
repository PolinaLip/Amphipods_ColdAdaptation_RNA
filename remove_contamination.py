#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description='filter out contigs in assembly which are not from desired taxon')
    parser.add_argument('-i', type=argparse.FileType(), help='input tsv file with needed contigs')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output fasta file with filtered assembly')
    parser.add_argument('-a', type=argparse.FileType(), help='path to non-filtered assembly')

    args = parser.parse_args()

    outp = args.output

    names = set()
    for contig_name in args.i:
        contig_name = contig_name.split('\t')[1]
        contig_name = contig_name.split('.')[0]
        names.add(contig_name)
    print(len(names))
    for contig in args.a:
        contig = contig.strip()
        if contig[0] == '>':
            pickup = False
            contig = contig.split('>')[1]
            if contig in names:
                outp.write('>%s\n' % (contig))
                pickup = True
        elif pickup:
            outp.write('%s\n' % (contig))

if __name__ == '__main__':
    main()
~           
