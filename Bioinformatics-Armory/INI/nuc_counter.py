#!/usr/bin/env python3

"""ROSALIND Intro to bioinformatics armory

This script takes a DNA sequence string and prints nucleotide counts.

The DNA string s length is at most 1000bp

The space-separated output contains four ints representing nucleotides A, C, G, T.
"""

from Bio.Seq import Seq

import argparse

def parse_args():
    parser = argparse.ArgumentParser(prog = 'nuc_counter.py', conflict_handler = 'resolve')
    parser.add_argument('-seq', type = str, required = False, help = '=> DNA sequence string')
    parser.add_argument('-seqfile', type = str, required = False, help = '=> path/to/DNA sequence string file')
    return(parser.parse_args())

def set_input(parsed_args) -> str:
    """If at least one arg exists, return sequence str.

    Note: if both exist, override with seqfile str.
    """
    if parsed_args.seq is None and parsed_args.seqfile is None:
        return(None)
    if parsed_args.seqfile is not None:
        with open(parsed_args.seqfile, 'r') as seqfile:
            return(seqfile.read().replace('\n',''))
    else:
        return(parsed_args.seq)

def count_nucs(DNA_string: str) -> int:
    nuc_tup = ("A", "C", "G", "T")
    count_str = ""
    for ind, nuc in enumerate(nuc_tup):
        if ind != len(nuc_tup)-1:
            count_str += f"{DNA_string.count(nuc)} "
        else:
            count_str += f"{DNA_string.count(nuc)}"
    return(count_str)

def main():
    # organize inputs from command line
    args = parse_args()
    input_seq = set_input(args)
    DNA_seq = Seq(input_seq)
    
    # impose limits
    if len(DNA_seq)-1 > 1000:
        print("Size limit of 1000 exceeded. Exiting.")
        return(None)

    # print counts
    print(count_nucs(DNA_seq))

if __name__ == '__main__':
    main()