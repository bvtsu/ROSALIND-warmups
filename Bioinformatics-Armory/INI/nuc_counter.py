"""ROSALIND Intro to bioinformatics armory

This script takes a DNA sequence string and prints nucleotide counts.

The DNA string s length is at most 1000bp

The space-separated output contains four ints representing nucleotides A, C, G, T.
"""

from Bio.Seq import Seq

import argparse

# Global, unchanging
nuc_tup = ("A", "C", "G", "T")

def parse_args():
    parser = argparse.ArgumentParser(prog = 'nuc_counter.py', conflict_handler = 'resolve')
    # parser.add_argument('-bl', type = str, required = True, help = '=> .txt with organism blacklist e.g. mm10')
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

def main():
    """Prints counts of A, C, G, T nucleotides in a DNA sequence string"""

    args = parse_args()
    input_seq = set_input(args)
    DNA_seq = Seq(input_seq)
    if len(DNA_seq)-1 > 1000:
        print("Size limit of 1000 exceeded. Exiting.")
        return(None)
    count_str = ""
    for ind, nuc in enumerate(nuc_tup):
        if ind != len(nuc_tup)-1:
            count_str += f"{DNA_seq.count(nuc)} "
        else:
            count_str += f"{DNA_seq.count(nuc)}"
    print(count_str)

if __name__ == '__main__':
    main()