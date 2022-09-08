#!/usr/bin/env python3

"""ROSALIND Data Formats - bioinformatics armory

This script returns the shortest of the strings associated with the IDs in FASTA format.

The input is a collection of space-separated n (n ≤ 10) GenBank entry IDs.

Please keep in mind NCBI's request requirements: 
https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
"""

from Bio import Entrez, SeqIO

import argparse

#help(Entrez.esearch)

def parse_args():
    parser = argparse.ArgumentParser(prog = 'shortest_fa.py', conflict_handler = 'resolve')
    parser.add_argument('-email', type = str, required = True, help = '=> email for Entrez queries')
    parser.add_argument('-ids', required = False, nargs='+', default=[], help = '=> space-separated GenBank IDs')
    parser.add_argument('-idfile', type = str, required = False, help = '=> path/to/GenBank IDs file')
    return(parser.parse_args())

def set_input(parsed_args) -> list:
    """If at least one arg exists, return list of IDs.

    Note: if both exist, override with idfile list of IDs.
    """
    
    if parsed_args.ids is None and parsed_args.idfile is None:
        return(None)
    if parsed_args.idfile is not None:
        with open(parsed_args.idfile, 'r') as idfile:
            return(idfile.read().replace('\n','').split())
    else:
        return(parsed_args.ids)

def get_entrez_fa(email: str, id_list: list[str]) -> list:
    """Searches GenBank "nucleotide for ID-associated FASTAs"""

    Entrez.email = email
    print(id_list)
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    return(records)

def shortest_seq(record_list: list) -> tuple[str, str]:
    """Enumerate through records for shortest FASTA sequence length"""

    # initiate len comparison vars
    desc_tracker = ""
    seq_tracker = "" 
    len_tracker = 0

    # initialize with first ind values
    for ind, rec in enumerate(record_list):
        if ind == 0:
            desc_tracker = rec.description
            seq_tracker = rec.seq
            len_tracker = len(rec.seq)
        
        # only keep those that are shorter than the prev shortest
        else:
            if len(rec.seq) < len_tracker:
                desc_tracker = rec.description
                seq_tracker = rec.seq
                len_tracker = len(rec.seq)
    return(f">{desc_tracker}\n{seq_tracker}")

def main():
    # organize inputs from command line
    args = parse_args()
    input_ids = set_input(args)

    # impose limits
    if len(input_ids) > 10:
        print("Size limit of 10 exceeded. Exiting.")
        return(None)

    # print counts
    GenBank_records = (get_entrez_fa(args.email, input_ids))
    print(shortest_seq(GenBank_records))


if __name__ == '__main__':
    main()