#!/usr/bin/env python3

"""ROSALIND GenBank intro - bioinformatics armory

This script returns the number of nucleotide records PUBLISHED for a genus between a start and end date.

The input date format is YYYY/M/D.

Please keep in mind NCBI's request requirements: 
https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
"""

from Bio import Entrez

import argparse

#help(Entrez.esearch)

def parse_args():
    parser = argparse.ArgumentParser(prog = 'gbk_counter.py', conflict_handler = 'resolve')
    parser.add_argument('-db', type = str, required = True, help = '=> db type: choose "nucleotide" or "pubmed"')
    parser.add_argument('-email', type = str, required = True, help = '=> email for Entrez queries')
    parser.add_argument('-genus', type = str, required = True, help = '=> genus of interest')
    parser.add_argument('-start', type = str, required = True, help = '=> start date to look for pubs')
    parser.add_argument('-end', type = str, required = True, help = '=> end date to look for pubs')
    return(parser.parse_args())

def count_entrez_pubs(db_type: str, email: str, genus: str, start: str, end: str) -> int:
    """Searches either GenBank "nucleotide" or Pubmed "pubmed" entries for publication counts within a date range"""
    Entrez.email = email
    handle = Entrez.esearch(db=db_type, term=f"{genus}[Organism] AND {start}[PDAT]:{end}[PDAT]")
    record = Entrez.read(handle)
    return(record["Count"])

def main():
    # organize inputs from command line
    args = parse_args()

    # print counts
    print(count_entrez_pubs(args.db, args.email, args.genus, args.start, args.end))

if __name__ == '__main__':
    main()