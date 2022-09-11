#!/usr/bin/env python3

"""ROSALIND New motif discovery - bioinformatics armory

This script returns a regular expression associated with a file of provided PROTEIN FASTAs.

The input is a set of protein strings in FASTA format that share a motif with minimum length 20.
"""

# from Bio import motifs
# help(meme) # doesn't give us regex

import argparse
import subprocess
import xml.etree.ElementTree as ET

def parse_args():
    parser = argparse.ArgumentParser(prog = 'MEME_regex.py', conflict_handler = 'resolve')
    parser.add_argument('-i', type = str, required = True, help = '=> path/to/fasta file')
    return(parser.parse_args())

def get_prot_meme(input_fa: str):
    """Takes a path to an input fasta and runs MEME tool via shell
    
    Minimum width of motif is 20

    Returns a statement about the error state
    """

    shCommand = f"meme {input_fa} -protein -minw 20"
    process = subprocess.Popen(shCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return(error)

def get_meme_regex() -> dict[str, str]:
    """Parses a meme.xml output to identify motifs>motif>regex output
    
    Returns a {motif_name: regex} dict
    """

    etree = ET.parse('meme_out/meme.xml')
    eroot = etree.getroot()
    regex_dict = {}
    for motifs in eroot.findall('motifs'):
        for motif in motifs.findall('motif'):
            motif_name = motif.attrib['id']
            regex_dict[motif_name] = motif.find('regular_expression').text.strip()
    return(regex_dict)


def main():
    args = parse_args()
    stderr = get_prot_meme(args.i)
    if stderr == None:
        print(get_meme_regex())
    else:
        print(f"MEME Suite error: {stderr}")

if __name__ == '__main__':
    main()