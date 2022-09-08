# shortest_fa.py
Print FASTA record with shortest FASTA sequence length from a list of GenBank IDs.

# Usage
*   Print shortest fasta (-rettype) from a list of accession numbers (-id, -idfile) for the GenBank "nucleotide" db (-db)

    *   Setting -id
    
        ```python shortest_fa.py -email bvtsu@ucsd.edu -ids FJ817486 JX069768 JX469983```

    *   Setting -idfile
        
        ```python shortest_fa.py -email bvtsu@ucsd.edu -idfile 'tests/rosalind_frmt.txt'```