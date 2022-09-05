# nuc_counter.py
Print space-separated counts of nucleotides in a given string or text file containing a string of DNA.

# Usage
* DNA sequence - string input

    ```python nuc_counter.py -seq 'ATCG'```

* DNA sequence - file containing a string input

    ```python nuc_counter.py -seqfile 'tests/rosalind_ini.txt'```

* DNA sequence - string and file inputs -> <b>Override to use file only</b>

    ```python nuc_counter.py -seq 'ATCG' -seqfile 'tests/rosalind_ini.txt'```