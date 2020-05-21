# coding: utf-8

'''Utility functions
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

def translate(seq): 
    '''Original source code from https://www.geeksforgeeks.org/dna-protein-python-3/
    Translates DNA to amino acid based on standard codon table
    '''

    # If the input sequence is None, then return nothing
    if seq is None:
        return None

    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein = '' 
    if len(seq) % 3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            if codon in table.keys():
                protein += table[codon] 
            else:
                protein += '?'
    # If the sequence is not a multiple of 3, then just return a wildcard '*' character
    else:
        protein = '*'

    return protein
