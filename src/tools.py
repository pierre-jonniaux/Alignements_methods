#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:54:21 2022

@author: pierre

Misceleanous Tools for
- alignement_methods.py
- dna_pack.py
"""




import sys
from time import time

def timer(func):
    """
    A function timer decorator

    """
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        stop = time()
        print("{} executed in {:.4f} s\n".format(func.__name__, (stop-start)))
        return result
    return wrapper

def read_fasta(file_path):
    """
    - Open and read a fasta formated file (!ONE LINE PER SEQUENCE!)
    - Return a dictionary where k,v are seq. names and seq.data respectively
    """
    sequences = {}
    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if line[0] != ">":
                sys.exit("Not a fasta file")
            seq_name = line
            line = f.readline().rstrip()
            sequences[seq_name] = line
    return sequences    
        
def print_fasta(fasta_sequences):
    for name, seq in fasta_sequences.items():
        print(name, "\n", seq[:len(name)], sep="")

def sequence_is_valid(sequence):
    assert len(sequence) > 0, "Sequence of lenght 0" 
    for car in sequence: 
        if car not in ["A", "C", "G", "T"]:
            return False
    return True   

def reverse(string):
    """
    Return a reversed string
    """
    string = "".join(reversed(str(string)))
    return string

def complement(seq):
    """
    Take a nucleotide sequence and returns the complement sequence.
    """
    assert sequence_is_valid(seq), "Sequence format invalid"
    comp =  ""
    nuc_comp = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}
    for nuc in seq:
        comp += nuc_comp[nuc]
    return comp
