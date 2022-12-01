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
        print("{} executed in {:.4f} s".format(func.__name__, (stop-start)))
        return result
    return wrapper

def read_fasta(file_path):
    """
    - Open and read a fasta formated file (!ONE LINE PER SEQUENCE!)
    - Return a dictionary where k,v are seq. names and seq. data respectively
    """
    sequences = {}
    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip() 
            if line != "":
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
        if car not in ["A", "C", "G", "T", "N", "-"]:
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

    
def split_seq(seq, n = 1):
    """
    Split a sequence in parts of equal lenght. 
    Last part gets the rest (is longer) if split is uneven.
    - seq = the sequence to split
    - n = the number of split (will return n+1 parts)
    Returns a list of tuples : (string, integer) where 
      - the string is the sequence part
      - the integer is the part's starting position relative to the sequence  
    """
    part_len = round(len(seq)/(n + 1))
    return [(i, seq[i:i+part_len]) for i in range(0, len(seq), part_len)]   

def adjust_seq(refr_seq, read_seq, pos):
    """
    Take two sequences, "align" them according to the starting position of the
    second within the first. Where "align" means filling indels "-" where needed.
    
    refr_seq : first seq, the longer one, usually a reference genome part
    read_seq : second seq, the shorter one, usually a sequenced read
    pos : starting point of the read_seq within the read_seq (can be negative)
    
    Return a tuple with the two edited sequences
    
    Ex: 
    normalise_seq(ACGTACGT, TA, 3)
    returns
    ("ACGTACGT", "---TA---")
    """
    assert sequence_is_valid(refr_seq), "Sequence format invalid"
    assert sequence_is_valid(read_seq), "Sequence format invalid"
    assert len(refr_seq) >= len(read_seq), "Target sequence longer than the reference sequence"
    seq1 = ""
    seq2 = ""
    if pos < 0:
        seq1 += "-"*abs(pos) + refr_seq
        seq2 += read_seq + "-"*(len(refr_seq)-len(read_seq)-pos)
    elif pos + len(read_seq) > len(refr_seq):
        seq1 += refr_seq + "-"*(pos+len(read_seq)-len(refr_seq))
        seq2 += "-"*(pos) + read_seq
    else:
        seq1 += refr_seq
        seq2 += "-"*pos + read_seq + "-"*(len(refr_seq)-len(read_seq)-pos)
    return (seq1, seq2)

def edist_same_len(seq1 , seq2):
    """
    Compute the edit distance between two sequence of the same length.
    - seq1 = seq of length n
    - seq2 = seq of length n
    Returns an integer
    """
    assert len(seq1) == len(seq2), "Sequence are of different lengths" 
    edist = [False for i in range(len(seq1)) if seq1[i] != seq2[i]]
    return len(edist)


class AlignedPair:
    """
    Used to store two aligned DNA sequence of different length with printing 
    and edit distance methods. Inserts "-" in non-overllaping parts. Can trim 
    the sequences in various fashion.
    
    refr_seq : first seq, the longer one, usually a reference genome part
    read_seq : second seq, the shorter one, usually a sequenced read
    pos : starting point of the read_seq within the read_seq (can be negative)
    
    Exemple :
    Given "ACGTACGT" and "TA", possible alignments are at -1, 3 7
    -ACGTACGT ACGTACGT ACGTACGT-
    TA------- ---TA--- -------TA
        
    possible triming and associated edit distances are
    0 : no_trim
    -ACGTACGT ACGTACGT ACGTACGT-
    TA------- ---TA--- -------TA
    with edist 8, 6, 8 respectively

    1 : target_trim
    -A TA T-
    TA TA TA
    with edist 1, 0, 1 respectively
    
    2 : big_trim
    A TA T
    A TA T
    with edist 0 in all cases

    """
    def __init__(self, refr_seq, read_seq, pos, trim = 0):
        assert sequence_is_valid(refr_seq), "Sequence format invalid"
        assert sequence_is_valid(read_seq), "Sequence format invalid"
        assert len(refr_seq) >= len(read_seq), "Target sequence longer than the reference sequence"
        self.seq1 = ""
        self.seq2 = ""
        self.trim = trim
        if pos < 0:
            self.seq1 += "-"*abs(pos) + refr_seq
            self.seq2 += read_seq + "-"*(len(refr_seq)-len(read_seq)-pos)
        elif pos + len(read_seq) > len(refr_seq):
            self.seq1 += refr_seq + "-"*(pos+len(read_seq)-len(refr_seq))
            self.seq2 += "-"*(pos) + read_seq
        else:
            self.seq1 += refr_seq
            self.seq2 += "-"*pos + read_seq + "-"*(len(refr_seq)-len(read_seq)-pos)
        #return (self.seq1, self.seq2)
        if trim == 1:
            start = pos if pos >= 0 else 0
            self.seq1 = self.seq1[start:start+len(read_seq)]
            self.seq2 = read_seq
        if trim == 2:
            sys.exit("Trim 2 not implemented yet...")
         
    def __str__(self):
        return self.seq1 + "\n" + self.seq2    
        
    def get_edist(self):
        edist = 0
        for i, nuc1 in enumerate(self.seq1):
            if nuc1 != self.seq2[i]:
                edist += 1
        return edist

class AlignRes:
    """
    Used to store printable alignment results with their edition distance
    - a : seq la plus longue (reference sequence/genome)
    - b : seq qu'on cherche a trouver dans a (sequenced read)
    - match_pos : une liste de position (dans a evidemment) ou y'a match
    """
    def __init__(self, a, b, match_pos : list):
        assert sequence_is_valid(a), "Sequence format invalid"
        assert sequence_is_valid(b), "Sequence format invalid"
        self.a = a
        self.b = b
        self.positions = match_pos if match_pos else []
        # A dic where keys are the match positions relative to a/genref and
        # values are tuples containing 2 sequences filled with indels "-"
        self.align = {}
        for pos in self.positions:
            protuding_left = False
            protuding_right = False
            if pos < 0:
                protuding_left = True
            if pos + len(self.b) > len(self.a):
                protuding_right = True
            if protuding_left and protuding_right:
                sys.exit("Target sequence longer than the reference sequence")
            refr_seq = ""
            read_seq = ""
            if protuding_left:
                refr_seq += "-"*abs(pos) + self.a
                read_seq += self.b + "-"*(len(self.a)-len(self.b)-pos)
            elif protuding_right:
                refr_seq += self.a + "-"*(pos+len(self.b)-len(self.a))
                read_seq += "-"*(pos) + self.b
            else:
                refr_seq += self.a
                read_seq += "-"*pos + self.b + "-"*(len(self.a)-len(self.b)-pos)
            self.align[pos] = (refr_seq, read_seq)
            
                
    def __str__(self):
        output = "ALIGNMENT RESULTS :\n"
        if not self.positions:
            output += "No results"
            return output
        for refr_seq, read_seq in self.align.values():
            output += refr_seq + "\n" + read_seq + "\n"
        return output    
    
    def edist(self, pos):
        assert pos in self.align.keys(), "No alignement at this position"
        score = 0
        seq1, seq2 = self.align[pos]
        for i, nuc1 in enumerate(seq1):
            if nuc1 == seq2[i]:
                score += 1
        return score
    
#a = "ACGTACGT"
#b = "TA"
#pos = [-1, 3, 7]
#
#aaa = AlignRes(a, b, pos)
#
##print(aaa)
#
#res =  aaa.edist(-1)

#print("la distance pour alignement a {} est {}".format(3, res))

#edist
    
#    def split_seq(seq, n = 1):
#    """
#    Return indexes indicating positions where to cut in order to split
#    a sequence in parts of equal lenght. 
#    Note: Last part is of different length if split is uneven.
#    - seq = the sequence to split
#    - n = the number of split (will return n+1 part indexes)
#    Returns a list integer : the part's starting position relative to the sequence 
#    """
#    part_len = round(len(seq)/(n + 1))
#    return [i for i in range(0, len(seq), part_len)]
#    
    
    
    