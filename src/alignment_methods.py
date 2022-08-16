# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 17:06:27 2022

@author: pierre
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

@timer
def naive_exact_matching(a, b):
    """
    Align a sequence along a longer one using naive method (i.e no skip)
    - a is the longest sequence
    - return all occurences of exact matching as positions
    """
    assert sequence_is_valid(a), "Sequence format invalid"
    assert sequence_is_valid(b), "Sequence format invalid"
    matches_pos = []
    # on scroll b sur a
    for i in range(len(a)-len(b)+1):
        for j in range(len(b)):
            # tant que b correspond a a on continue
            if a[i+j] != b[j]:
                break
            # si on arrive au bout de b -> on a un match !
            if j == len(b)-1:
                matches_pos.append(i)
    return matches_pos
    
class GoodSuffixIndex:
    """
    - Used for Boyer-Moore good suffix rule
    
    - Create a list in which :
      - each item of index i is a dictionary of substring of lenght i where
          - key is a substring of length i
          - value is a list of the positions of the substring within the sequence
    
    - Has method "substring_pos" that gives the list of the positions
      for a specified substring
      
    - lmax cap the length of substring stored. 
      Set to -1 to generate all possible substring
    """
    def __init__(self, sequence):
        assert sequence_is_valid(sequence), "Sequence format invalid"
        self.sequence = sequence
        
        # initialise data structure
        self.sub_str_index = [{} for l in range(len(self.sequence)+1)]
        
        # scan the sequence and generate the data structure
        for start in (range(len(self.sequence)+1)):
            for stop in (range(len(self.sequence)+1)):
                if start < stop:
                    sub_str = self.sequence[start:stop]
                    if sub_str not in self.sub_str_index[len(sub_str)]:
                        self.sub_str_index[len(sub_str)][sub_str] = [start]
                    else:
                        self.sub_str_index[len(sub_str)][sub_str].append(start)
   
    def __str__(self):
        output = self.sequence+"\n"
        for i, substrings in enumerate(self.sub_str_index):
            output += "Substring size {}:\n".format(i)
            for k, v in substrings.items():
                output += "{} at {}\n".format(k, v)
        return output
    
    def substring_pos(self, target):
        """
        Returns
        -------
        LIST
            Positions of target substring within the sequence object.
        """
        try:
            return self.sub_str_index[len(target)][target]
        except:
            return []

class BadCarIndex:
    """
    - Used for Boyer-Moore bad character rule
    
    - Create a table-like dic in which :
      - key is the index considered in the sequence
      - value is a dic. of nucleotide(key) and its distance(val) from the index
    
    - Has method "dist_to_next" that gives the distance to next downstream 
      specified nucleotide
    """
    def __init__(self, sequence):
        assert sequence_is_valid(sequence), "Sequence format invalid"
        self.seq = sequence
        self.rev_seq = reverse(self.seq)
        self.bad_car_dist = {}
        for i in range(len(self.rev_seq)):
            if self.rev_seq[i] not in ["A", "C", "G", "T"]:
                sys.exit("ERROR. Not a valid sequence : Unknown nucleotide")
            self.bad_car_dist[i] = {}
            for j in range(i+1, len(self.rev_seq)):
                nxt_car = self.rev_seq[j]
                if nxt_car not in self.bad_car_dist[i].keys():
                    self.bad_car_dist[i][nxt_car] = int(j-i)
    
    def __str__(self):
        output = "BAD CAR INDEX\n"
        for k, v in self.bad_car_dist.items():
            output += "{} : {}\n".format(k, v)
        return output
    
    def dist_to_next(self, i, car):
        """
        Parameters
        ----------
        i : integer
            the mismatch's position in the read/short_sequence
        car : string of lenght 1
            the mismatched car in the ref_genome/lon_sequence

        Returns
        -------
        distance to the next "car" caracter in the read/short_sequence
        """
        assert i >= 0 and i <= len(self.rev_seq), "ERROR. position out of range"
        assert car in ["A", "C", "G", "T"], "ERROR. Not a valid nucleotide"
        # we work on the rev_seq so we invert the index
        i =  len(self.rev_seq) - 1 - i
        if car in self.bad_car_dist[i]:
            return self.bad_car_dist[i][car]
        else:
            return None           

class AlignRes:
    """
    Used to print the results of an alignment
    - a : seq la plus longue (genref)
    - b : seq qu'on cherche a trouver dans a (read)
    - match_pos : une liste de position (dans a evidemment) ou y'a match
    """
    def __init__(self, a, b, match_pos : list):
        assert sequence_is_valid(a), "Sequence format invalid"
        assert sequence_is_valid(b), "Sequence format invalid"
        self.a = a
        self.b = b
        self.pos = match_pos if match_pos else []
    def __str__(self):
        output = "ALIGNMENT RESULTS :\n"
        if not self.pos:
            output += "No results"
            return output
        for i, pos in enumerate(self.pos):
            assert pos > 0 , "Match position error (negative)"
            assert pos <= len(self.a) - len(self.b) + 1 , "Match position error (out of bound)"
            output += "# Alignment {} at position {} :\n".format(i+1, pos)
            output += self.a + "\n" + "-"*pos + self.b + "-"*(len(self.a)-len(self.b)-pos)+"\n"
        return output

@timer
def bad_car_exact_matching(a, b):
    """
    Use the bad caracter rule from boyer moore to skip some alignment check
    """
    assert sequence_is_valid(a), "Sequence format invalid"
    assert sequence_is_valid(b), "Sequence format invalid"
    match_positions = []
    bad_car = BadCarIndex(b)
    i = 0
    while True:
        if i > (len(a) - len(b)):
            break
        for j in reversed(range(len(b))):
            if a[i+j] != b[j]:
                to_skip = bad_car.dist_to_next(j, a[i+j])
                to_skip = to_skip if to_skip else j
                i += to_skip
                break
            if j == 0:
                match_positions.append(i)
                i += 1
    return match_positions

@timer
def boyer_moore(a, b):
    """
    Boyer Moore Exact Matching: skips useless alignment check by using
    - bad caracter rule
    - good suffix rule
    - a = "long" sequence (ref genome)
    - b = "short" sequence (read)
    """
    assert sequence_is_valid(a), "Sequence format invalid"
    assert sequence_is_valid(b), "Sequence format invalid"
    match_positions = []
    bad_car = BadCarIndex(b)
    good_suffix = GoodSuffixIndex(b)
    i = 0
    while True:
        if i > (len(a) - len(b)):
            break
        for j in reversed(range(len(b))):
            if a[i+j] != b[j]:
                # BAD CAR
                bc_skip = bad_car.dist_to_next(j, a[i+j])
                bc_skip = bc_skip if bc_skip else j+1
                
                # GOOD SUFFIX
                # - search if there is a substring located 5' of the mismatch
                #   that correspond to a substring of the 3' matching part 
                #   (a.k.a.the "good suffix") (the longest the better).
                # - skip is set to max and will stay if no prefix corresponding 
                #   to a suffix within the matching part
                gs_skip = len(b)                            
                
                # go through each possible substring of the matching part
                # scan the matching part from longest suffixes to shortest
                matching_part = b[j+1:]
                if matching_part != "":
                    for s in range(len(matching_part)):
                        gs_pos = good_suffix.substring_pos(matching_part[s:])
                        # check if a identical seq occurs in the part before the mismatch
                        # record the last one because its the closest from the split
                        for pos in gs_pos:
                            if pos <= j:
                                gs_skip = len(b) - len(matching_part[s:]) - pos
                        # 1st correspondance to occurs is the longest gs possible  
                        if gs_skip != len(b):
                            break     
                # if there is no matching part, i.e. mismatch immediatly on the
                # last caracter of b no gs rule we leave the job to the bc rule
                else:
                    gs_skip = 0
                i += max(bc_skip, gs_skip)
                break
            # If we arrive here : no mismatch for the whole of b
            if j == 0:    
                match_positions.append(i)
                i += 1
    return match_positions          

def main():
    PATH = "C:/Users/pierre/Documents/Jupyter_workspace/python_script/training_align/test.fasta"
    myfasta = read_fasta(PATH)
    # a = myfasta[">seq_3A"]
    # b = myfasta[">seq_3B"]
    
    # print("Naive exact matching")
    # print(AlignRes(a, b, naive_exact_matching(a, b)))

    # print("Bad caracter rule exact matching")
    # print(AlignRes(a, b, bad_car_exact_matching(a, b)))
    
    #a = myfasta[">seq_4A"]
    #b = myfasta[">seq_4B"]
    #print("Naive exact matching")
    #print(AlignRes(a, b, naive_exact_matching(a, b)))
    #print("Bad caracter rule exact matching")
    #print(AlignRes(a, b, bad_car_exact_matching(a, b)))
    
    a = myfasta[">seq_6A"]
    b = myfasta[">seq_6B"]
    #print(AlignRes(a, b, bad_car_exact_matching(a, b)))
    #print(a[8])
    #print(b)
    
    print("Boyer Moore Exact Matching")
    print(AlignRes(a, b, boyer_moore(a, b)))
    print(AlignRes(a, b, []))

if __name__ == "__main__":
    main()
    

    
    
    
# NEXT

# - implement test including 
#   - empty sequence/file
#   - invalid nucleotides
#   - match on first/last
#   - overlapping matches

# - implement indexing (kmer...) for Offline algo

# - implement approximate matching with edit distance

# - implement approximate matching with pigeon hole

# - a gui with kivy ?
    
    