# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 17:06:27 2022

@author: pierre
"""

import sys
from tools import timer, read_fasta, sequence_is_valid, reverse 
from tools import split_seq, AlignedPair



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
def boyer_moore(a, b, verbose = False):
    """
    Boyer Moore Exact Matching: skips useless alignment check by using
    - bad caracter rule
    - good suffix rule
    - a = "long" sequence (ref genome)
    - b = "short" sequence (read)
    - returns a list of numbers : the starting indexes of matches within "a"  
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
 
def inexact_matching(ref_gen, read_sq, mismatch_nbr = 1, verbose = False):
    """
    Inexact matching by applying pigeon hole principle to Boyer Moore algorithm 
    - ref_gen  = "long" sequence (reference genome)
    - read_sq = "short" sequence (read)
    - mismatch_nbr = number of allowed mismatch
    """
    assert sequence_is_valid(ref_gen), "Sequence format invalid"
    assert sequence_is_valid(read_sq), "Sequence format invalid"
    match_positions  = set()
    match_candidates = set()
    # pigeon hole : for inexact matching. Split the read into n parts
    # where n = mismatch_nbr + 1. eg. for 2 mismatches allowed split in 3 parts
    for split_pos, part in  split_seq(read_sq, mismatch_nbr):
        # find an exact match on a part of ref_sq -> its a possible inexact match
        for matching_part_pos in boyer_moore(ref_gen, part):
            # record the starting point of the whole read within the ref_gen 
            match_candidates.add(matching_part_pos - split_pos)
    # check whether the edit distance for the complete read is still within bounds
    for pos in match_candidates:
        res = AlignedPair(ref_gen, read_sq, pos, trim = 1)
        if res.get_edist() <= mismatch_nbr:
            match_positions.add(pos)
    if verbose:
        print("{} matches found at positions {}\n".format(len(match_positions), sorted(match_positions)))
    return sorted(match_positions)
        
def do_tests():
    PATH = "/home/pierre/git/Alignements_methods/test/test.fasta"
    myfasta = read_fasta(PATH)
    a = myfasta[">seq_7A"]
    b = myfasta[">seq_7B"]
    
    #print(myfasta)
    print("sequences are :\n{}\n{}".format(a, b))
    #res = inexact_matching(a, b, mismatch_nbr = 2, verbose = True)
    #print("res = ", res)
    #for pos in res:
    #    print("Match at position {}:\n{}".format(pos, AlignedPair(a, b, pos)))
    

def main():
    do_tests()

if __name__ == "__main__":
    main()
       
###############    
# TO DO NEXT ->
###############

# - push to git
    
# - check for and comment out unused tools from tools.py

# - implement inexact matching with gaps

# - implement indexing

# - gerer le @timer avec verbose
    
# - gerer la verbosite
    
# - gerer l'impression de resultats/alignements

# - gerer l'edist quand les read debordent du gen_ref

# - check what to do when a match position results in two seq without overlap

# - implement trim = 2 dans class AlignedPair

# - implement test including 
#   - empty sequence/file
#   - invalid nucleotides
#   - match on first/last
#   - overlapping matches
#   - different lengths seq
#   - sticking out alignment/matches
    
# - implement "repeat" compression format for DNA (such as A5T3CG12)
#   use it to store and search for substring
#   store those substring in alphabetical order
    
# - binary compression with two "overlaping" 16 bit 0101 0101 0101 0001
#   - first define if is it a purine or a pyrimidine
#   - second define if it is default pur/pyr (defautlt pur = A default pyr = T 
#     for ex while G and C are 1 pur/pyr)
#   - A 16 bit encoded car will take 16 bit space (duh) -> 2 nucleotides encoded on 32 bits
#   - With this system 32 bits bit will encode for 16 nucleotides (so 8 times compression)
#   - but in case of long repeats coding for them as two car (e.g. A9) so 32 bits
#     lead to better compression.
#   - a first (third) byte could be used for storing repeated position

#   - with 16 bit = a byte
#     in case of repeated sequence
#         REPEAT       BINARY
# 1 nuc   1 byte       2 byte 
# 9 nuc   2 byte       2 byte
# 16 nuc  3 byte       2 byte
# 99 nuc  3 byte       28 bytes
#     the longest the repeat the worst performs binary compression compared to
#     repeat

# - implement indexing (kmer...) for Offline algo

# - implement approximate matching with edit distance

# - a gui with kivy ? 

# - check existence of interactions GODOT-PYTHON for gui

    
##############
# DONE \(^_^)/
##############
    
# DONE : - read fasta bug si ligne vide dans le fichier (a la fin aussi)
# DONE : - implement approximate matching with pigeon hole
# DONE : - implement different lenght edit distance computation
# DONE : - fix the import
# DONE : gerer quand les read debordent dans inexact matching
# DONE : - take the tools from tools.py (reverse, split seq...)    
# DONE : - implement same lenght edit distance computation
# DONE : - fix the AlignRes class so that reads can stick out from the ref genome   
# DONE : - gerer l'affichage des align quand les read debordent
