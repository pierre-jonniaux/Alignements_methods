#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 17:08:56 2022

@author: pierre
"""

#import sys
#from bitstring import BitArray
import bitstring

from tools import sequence_is_valid

# FINALEMENT
# 1 - cut sequence data in 8 nucleotides substrings
# 2 - convert those 8 nuc into 2 8-bit sequence
# 3 - convert and store it as caracter 


def _byte_conv(seq):
    """
    yield 2 bytes/octet for a 8*n length nucleotide sequence
    OR
    return 2 (3) lists where sequence data is compressed and stored in byte format
    one byte = 8 sequence nucleotide, one per bit 
    - pur vs pyr : bit 0 means a purine, 1 a pyrimidine
    - A vs G or C vs T : bit 0 means either (A/G) and 1 means (C/T)
    (- location of repeats)
    
    """
    assert len(seq)%8 == 0, "Sequence length not a multiple of 8"
    assert sequence_is_valid(seq), "Sequence invalid"
    
    base = "".join(["0" if nuc in ["A", "G"] else "1" for nuc in seq])
    sort = "".join(["0" if nuc in ["A", "C"] else "1" for nuc in seq])
    base = bitstring.BitArray(bin=base).tobytes()
    sort = bitstring.BitArray(bin=sort).tobytes()
    
    print(type(base))
    return (base, sort)

print(_byte_conv("ACGTACGT"))

def _dna_seq(byte):
    #assert des trucs 
    pass


def bin_pack(sequence):
    """
    Use DNA to byte conversion to store nucleotides sequence in binary file
    for low space storage purpose and fast susbtring search
    
    Format:
    - Only A,C,G,T nucleotides are accepted
    - Only one sequence a time
    - raw sequence or fasta (in which case header is ignored)
    
    
    pros :
    - compact as 32 bits store 16 nuc (not contemplating repeat data)
    - manipulation of binary data , substring binary search using trees
    - 2(3) kinds of data : if parrallel search can be implemented, then
      searching simultaneously into the different data can allows
      each search to be hastened (chunks of concurent trees are removed)
    - could use boolean operations on bytes to compare/search substring
      during alignement
    cons :
    - no coding for unknown nucleotides (N)
    - unprintable (use of the first of the 8 bits, control caracters...)
    """
    
    # TODO
    # - load sequence
    # - cut it into 8nuc pieces
    # - send to _byte_conv
    # - see how to save the output of _byte_conv, of class "bytes", to a file
    pass

def bin_unpack(file):
    """
    reverse of the bin_pack fucntion
    """
    # TODO
    # - will take a binary file as input
    # - will use yield to cut into bytes
    # - will use dnaseq_conv
    pass

def r_pack():
    """
    return a single list where sequence data is compressed and stored in byte format
    - 4 bits per at least one nucleotide 3 at most
    - 2 bit indicate if the nucleotid is repeated (from 0 to 3 times)
    - 2 bit code for the nucleotide (00 : A, 0I : C, I0 : G, II : T fro example) 
    - a 16 bit byte potentialy codes for 4nuc at min and 12 at most
    """
    pass

seq = "ATCTAGCC"
#print("size of string seq : {}".format(sys.getsizeof(seq)))
bitseq = ''.join(format(i, '08b') for i in bytearray(seq, encoding ='utf-8'))
#print(bitseq)



        
        
#print(d1,"\n",d2, sep="")

#print(int(d1,2),"\n",int(d2,2), sep="")

#print("size after {}".format(sys.getsizeof(int(d1))+sys.getsizeof(int(d2))))



#seq = "A"
#print("size of bit seq 1 and 2 :\n")
#print(d1)
#print(d2)
#
#print(sys.getsizeof(bitstring.BitArray(bin=d1).tobytes()))
#print(sys.getsizeof(d2))

#bitseq = ''.join(format(i, '08b') for i in bytearray(seq, encoding ='utf-8'))
  
#print(bitseq)

#v1 = bitstring.BitArray(bin=bitseq)

#my_bytes = v1.tobytes()

#print(sys.getsizeof(my_bytes))


print("#################################")
#string = "s"
##string = seq
#print("String : ", string)
#bitseq = ''.join(format(i, '08b') for i in bytearray(string, encoding ='utf-8'))
#print("ByteString : ", bitseq)
#my_bytes = bitstring.BitArray(bin=bitseq).tobytes()
#print("Mybytes : ", my_bytes)
#print("Size in bytes : {}".format(sys.getsizeof(my_bytes)))
#print("Size in bytes : {}".format(sys.getsizeof(string)))
#print("len", len(my_bytes))
#print("\n")
#
#
#print("Pour :{}\n".format(string))
#print("taille utf 8 :", len(string.encode('utf-8')))
#print("taille utf16 :", len(string.encode('utf-16')))
#print("taille ASCII :", len(string.encode('ASCII')))
#print("\n")
#
#string = d1
#
#print("Pour :{}\n".format(string))
#string = bitstring.BitArray(bin=string).tobytes()
#print("taille utf 8 :", sys.getsizeof(string))
#print(string)
#print(type(string))














