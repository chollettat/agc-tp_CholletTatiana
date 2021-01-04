#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw

__author__ = "Tatiana Chollet"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Tatiana Chollet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Tatiana Chollet"
__email__ = "chollettat@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()



## 1- Dé-duplication en séquence “complète”
def read_fasta(amplicon_file, minseqlen):
    """
    Read the file
    :Parameters:
        amplicon_file: file fasta.gz
        minseqlen: minimum length of sequences
    Returns: sequences generator
    """
    file = gzip.open(amplicon_file, "rt") #gzip because of gz extension
    sequence = ""
    for line in file:
        if line.startswith(">"):
            if len(sequence)>= minseqlen:
                yield sequence
            sequence = ""
        else:
            sequence = sequence + line.strip() #remove whitespace
    yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Create a dictionary containing the unique amplicons
    :Parameters:
        amplicon_file: file fasta.gz
        minseqlen: minimum length of sequences
        mincount: minimum sequence count
    Returns: generator of unique sequences and their occurrence
    """
    dict_sequence = {} #We use {} to use a dictionnary
    list_sequence = read_fasta(amplicon_file, minseqlen)

    for sequence in list_sequence:
        if sequence in dict_sequence: #if the sequence is in the dico, then we count +1 occurence
            dict_sequence[sequence] += 1
        else :
            dict_sequence[sequence] = 1
    for seq, occ in sorted(dict_sequence.items(), key=lambda item: item[1], reverse = True): #item: item[1] sorted on occurence in dict_sequence, reverse = True return in descending order
        if occ >= mincount:
            yield [seq, occ]



## 2- Recherche de séquences chimériques par approche “de novo”
def get_chunks(sequence, chunk_size):
    """
    ----------------------
    :Parameters:
        sequence: sequence
        chunk_size: size of sequence
    Returns: list of sub-sequences of size l not overlapping, at least 4 segments must be obtained per sequence
    """
    chunks = []
    for i in (range(0, len(sequence), chunk_size)):
        if i+chunk_size<=len(sequence):
            chunks.append(sequence[i:i+chunk_size])
    if len(chunks)>=4: #4 segments per sequence
        return chunks

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    """
    ----------------------
    :Parameters:
        sequence: sequence
        kmer_size: size of kmer
    Returns: k-mer generator
    """
    for comp in range ( (len(sequence)+1) - kmer_size):
        kmer = sequence[comp : comp + kmer_size]
        yield kmer


def get_identity(alignment_list):
    """
    Calculate the percentage of identity between the two sequences
    :Parameters:
        alignment_list: alignment
    Returns: percentage of identify with 2 digits after the decimal point
    """
    count = 0
    len_list = len(alignment_list[0])
    for i in range(len_list):
        if alignment_list[0][i] == alignment_list[1][i]:
            count += 1
    id = count/len_list * 100
    return id


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

## 3 - Regroupement glouton
def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    file = open(output_file, "w")
    index = 1
    for OTU in enumerate(OTU_list):
        file.write(">OTU_"+str(index)+" occurence:"+str(OTU[1])+"\n")
        file.write(fill(str(OTU[0])+"\n"))
        index += 1


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()
