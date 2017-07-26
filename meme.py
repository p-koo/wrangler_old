from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys


def shuffle(fasta_path, output_path, kmer=2, verbose=1):
    """perform dinucleotide shuffle or random shuffle (kmer=1) on sequences in fasta file"""
    cmd = 'fasta-shuffle-letters -kmer '+str(kmer)+' -line 2000 '+fasta_path+' '+output_path 

    if verbose:
        print('>>' + cmd)
    os.system(cmd)    


