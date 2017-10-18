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


def find_motifs(fasta_path, output_path, mod='zoops', alphabet='rna', minw=6, maxw=25, options=None):
    """find motifs from fasta file"""

    cmd = 'meme '+fasta_path+' -'+alphabet+' -mod '+mod+' -nmotifs '+str(nmotifs)+' -minw '+str(minw)+' -maxw '+str(maxw)+' -oc '+output_path+' '+options
    if verbose:
        print('>>' + cmd)
    os.system(cmd)
