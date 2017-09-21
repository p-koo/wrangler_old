#!/bin/python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import pandas as pd
import numpy as np
from . import meme



def convert_one_hot(sequence, max_length):
    one_hot_seq = []
    for seq in sequence:
        seq = seq.upper()
        seq_length = len(seq)        
        one_hot = np.zeros((4,seq_length))
        index = [j for j in xrange(seq_length) if seq[j] == 'A']
        one_hot[0,index] = 1
        index = [j for j in xrange(seq_length) if seq[j] == 'C']
        one_hot[1,index] = 1
        index = [j for j in xrange(seq_length) if seq[j] == 'G']
        one_hot[2,index] = 1
        index = [j for j in xrange(seq_length) if (seq[j] == 'U') | (seq[j] == 'T')]
        one_hot[3,index] = 1
        offset1 = ((max_length - seq_length)/2).astype(int)
        offset2 = max_length - seq_length - offset1

        if offset1:
            one_hot = np.hstack([np.zeros((4,offset1)), one_hot])
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((4,offset2))])

        one_hot_seq.append(one_hot)

    # convert to numpy array
    one_hot_seq = np.array(one_hot_seq)

    return one_hot_seq


def filter_nonsense_sequences(sequences):
    """Parse fasta file for sequences"""

    # parse sequence and chromosome from fasta file
    good_index = [] 
    filter_sequences = []
    for i, seq in enumerate(sequences):
        if 'N' not in seq.upper():
            good_index.append(i)
            filter_sequences.append(seq)
    return np.array(filter_sequences), np.array(good_index)
    


def process_background_sequences(fasta_path, window=None, background='dinuc', verbose=0):

    # generate different backgrounds
    #backgrounds = ['genome', 'dinuc', 'random']

    if background == 'nonspecific':
        # background sequences derived from genomic sequences
        background_path = fasta_path

    elif background == 'dinuc':
        # background sequences derived from dinucleotide shuffle
        index = fasta_path.index('.')
        background_path = fasta_path[:index] + '_dinuc.fa'
        meme.shuffle(fasta_path, background_path, kmer=2, verbose=verbose)

    elif background == 'random':
        # background sequences derived from random shuffle
        background_path = fasta_path[:index] + '_random.fa'
        meme.shuffle(fasta_path, background_path, kmer=1, verbose=verbose)
        
    # process background sequences
    neg_sequences = fasta.parse_sequences(background_path)

    return neg_sequence

