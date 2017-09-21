#!/bin/python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import pandas as pd
import numpy as np


def generate_fasta(sequences, fasta_path):
    """generate fasta file from an array of sequences
    """

    with open(fasta_path, 'w+') as f:
        for i in xrange(len(sequences)):
            f.write('>seq '+str(i))
            f.write('\n')
            f.write(sequences[i])
            f.write('\n')


def parse_sequences(seq_path):
    """Parse fasta file for sequences"""

    # parse sequence and chromosome from fasta file
    num_data = np.round(sum(1 for line in open(seq_path))/2).astype(int)
    fin = open(seq_path, "r")
    sequences = []
    for j in range(num_data):
        coord = fin.readline()
        line = fin.readline()[:-1].upper()
        sequences.append(line)
    sequences = np.array(sequences)
    return sequences


def count_fasta_entries(file_path, factor=2):
    with open(file_path, 'r') as f:
        counts = 0
        for line in f:
            counts += 1
    return counts/factor
