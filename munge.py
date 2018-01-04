#!/bin/python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import pandas as pd
import numpy as np
from . import meme

def convert_one_hot(sequence, max_length=None):
    """convert DNA/RNA sequences to a one-hot representation"""

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

        # handle boundary conditions with zero-padding
        if max_length:
            offset1 = int((max_length - seq_length)/2)
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
    """generate background sequences"""

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



def split_dataset(one_hot, labels, valid_frac=0.1, test_frac=0.2):
	"""split dataset into training, cross-validation, and test set"""

	def split_index(num_data, valid_frac, test_frac):
		# split training, cross-validation, and test sets

		train_frac = 1 - valid_frac - test_frac
		cum_index = np.array(np.cumsum([0, train_frac, valid_frac, test_frac])*num_data).astype(int)
		shuffle = np.random.permutation(num_data)
		train_index = shuffle[cum_index[0]:cum_index[1]]
		valid_index = shuffle[cum_index[1]:cum_index[2]]
		test_index = shuffle[cum_index[2]:cum_index[3]]

		return train_index, valid_index, test_index


	# split training, cross-validation, and test sets
	num_data = len(one_hot)
	train_index, valid_index, test_index = split_index(num_data, valid_frac, test_frac)

	# split dataset
	train = (one_hot[train_index], labels[train_index,:])
	valid = (one_hot[valid_index], labels[valid_index,:])
	test = (one_hot[test_index], labels[test_index,:])
	indices = [train_index, valid_index, test_index]

	return train, valid, test, indices



def save_dataset_hdf5(save_path, train, valid, test):
	"""save datasets as hdf5 file"""
	with h5py.File(save_path, "w") as f:
		dset = f.create_dataset("X_train", data=train[0].astype(np.float32), compression="gzip")
		dset = f.create_dataset("Y_train", data=train[1].astype(np.float32), compression="gzip")
		dset = f.create_dataset("X_valid", data=valid[0].astype(np.float32), compression="gzip")
		dset = f.create_dataset("Y_valid", data=valid[1].astype(np.float32), compression="gzip")
		dset = f.create_dataset("X_test", data=test[0].astype(np.float32), compression="gzip")
		dset = f.create_dataset("Y_test", data=test[1].astype(np.float32), compression="gzip")


def split_rnacompete_dataset(data, targets, experiment, valid_frac):
    """split RNA compete dataset into set A sequences for training
        and validation and set B sequences for testing."""

	index = np.where(experiment == 'A')[0]
	num_seq = len(index)
	num_valid = int(num_seq*valid_frac)

	shuffle = np.random.permutation(num_seq)

	X_train = data[shuffle[num_valid:]]
	Y_train = targets[[shuffle[num_valid:]]]
	train = (X_train, Y_train)

	X_valid = data[shuffle[:num_valid]]
	Y_valid = targets[[shuffle[:num_valid]]]
	valid = (X_valid, Y_valid)

	test_index = np.where(experiment == 'B')[0]
	X_test = data[test_index]
	Y_test = targets[test_index]
	test = (X_test, Y_test)

	return train, valid, test
