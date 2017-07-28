#!/bin/python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import pandas as pd
import numpy as np
from . import bedtools, meme


def make_directory(base_path, dir_name):
    dir_path = os.path.join(base_path, dir_name)
    if not os.path.isdir(dir_path):
        print('making directory: '+dir_path)
        os.mkdir(dir_path)
    return dir_path
    

def parse_fasta_sequences(seq_path):
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


def bedfile_enforce_constant_size(bed_path, output_path, window):
    # load bed file

    f = open(bed_path, 'rb')
    df = pd.read_table(f, header=None)
    chrom = df[0].as_matrix().astype(str)
    start = df[1].as_matrix()
    end = df[2].as_matrix()

    # calculate center point and create dataframe
    middle = np.round((start + end)/2).astype(int)
    half_window = np.round(window/2).astype(int)

    # calculate new start and end points
    start = middle - half_window
    end = middle + half_window

    # filter any negative start positions
    index = np.where(middle - half_window > 0)[0]
    data = {}
    for i in range(len(df.columns)):
        data[i] = df[i].as_matrix()[index]
    data[1] = start[index]
    data[2] = end[index]

    # create new dataframe
    df_new = pd.DataFrame(data);

    # save dataframe with fixed width window size to a bed file
    df_new.to_csv(output_path, sep='\t', header=None, index=False)



def extract_fasta_sequences(bed_path, fasta_path, genome_path, window=None):
    """process bed file to a constand window and then extract sequences
        from reference genome 
    """

    # extract sequences from reference genome and store in fasta file
    bedtools.getfasta(bed_path, fasta_path, genome_path)

    # parse sequence and chromosome from fasta file
    sequences = parse_fasta_sequences(fasta_path)    

    return sequences


def process_background_sequences(bed_file_path, genome_path, window=None, background='dinuc', verbose=0):

    # generate different backgrounds
    #backgrounds = ['genome', 'dinuc', 'random']

    if background == 'nonspecific':
        # background sequences derived from genomic sequences
        neg_sequences = extract_fasta_sequences(bed_file_path, genome_path)

    elif background == 'dinuc':
        # background sequences derived from dinucleotide shuffle
        background_path = bed_file_path + '_dinuc.fa'
        meme.shuffle(bed_file_path+'.fa ', background_path, kmer=2, verbose=verbose)

        # process background sequences
        neg_sequences = parse_fasta_sequences(background_path)

    elif background == 'random':
        # background sequences derived from random shuffle
        background_path = bed_file_path + '_random.fa'
        meme.shuffle(bed_file_path+'.fa ', background_path, kmer=1, verbose=verbose)
        
        # process background sequences
        neg_sequences = parse_fasta_sequences(background_path)

    return neg_sequence
    

def convert_sequences_one_hot(sequences):
    """convert a sequence into a 1-hot representation"""
    def convert_one_hot(seq):
        nucleotide = 'ACGT'
        N = len(seq)
        one_hot_seq = np.zeros((4,N))
        for i in range(N):         
            index = [j for j in range(4) if seq[i] == nucleotide[j]]
            one_hot_seq[index,i] = 1
        return one_hot_seq

    # convert sequences to one-hot representation
    one_hot = []
    for seq in sequences:
        one_hot.append(convert_one_hot(seq))
    one_hot = np.array(one_hot)
    return one_hot



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

    return train, valid, test


def save_dataset_hdf5(savepath, train, valid, test):
    # save datasets as hdf5 file
    f = h5py.File(savepath, "w")
    dset = f.create_dataset("X_train", data=train[0].astype(np.float32), compression="gzip")
    dset = f.create_dataset("Y_train", data=train[1].astype(np.float32), compression="gzip")
    dset = f.create_dataset("X_valid", data=valid[0].astype(np.float32), compression="gzip")
    dset = f.create_dataset("Y_valid", data=valid[1].astype(np.float32), compression="gzip")
    dset = f.create_dataset("X_test", data=test[0].astype(np.float32), compression="gzip")
    dset = f.create_dataset("Y_test", data=test[1].astype(np.float32), compression="gzip")
    f.close()
    


def load_dataset_hdf5(file_path, dataset_name=None):

    dataset = h5py.File(file_path, 'r')
    
    if not dataset_name:
        # load set A data
        X_train = np.array(dataset['X_train']).astype(np.float32)
        Y_train = np.array(dataset['Y_train']).astype(np.float32)
        X_valid = np.array(dataset['X_valid']).astype(np.float32)
        Y_valid = np.array(dataset['Y_valid']).astype(np.float32)    
        X_test = np.array(dataset['X_test']).astype(np.float32)
        Y_test = np.array(dataset['Y_test']).astype(np.float32)

    else:
        X_train = np.array(dataset['/'+dataset_name+'/X_train']).astype(np.float32)        
        y_train = np.array(dataset['/'+dataset_name+'/Y_train']).astype(np.float32)
        X_valid = np.array(dataset['/'+dataset_name+'/X_valid']).astype(np.float32)       
        y_valid = np.array(dataset['/'+dataset_name+'/Y_valid']).astype(np.float32)
        X_test = np.array(dataset['/'+dataset_name+'/X_test']).astype(np.float32)
        y_test = np.array(dataset['/'+dataset_name+'/Y_test']).astype(np.float32)


    # add another dimension to make a 4d tensor
    X_train = np.expand_dims(X_train, axis=3).transpose([0, 2, 3, 1])
    X_test = np.expand_dims(X_test, axis=3).transpose([0, 2, 3, 1])
    X_valid = np.expand_dims(X_valid, axis=3).transpose([0, 2, 3, 1])

    # dictionary for each dataset
    train = {'inputs': X_train, 'targets': y_train}
    valid = {'inputs': X_valid, 'targets': y_valid}
    test = {'inputs': X_test, 'targets': y_test}

    return train, valid, test


def dataset_keys_hdf5(file_path):

    dataset = h5py.File(file_path, 'r')

    return dataset.keys()
