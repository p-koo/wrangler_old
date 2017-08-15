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
    data = {}
    for i in range(len(df.columns)):
        data[i] = df[i].as_matrix()
    data[1] = start
    data[2] = end

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


def generate_fasta(sequences, fasta_path):
    """generate fasta file from an array of sequences
    """

    with open(fasta_path, 'w+') as f:
        for i in xrange(len(sequences)):
            f.write('>seq '+str(i))
            f.write('\n')
            f.write(sequences[i])
            f.write('\n')


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
    

def convert_sequences_one_hot(sequences, nucleotide='ACGT'):
    """convert a sequence into a 1-hot representation"""
    def convert_one_hot(seq):
        N = len(seq)
        one_hot_seq = np.zeros((4,N))
        for i in range(N):         
            index = [j for j in range(4) if seq[i].upper() == nucleotide[j]]
            one_hot_seq[index,i] = 1
        return one_hot_seq

    # convert sequences to one-hot representation
    one_hot = []
    for seq in sequences:
        one_hot_seq = convert_one_hot(seq)
        if np.sum(one_hot_seq) == len(seq):
            one_hot.append(one_hot_seq)
    one_hot = np.array(one_hot)
    return one_hot


def convert_one_hot_sequences(one_hots, nucleotide='ACGT'):
    """convert 1-hot representation to sequence"""

    seq_length = one_hots.shape[1]
    seq = []
    for one_hot in one_hots:
        sequence = np.array(['']*seq_length)    
        for i in range(4):
            index = np.where(one_hot[:,i] == 1)[0]
            sequence[index] = nucleotide[i]

        seq.append(''.join(sequence))
    return seq




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
        y_train = np.array(dataset['Y_train']).astype(np.float32)
        X_valid = np.array(dataset['X_valid']).astype(np.float32)
        y_valid = np.array(dataset['Y_valid']).astype(np.float32)    
        X_test = np.array(dataset['X_test']).astype(np.float32)
        y_test = np.array(dataset['Y_test']).astype(np.float32)

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
    keys = []
    for key in dataset.keys():
        keys.append(str(key))

    return np.array(keys)


def conservation_bed(bed_path, conservation_path):
    
    df = pd.read_csv(bed_path, sep='\t', header=None)
    chrom = df[0].as_matrix()
    start = df[1].as_matrix()
    end = df[2].as_matrix()
    strand = df[3].as_matrix()

    dataset = h5py.File(conservation_path, 'r')

    conservation = []
    good_index = []
    for i in range(len(chrom)):
        if str(chrom[i]) in dataset.keys():
            good_index.append(i)
            conservation.append(np.array(dataset[chrom[i]][start[i]:end[i]]))     
            
    good_index = np.array(good_index)
    conservation = np.array(conservation)
    return conservation, good_index



def conservation_bed_all2(bed_path, phylop_path, phastcons_path):
    
    df = pd.read_csv(bed_path, sep='\t', header=None)
    chrom = df[0].as_matrix()
    start = df[1].as_matrix()-1
    end = df[2].as_matrix()-1
    strand = df[3].as_matrix()

    dataset1 = h5py.File(phylop_path, 'r')
    dataset2 = h5py.File(phastcons_path, 'r')

    conservation1 = []
    conservation2 = []
    for i in range(len(chrom)):
        conservation1.append(np.array(dataset1[chrom[i]][start[i]:end[i]]))     
        conservation2.append(np.array(dataset2[chrom[i]][start[i]:end[i]]))     

    conservation1 = np.array(conservation1)
    conservation2 = np.array(conservation2)

    return conservation1, conservation2


def conservation_bed_all(bed_path, phylop_path, phastcons_path, good_index):
    
    df = pd.read_csv(bed_path, sep='\t', header=None)
    chrom = df[0].as_matrix()
    start = df[1].as_matrix()-1
    end = df[2].as_matrix()-1
    strand = df[3].as_matrix()

    dataset1 = h5py.File(phylop_path, 'r')
    dataset2 = h5py.File(phastcons_path, 'r')

    conservation1 = []
    conservation2 = []
    for j, i in enumerate(good_index):
        conservation1.append(np.array(dataset1[chrom[i]][start[i]:end[i]]))     
        conservation2.append(np.array(dataset2[chrom[i]][start[i]:end[i]]))     

    conservation1 = np.array(conservation1)
    conservation2 = np.array(conservation2)

    return conservation1, conservation2

    

def count_bed_entries(file_path):
    with open(file_path, 'r') as f:
        counts = 0
        for line in f:
            counts += 1
    return counts



def prepare_data(train, struct=None, conservation=None):
    seq = train['inputs'][:,:,:,:4]
    if struct == 'pu':     
        structure = train['inputs'][:,:,:,4:9]
        paired = np.expand_dims(structure[:,:,:,0], axis=3)
        unpaired = np.expand_dims(np.sum(structure[:,:,:,1:], axis=3), axis=3)
        seq = np.concatenate([seq, paired, unpaired], axis=3)
    elif struct == 'all':    
        structure = train['inputs'][:,:,:,4:9]
        paired = np.expand_dims(structure[:,:,:,0], axis=3)
        HIME = structure[:,:,:,1:]
        seq = np.concatenate([seq, paired, HIME], axis=3)

    if conservation == 'phylop':
        phylop = np.expand_dims(train['inputs'][:,:,:,9], axis=3)
        phylop /= np.std(phylop)
        seq = np.concatenate([seq, phylop], axis=3)
    elif conservation == 'phastcon':
        phastcon = np.expand_dims(train['inputs'][:,:,:,10], axis=3)
        seq = np.concatenate([seq, phastcon], axis=3)
    elif conservation == 'all':
        phylop = np.expand_dims(train['inputs'][:,:,:,9], axis=3)
        phylop /= np.std(phylop)
        phastcon = np.expand_dims(train['inputs'][:,:,:,10], axis=3)
        seq = np.concatenate([seq, phylop, phastcon], axis=3)

    train['inputs']  = seq
    return train

    