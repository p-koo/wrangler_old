from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys
import numpy as np
import pandas as pd

def bam_to_bed(bam_path, bed_path, verbose=1):
    """ Bam file to bed file conversion """

    command = 'bedtools bamtobed -i ' + bam_path + ' > ' + bed_path
    if verbose:
        print('>>' + command)

    os.system(command)


def sort(file_path, output_path, verbose=1):
    """sort bedfile"""

    cmd = 'sortBed -i '+file_path+' > '+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def merge(file_path, output_path, verbose=1):
    """merge bed file coordinates"""

    cmd = 'bedtools merge -i '+file_path+" -s | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\tMergedPeak\"NR\"\\t\"($3-$2)\"\\t\"$4}' > "+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def nonoverlap(file_path, rep_paths, ouput_path, options=['-wa', '-s'], verbose=1):
    """find non-overlapping peaks in file_path and rep_paths (can be multiple files)"""

    options_str = ''
    for option in options:
        options_str += option+ ' '

    # make string list of replicate paths
    rep_path = ''
    if isinstance(rep_paths, (list, tuple)):
        for path in rep_paths:
            rep_path = rep_path + path + ' '
    else:
        rep_path = rep_paths

    cmd = 'bedtools intersect '+options_str+'-a '+file_path+' -b '+rep_path+' -v > '+ouput_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def overlap(file_path, rep_paths, output_path, options=['-wa', '-wb', '-s'], verbose=1):
    """find overlapping peaks in file_path and rep_paths (can be multiple files)"""

    options_str = ''
    for option in options:
        options_str += option+ ' '

    # make string list of replicate paths
    rep_path = ''
    if isinstance(rep_paths, (list, tuple)):
        for path in rep_paths:
            rep_path = rep_path + path + ' '
    else:
        rep_path = rep_paths

    cmd = 'bedtools intersect '+options_str+'-a '+file_path+' -b '+rep_path+' > '+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def enforce_constant_size(bed_path, output_path, window, compression=None):
    """generate a bed file where all peaks have same size centered on original peak"""

    # load bed file
    f = open(bed_path, 'rb')
    df = pd.read_table(f, header=None, compression=compression)
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


def to_fasta(bed_path, output_fasta, genome_path):
    """extract sequences from bed files and save as fasta file """

    os.system('bedtools getfasta -s -fi '+genome_path+' -bed '+bed_path+' -fo '+output_fasta)


def count_bed_entries(file_path):
    with open(file_path, 'r') as f:
        counts = 0
        for line in f:
            counts += 1
    return counts
