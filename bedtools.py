from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys

def bam_to_bed(bam_path, bed_path, verbose=1):
    """ Bam file to bed file conversion """

    command = 'bedtools bamtobed -i ' + bam_path + ' > ' + bed_path
    if verbose:
        print('>>' + command)
    
    os.system(command)


def sort(file_path, output_path, verbose=1):    
    # merge bed file coordinates
    cmd = 'sortBed -i '+file_path+' > '+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def merge(file_path, output_path, verbose=1):    
    # merge bed file coordinates
    cmd = 'bedtools merge -i '+file_path+' -s > '+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def nonoverlap(file_path, rep_paths, ouput_path, options=['-wa', '-s'], verbose=1):
    # find CLIP-peaks that don't overlap with background peaks

    options_str = ''
    for option in options:
        options_str += option+ ' '

    # make string list of replicate paths
    rep_path = ''
    if isattribute(rep_paths, (list, tuple)):
        for path in rep_paths:
            rep_path = rep_path + ' '
    else:
        rep_path = rep_paths

    cmd = 'bedtools intersect '+options_str+'-a '+file_path+' -b '+
                rep_path+'-v > '+ouput_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def overlap(file_path, rep_paths, output_path, options=['-wa', '-wb', '-s'], verbose=1):
    # find CLIP-peaks that overlap between replicates 

    options_str = ''
    for option in options:
        options_str += option+ ' '

    # make string list of replicate paths
    rep_path = ''
    if isattribute(rep_paths, (list, tuple)):
        for path in rep_paths:
            rep_path = rep_path + ' '
    else:
        rep_path = rep_paths

    cmd = 'bedtools intersect '+options_str+'-a '+file_path+' -b '+rep_path+' > '+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def getfasta(bed_path, output_fasta, genome_path)    
    # extract sequences from bed files and save as fasta file 
    os.system('bedtools getfasta -s -fi '+genome_path+' -bed '+bed_path+' -fo '+output_fasta)




"""

def filter_overlap_signalValue(nonspecific_overlap_path, filter_path, threshold=1.0, verbose=1):
    # filter based on: (2*(rep1 - rep2)/(rep1 + rep2))**2 < 1.0
        
    df_nonspecific = pd.read_csv(nonspecific_overlap_path, sep='\t', header=None)
    df_nonspecific.head()
    rep1_nonspecific = df_nonspecific[4].as_matrix()
    rep2_nonspecific = df_nonspecific[11].as_matrix()
    chr_nonspecific = df_nonspecific[7].as_matrix()
    start_nonspecific = df_nonspecific[8].as_matrix()
    end_nonspecific = df_nonspecific[9].as_matrix()

    df_filtered = pd.read_csv(filter_path, sep='\t', header=None)
    df_filtered.head()
    rep1_filtered = df_filtered[4].as_matrix()
    rep2_filtered = df_filtered[11].as_matrix()
    chr_filtered = df_filtered[0].as_matrix()
    start_filtered = df_filtered[1].as_matrix()
    end_filtered = df_filtered[2].as_matrix()


    error = 2*(rep1_nonspecific-rep2_nonspecific)/(rep1_nonspecific+rep2_nonspecific)

    filtered_index = []
    for i in range(len(chr_filtered)):
        index = np.where((chr_filtered[i] == chr_nonspecific) & 
                         (start_filtered[i] == start_nonspecific) & 
                         (end_filtered[i] == end_nonspecific))[0]
        if index.any():
            for j in index:
                filtered_index.append(j)
    filtered_index = np.unique(filtered_index)
"""
