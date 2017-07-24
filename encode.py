from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, time, h5py
import numpy as np
import pandas as pd


def parse_metadata(meta_data_path, parse_list=None, genome=None, file_format=None):
    if not parse_list:
        parse_list = ['File accession', 'File format', 
                      'Experiment accession', 'Biological replicate(s)',
                      'Experiment target', 'Biosample term name']

    # open meta file
    df = pd.read_csv(meta_data_path, sep='\t')

    # parse target names
    values = []
    for item in parse_list:
        values.append(df[item].as_matrix())

    """
    targets = df['Experiment target'].as_matrix()
    cell_types = df['Biosample term name'].as_matrix()
    file_names = df['File accession'].as_matrix()
    file_format = df['File format'].as_matrix()
    experiments = df['Experiment accession'].as_matrix()
    replicates = df['Biological replicate(s)'].as_matrix()
    #antibody = df['Antibody accession'].as_matrix()
    #links = df['File download URL'].as_matrix()
    """

    if genome:
        # filter for GRCh38
        assembly = df['Assembly'].as_matrix()
        genome_index = np.where(assembly == genome)[0]
        for i in range(len(values)):
            values[i] = values[i][genome_index]
        """    
        targets = targets[genome_index]
        cell_types = cell_types[genome_index]
        file_names = file_names[genome_index]
        file_format = file_format[genome_index]
        experiments = experiments[genome_index]
        replicates = replicates[genome_index]
        #antibody = antibody[genome_index]
        #links = links[genome_index]
        """

    if file_format:
        # filter for bam files
        format_index = np.where(file_format == 'bam')[0]
        for i in range(len(values)):
            values[i] = values[i][format_index]
        """
        targets = rbp[format_index]
        cell_types = cell_types[format_index]
        file_names = file_names[format_index]
        experiments = experiments[format_index]
        replicates = replicates[format_index]
        #antibody = antibody[format_index]
        #links = links[format_index]
        """
    return values


def parse_list(values):
            parse_list = ['File accession', 'File format', 
                      'Experiment accession', 'Biological replicate(s)',
                      'Experiment target', 'Biosample term name']

    file_names = values[0]
    file_formats = values[1]
    experiments = values[2]
    replicates = values[3]
    targets = values[4]
    cell_types = values[5]

    return file_names, file_formats, experiments, replicates, targets, cell_types


def match_control_experiments(targets, cell_types, targets_control, cell_types_control):
    # find which experiments match control experiments via target and cell-type 

    index = np.where(targets_control[:targets_control.index(' ')] == targets)[0]
    match_index = index[np.where(cell_types_control == cell_types[index])[0]]

    return match_index







