from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, time, h5py
import numpy as np
import pandas as pd


def parse_metadata(meta_data_path, parse_list=None, genome=None, file_format=None):
    """Parse the metadata.tsv file that comes with a download data from ENCODE"""
    if not parse_list:
        parse_list = ['File accession', 'Experiment accession', 'Biological replicate(s)',
                      'Experiment target', 'Biosample term name']

    # open meta file
    df = pd.read_csv(meta_data_path, sep='\t')

    # parse target names
    values = {}
    for item in parse_list:
        values[item] = df[item].as_matrix()

    if genome:
        # filter for GRCh38
        assembly = df['Assembly'].as_matrix()
        genome_index = np.where(assembly == genome)[0]
        for item in parse_list:
            values[item] = values[item][genome_index]

    if file_format:
        values['File format'] = df['File format'].as_matrix()
        if genome:
            values['File format'] = values['File format'][genome_index]

        # filter for bam files
        format_index = np.where(values['File format'] == file_format)[0]
        for item in parse_list:
            values[item] = values[item][format_index]

    return values


def match_control_experiments(targets, cell_types, targets_control, cell_types_control):
    """find the indices of the control experiments that match the experiment"""

    # find which experiments match control experiments via target and cell-type

    index = np.where(targets_control[:targets_control.index(' ')] == targets)[0]
    match_index = index[np.where(cell_types_control == cell_types[index])[0]]

    return match_index
