from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys
import numpy as np
from Bio import SearchIO


def build(file_path, model_path, options=[], verbose=1):
    """build an hmmer model"""

    options_str = ''
    for option in options:
        if isinstance(option, (list, tuple)):
            options_str += option[0] + ' ' + option[1] + ' '
        else:
            options_str += option+ ' '

        rep_path = rep_paths

    cmd = 'hmmbuild ' + options_str + model_path + ' ' + file_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def search(file_path, model_path, output_path, options=[], verbose=1):
    """use a hmmer model to search for homologs across sequences in file_path"""

    options_str = ''
    for option in options:
        if isinstance(option, (list, tuple)):
            options_str += option[0] + ' ' + option[1] + ' '
        else:
            options_str += option+ ' '

    cmd = 'hmmsearch ' + options_str + ' -o ' + output_path + ' ' + model_path + ' ' + file_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def parse_results(file_path):
    """parse hmmer results for names, scores, and e-values"""

    with open(file_path, 'rU') as f:
        for record in SearchIO.parse(f, 'hmmer3-text'):
            hits = record.hits
            num_hits = len(hits)
            names = [[]]*num_hits
            scores = np.zeros(num_hits)
            evalues = np.zeros(num_hits)
            for i in range(num_hits):
                names[i] = hits[i].id
                scores[i] = hits[i].bitscore
                evalues[i] = hits[i].evalue
    return names, scores, evalues


def reformat(file_path, output_path, file_format='stockholm', verbose=1):
    """use easel to reformat a file into another format"""

    #esl-reformat -o PF00011_train.sto stockholm PF00011_full_train.fa
    cmd = 'esl-reformat -o ' + output_path + ' ' + file_format + ' ' + file_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)
