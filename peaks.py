from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys


def call_peaks(file_path, output_path, bin_size=20, p_value=0.01, verbose=1):
    """ Piranha peak calling """

    cmd = 'Piranha -s -b '+str(bin_size)+' -d ZeroTruncatedNegativeBinomial -p '+str(p_value)+' -o '+output_path+' '+file_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def call_covariate_peaks(file_path, output_path, bin_size=20, p_value=0.01, covariates_path='', verbose=1):
    """ Piranha peak calling with covariates"""

    cmd = 'Piranha -s -i '+str(bin_size)+' -b '+str(bin_size)+ \
            ' -d ZeroTruncatedNegativeBinomialRegression -p ' \
            +str(p_value)+' -o '+output_path+' '+file_path
    if covariates_path:
         cmd = cmd + ' ' + covariates_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


"""

def generate_covariate_bed(pirana_bed, tpm_bed, outfile, verbose=1):

    command = 'bedtools intersect -a ' + str(pirana_bed) + ' -b ' + str(tpm_bed)
        + ' -wa -wb -s | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$12"\t"$6}\' > ' + str(outfile)

    if verbose:
        print('>>' + command)

    os.system(command)

"""
