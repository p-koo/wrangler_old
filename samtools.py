
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys

def parse_first_pair(bam_path, output_path, verbose=1):
    """Filter paired-end reads bam file for just the first pair"""
    cmd = 'samtools view -hu -f 64 ' + bam_path + ' -o ' + output_path
    if verbose:
        print('>>' + cmd)
    
    os.system(cmd)


def flagstat(bam_path):
	cmd = 'samtools flagstat '+bam_path
	os.system(cmd)