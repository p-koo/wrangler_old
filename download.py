#!/bin/python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys

def reference_genome(genome, data_path):
    """download the reference genome"""

    if genome == 'hg19':
        genome_link = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit'
    elif genome == 'hg38':
        genome_link = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit'
    elif genome == 'mm10':
        genome_link = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit'

    # download file
    bit_path = os.path.join(data_path, genome)
    download_file(genome_link, bit_path+'.2bit')

    # convert 2bit to fasta
    genome_path = convert_2bit_to_fa(data_path, genome)


def download_file(link, output_path=None):
    """ download files with wget """

    if output_path:
        if not os.path.isfile(output_path):
            os.system('wget -O ' + output_path + ' ' + link)
    else:
        os.system('wget ' + link)


def convert_2bit_to_fa(data_path, genome):
    """Convert 2bit to fasta. Download converter program if not in data_path"""
    # download twoBitToFa file to convert genomes
    twobit_path = os.path.join(data_path, 'twoBitToFa')
    if not os.path.isfile(twobit_path):
        link = 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
        download_file(link, twobit_path)
        os.system('chmod +x ' + twobit_path)

    # convert 2bit files to fasta
    genome_path = os.path.join(data_path, genome)
    os.system(twobit_path + ' ' + genome_path+'.2bit ' + genome_path+'.fa')


def reference_transcriptome(transcriptome, data_path):
    """Download reference transcriptome if not in data_path"""
    if transcriptome == 'hg38':
        genome_link = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz'

    # download file
    genome_path = os.path.join(data_path, 'gencode.v26.annotation.gtf.gz')
    download_file(genome_link, genome_path)
