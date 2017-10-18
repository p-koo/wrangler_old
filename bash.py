from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys


def correlation_filter(bed_path, output_path, signals=[5, 12], equality='<', threshold=1.0, verbose=1):
    """filter based on: (2*(rep1 - rep2)/(rep1 + rep2))**2 < 1.0"""

    signal1 = str(signals[0])
    signal2 = str(signals[1])
    cmd = "cat "+bed_path+" | awk '{ if ((2*($"+signal1+"-$"+signal2+")/($"+signal1+"+$"+signal2+"))^2 "+equality+" "+str(threshold)+") print $0}' > "+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def positive_filter(bed_path, output_path, signal=6, threshold=3.0, verbose=1):
    """filter based on: signal > threshold"""

    cmd = "cat "+bed_path+" | awk '{ if ($"+str(signal)+" > "+str(threshold)+") print $0}' > "+output_path
    if verbose:
        print('>>' + cmd)
    os.system(cmd)


def rename_data(file_path, file_names, names, ext):
    """rename files """

    # generate new names
    new_names = []
    for i in range(len(file_names)):
        name = ''
        for j in range(len(names)):
            name = name + names[j][i] + '_'
        new_names.append(name[:-2])

    for i, file_name in enumerate(file_names):
        old_path = os.path.join(file_path, file_name+'.'+ext)
        new_path = os.path.join(file_path, new_names[i]+'.'+ext)
        os.command('mv ' + old_path + ' ' + new_path)

    return new_names


def move_files(file_names, old_path, new_path, ext=None):
    """move files to different directory"""

    # generate new names
    if ext:
   		for i in range(len(file_names)):
   			file_names[i] = file_names[i] + '.' + ext

    for i, file_name in enumerate(file_names):
        os.command('mv ' + os.path.join(old_path, file_name) + ' ' + os.path.join(new_path, file_name))


def convert_gtf_to_bed(file_path, verbose=1):
    cmd = """ awk '{if ($3 == "transcript") print $1"\t"$4"\t"$5"\t"$12"\t"($5-$4)"\t"$7}' gencode.v26.annotation.gtf | sed 's:"::g' | sed 's:;::g' >> gencode.v26.annotation.bed"""
    if verbose:
        print('>>' + cmd)
    os.command(cmd)
