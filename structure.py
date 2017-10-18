from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import pandas as pd
import numpy as np


def predict_structure(fasta_path, profile_path, window):
    """predict secondary structure profiles with RNAplfold modified scripts"""

    E_path = profile_path+'E_profile.txt'
    # predict external loops
    os.system('E_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+E_path)

    # predict hairpin loops
    H_path = profile_path+'H_profile.txt'
    os.system('H_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+H_path)

    # predict internal loops
    I_path = profile_path+'I_profile.txt'
    os.system('I_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+I_path)

    # predict multi-loops
    M_path = profile_path+ 'M_profile.txt'
    os.system('M_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+M_path)



def merge_structural_profile(profile_path, merged_path):
    """merge the secondary structure profiles into a single file"""
    def list_to_str(lst):
        ''' Given a list, return the string of that list with tab separators
        '''
        return reduce( (lambda s, f: s + '\t' + str(f)), lst, '')

    # external loop profile
    E_path = profile_path+'E_profile.txt'
    fEprofile = open(E_path)
    Eprofiles = fEprofile.readlines()

    # hairpin loop profiles
    H_path = profile_path+'H_profile.txt'
    fHprofile = open(H_path)
    Hprofiles = fHprofile.readlines()

    # internal loop profiles
    I_path = profile_path+'I_profile.txt'
    fIprofile = open(I_path)
    Iprofiles = fIprofile.readlines()

    # multi-loop profiles
    M_path = profile_path+ 'M_profile.txt'
    fMprofile = open(M_path)
    Mprofiles = fMprofile.readlines()

    num_seq = int(len(Eprofiles)/2)

    # parse into a single file
    fhout = open(merged_path, 'w')
    for i in xrange(num_seq):
        id = Eprofiles[i*2].split()[0]
        fhout.write(id+'\n')
        H_prob =  Hprofiles[i*2+1].split()
        I_prob =  Iprofiles[i*2+1].split()
        M_prob =  Mprofiles[i*2+1].split()
        E_prob =  Eprofiles[i*2+1].split()
        P_prob = map( (lambda a, b, c, d: 1-float(a)-float(b)-float(c)-float(d)), H_prob, I_prob, M_prob, E_prob)
        fhout.write(list_to_str(P_prob[:len(P_prob)])+'\n')
        fhout.write(list_to_str(H_prob[:len(P_prob)])+'\n')
        fhout.write(list_to_str(I_prob[:len(P_prob)])+'\n')
        fhout.write(list_to_str(M_prob[:len(P_prob)])+'\n')
        fhout.write(list_to_str(E_prob[:len(P_prob)])+'\n')
    fhout.close()

    return num_seq


def extract_structural_profile(merged_path, num_seq, window):
    """extract secondary structure profiles from a merged file and return a
       numpy array """

    # parse further and load structural profile as np.array
    f = open(merged_path, 'r')
    structure = []
    for i in xrange(num_seq):
        seq = f.readline()
        paired = f.readline().strip().split('\t')
        hairpin = f.readline().strip().split('\t')
        internal = f.readline().strip().split('\t')
        multi = f.readline().strip().split('\t')
        external = f.readline().strip().split('\t')

        paired = np.array(paired).astype(np.float32)
        hairpin = np.array(hairpin).astype(np.float32)
        internal = np.array(internal).astype(np.float32)
        multi = np.array(multi).astype(np.float32)
        external = np.array(external).astype(np.float32)

        # pad sequences
        seq_length = len(paired)
        offset1 = int((window - seq_length)/2)
        offset2 = window - seq_length - offset1
        struct = np.array([paired, hairpin, internal, multi, external])
        num_dims = struct.shape[0]
        if offset1:
            struct = np.hstack([np.zeros((num_dims,offset1)), struct])
        if offset2:
            struct = np.hstack([struct, np.zeros((num_dims,offset2))])
        structure.append(struct)

    return np.array(structure)



def RNAplfold_profile(fasta_path, profile_path, window):
    """predict secondary structure profiles for sequences in a fasta file"""

    # predict secondary structural profiles
    predict_structure(fasta_path, profile_path, window)

    # generate merged secondary structure profile
    merged_path = profile_path+'_structure_profiles.txt'
    num_seq = merge_structural_profile(profile_path, merged_path)

    # extract secondary structure profiles
    structure = extract_structural_profile(merged_path, num_seq, window)

    return structure
