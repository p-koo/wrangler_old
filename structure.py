from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import pandas as pd
import numpy as np


def RNAplfold_profile(fasta_path, output_path, window):

    def list_to_str(lst):
        ''' Given a list, return the string of that list with tab separators 
        '''
        return reduce( (lambda s, f: s + '\t' + str(f)), lst, '')
    
    # external loop profile
    E_path = os.path.join(output_path, 'E_profile.txt')
    os.system('E_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+E_path)
    fEprofile = open(E_path)
    Eprofiles = fEprofile.readlines()

    # hairpin loop profiles
    H_path = os.path.join(output_path, 'H_profile.txt')
    os.system('H_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+H_path)
    fHprofile = open(H_path)
    Hprofiles = fHprofile.readlines()

    # internal loop profiles
    I_path = os.path.join(output_path, 'I_profile.txt')
    os.system('I_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+I_path)
    fIprofile = open(I_path)
    Iprofiles = fIprofile.readlines()

    # multi-loop profiles
    M_path = os.path.join(output_path, 'M_profile.txt')
    os.system('M_RNAplfold -W '+str(window)+' -u 1 <'+fasta_path+' >'+M_path)
    fMprofile = open(M_path)
    Mprofiles = fMprofile.readlines()

    num_seq = int(len(Eprofiles)/2)
    
    # parse into a single file
    profile_path = os.path.join(output_path,'structure_profiles.txt')
    fhout = open(profile_path, 'w')
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

    # parse further and load structural profile as np.array
    f = open(profile_path, 'r')
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

        structure.append(np.array([paired, hairpin, internal, multi, external]))

    return np.array(structure)