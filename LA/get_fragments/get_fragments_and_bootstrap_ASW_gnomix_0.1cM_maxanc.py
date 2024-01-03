#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import sys
import scipy
import os
import pydot
import pydotplus
import random
import pickle
import fileinput


SCALE=100000
chr = sys.argv[1]

def get_fragments_lengths(prov_vector):
    lengths_counts_3anc = []
    lengths_total_3anc = []
    for anc in range(3):
        startpos=pos_chr[0]
        oldpos = 0
        fragment = 0
        fragmentslengths = []
        #chr=1
        for pos in range(len(pos_chr)):
            if((np.argmax(prov_vector[pos])==anc)&(np.argmax(prov_vector[oldpos])!=anc)):
            #if ((prov_vector[pos]>=0.9)&(prov_vector[oldpos]<0.9)):
                startpos=pos_chr[pos]
            if((np.argmax(prov_vector[pos])!=anc)&(np.argmax(prov_vector[oldpos])==anc)):
            #if ((prov_vector[pos]<0.9)&(prov_vector[oldpos]>=0.9)):
                endpos=pos_chr[pos]
                fragment = endpos-startpos;
                fragmentslengths.append(fragment)
            if((np.argmax(prov_vector[pos])==anc)&(pos_chr[pos]==pos_chr[len(pos_chr)-1])):
            #if ((prov_vector[pos]>=0.9)&(pos_chr[pos]==pos_chr[len(pos_chr)-1])):
                endpos=pos_chr[pos]
                fragment = endpos-startpos;
                fragmentslengths.append(fragment)
            ##if (pos_chr[pos]==pos_chr[len(pos_chr)-1]):
                #chr=chr+1
            oldpos=pos
        if (len(fragmentslengths)==0):
             fragmentslengths.append(0.0)
        edges = np.concatenate([[0,],np.power(2,np.arange(1,11.5,0.5))])*SCALE
        windows=np.apply_along_axis(sum, 1, np.array([(i>=edges) for i in fragmentslengths]))
        lengths_counts=np.apply_along_axis(sum, 1,np.array([(windows==i) for i in range(1,(len(edges)+1))]))
        lengths_total=np.sum(fragmentslengths)
        lengths_counts_3anc.append(lengths_counts)
        lengths_total_3anc.append(lengths_total)
    return [lengths_counts_3anc, lengths_total_3anc]

target_panel=np.loadtxt("/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/inds_target_rename_rfmix1.5.4_AFRAME.txt",delimiter='_',dtype="str")
all_inds=target_panel

inds_1=np.where(all_inds[:,1]=="ASW")[0]*6
inds_2=np.where(all_inds[:,1]=="ASW")[0]*6+1
inds_3=np.where(all_inds[:,1]=="ASW")[0]*6+2
inds_4=np.where(all_inds[:,1]=="ASW")[0]*6+3
inds_5=np.where(all_inds[:,1]=="ASW")[0]*6+4
inds_6=np.where(all_inds[:,1]=="ASW")[0]*6+5
indices=np.sort(np.concatenate([inds_1,inds_2,inds_3,inds_4,inds_5,inds_6]))

filename="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/gnomix/query_results_0.1cM_"+str(chr)+"_AFRAME.fb"
gnomix = pd.read_csv(filename,sep="\t",skiprows=2,header=None)
gnomix_probs = np.array(gnomix.drop(gnomix.columns[range(4)], axis=1))
pos_chr=np.array(gnomix)[:,1]
gnomix_chr_pop=gnomix_probs[:,indices]
n_anc=3
n_cells=22
n_inds=int(len(indices)/n_anc)
n_windows=np.shape(gnomix_chr_pop)[0]
gnomix_chr_pop_3d=gnomix_chr_pop.reshape([n_windows,n_inds,n_anc])
gnomix_chr_pop_3d_rearr = np.moveaxis(gnomix_chr_pop_3d, 1, 0)
output_get_frag = [get_fragments_lengths(i) for i in gnomix_chr_pop_3d_rearr]
lengths = np.array([output_get_frag[i][0] for i in range(len(output_get_frag))])
anc_total = np.array([output_get_frag[i][1] for i in range(len(output_get_frag))])
lengths_counts=np.moveaxis(lengths.reshape(n_inds*n_anc,n_cells),[0,1],[1,0])

with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_gnomix_maxanc_ASW_chr_"+str(chr), 'wb') as f: 
    pickle.dump(lengths_counts,f)

with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_gnomix_ASW_maxanc_chr_"+str(chr), 'wb') as f:
    pickle.dump(anc_total,f)

