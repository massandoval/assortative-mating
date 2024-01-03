#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras.layers import Reshape, Input, Flatten, concatenate, Dense, BatchNormalization ,Conv1D,Conv2D,Conv3D,MaxPool1D, MaxPool2D,MaxPool3D, GlobalMaxPool2D, Dropout
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.models import Model
from tensorflow.keras.utils import plot_model
import tensorflow.math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import scipy
import sklearn
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import os
import plotnine as pn
import pydot
import pydotplus
import random
import pickle
import fileinput
from math import log

SCALE=100000

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

target_panel=np.loadtxt("/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/target_panel.txt",delimiter='\t',dtype="str")
all_inds=target_panel
rfmix_chr=np.loadtxt("/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/HGDP_1000gNYC_3pop_target_ref_output_rfmix154_short_B_X.1.ForwardBackward.txt")
females_1=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="F"))[0]*6
females_2=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="F"))[0]*6+1
females_3=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="F"))[0]*6+2
females_4=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="F"))[0]*6+3
females_5=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="F"))[0]*6+4
females_6=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="F"))[0]*6+5
males_1=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="M"))[0]*6
males_2=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="M"))[0]*6+1
males_3=np.where((all_inds[:,1]=="CLM")&(all_inds[:,2]=="M"))[0]*6+2
indices=np.sort(np.concatenate([females_1,females_2,females_3,females_4,females_5,females_6]))
pos_chr=np.loadtxt("/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/snps_pos_rfmix_short_B_X")
rfmix_chr_pop=rfmix_chr[:,indices]
n_anc=3
n_cells=22
n_inds=int(len(indices)/n_anc)
n_windows=np.shape(rfmix_chr_pop)[0]
rfmix_chr_pop_3d=rfmix_chr_pop.reshape([n_windows,n_inds,n_anc])
rfmix_chr_pop_3d_rearr = np.moveaxis(rfmix_chr_pop_3d, 1, 0)
output_get_frag = [get_fragments_lengths(i) for i in rfmix_chr_pop_3d_rearr]
lengths = np.array([output_get_frag[i][0] for i in range(len(output_get_frag))])
anc_total = np.array([output_get_frag[i][1] for i in range(len(output_get_frag))])
lengths_counts=np.moveaxis(lengths.reshape(n_inds*n_anc,n_cells),[0,1],[1,0])

with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_CLM_chr_X_onlyfemales", 'wb') as f: 
    pickle.dump(lengths_counts,f)

with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_rfmix_CLM_maxanc_chr_X_onlyfemales", 'wb') as f:
    pickle.dump(anc_total,f)

