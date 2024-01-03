#!/usr/bin/env python
# coding: utf-8

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
from math import log

bins=len(np.arange(1,11.5,0.5))+1	 #fragment windows
ancestries=3 # 
samplePOP=100
TOTALPOP=1000
GEN=19
scale=1000
SCALE=scale

POPname=str(sys.argv[1])

def resample_array(myarray):
    myarray_indices=np.array(range(len(myarray)))
    my_array_reps = np.repeat(myarray_indices,myarray)
    my_array_resampled=np.random.choice(my_array_reps,size=len(my_array_reps),replace=True)
    resampled_array_counts=[]
    for i in myarray_indices:
        resampled_array_counts.append(np.count_nonzero(my_array_resampled==i))
    return(resampled_array_counts)

with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_"+str(POPname)+"_chr_X_onlyfemales", 'rb') as f:
    all_lengths_counts_X=pickle.load(f)

X_matrix=all_lengths_counts_X

bootstrap_X=[]
for i in range(1000):
    print(i)
    bootstrap_X.append(np.apply_along_axis(resample_array, 0, X_matrix))

with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/Bootstrap_"+str(POPname)+"_log_rfmix_maxanc_X_onlyfemales", 'wb') as f:    
    pickle.dump([bootstrap_X], f)

