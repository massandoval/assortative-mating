#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import scipy
import os
import pydot
import pydotplus
import random
import pickle


# In[2]:




# In[3]:


bins=len(np.arange(1,11.5,0.5))+1	 #fragment windows
ancestries=3 # 
samplePOP=100
TOTALPOP=1000
GEN=19
scale=1000
SCALE=scale


# In[4]:


bins


# In[5]:


POPname=str(sys.argv[1])


# In[6]:


def resample_array(myarray):
    myarray_indices=np.array(range(len(myarray)))
    my_array_reps = np.repeat(myarray_indices,myarray)
    my_array_resampled=np.random.choice(my_array_reps,size=len(my_array_reps),replace=True)
    resampled_array_counts=[]
    for i in myarray_indices:
        resampled_array_counts.append(np.count_nonzero(my_array_resampled==i))
    return(resampled_array_counts)


# In[8]:


all_lengths_counts =[]
for i in range(22):
    print(i)
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_"+str(POPname)+"_chr_"+str(i+1), 'rb') as f:
        lengths_counts=pickle.load(f)
    all_lengths_counts.append(lengths_counts)
with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_"+str(POPname)+"_chr_X", 'rb') as f:
    all_lengths_counts_X=pickle.load(f)


# In[9]:


Aut_matrix=np.apply_along_axis(sum,0,all_lengths_counts)
X_matrix=all_lengths_counts_X


# In[10]:


bootstrap_Aut=[]
bootstrap_X=[]
for i in range(1000):
    print(i)
    bootstrap_Aut.append(np.apply_along_axis(resample_array, 0, Aut_matrix))
    bootstrap_X.append(np.apply_along_axis(resample_array, 0, X_matrix))


# In[ ]:


with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/Bootstrap_"+str(POPname)+"_log_rfmix_maxanc", 'wb') as f:    
    pickle.dump([bootstrap_Aut], f)
with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/Bootstrap_"+str(POPname)+"_log_rfmix_maxanc_X", 'wb') as f:    
    pickle.dump([bootstrap_X], f)

