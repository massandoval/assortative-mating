
import numpy as np
import pandas as pd
import sys
import scipy
import os
import pydot
import pydotplus
import random
import pickle

#for POPname in ["ASW"]:
for POPname in ["CLM","MXL","PEL","PUR","ACB","ASW"]:
    all_lengths_counts =[]
    for i in range(22):
        print(i)
        with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_"+str(POPname)+"_chr_"+str(i+1), 'rb') as f:
            lengths_counts=pickle.load(f)
        all_lengths_counts.append(lengths_counts)
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_"+str(POPname)+"_chr_X", 'rb') as f:
        all_lengths_counts_X=pickle.load(f)
    Aut_matrix=np.apply_along_axis(sum,0,all_lengths_counts)
    X_matrix=all_lengths_counts_X
    Aut_matrix_pop=np.apply_along_axis(np.mean,1,Aut_matrix)
    X_matrix_pop=np.apply_along_axis(np.mean,1,X_matrix)
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_rfmix_maxanc_Aut_pop_'+str(POPname)+'.txt', np.round(Aut_matrix_pop,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_rfmix_maxanc_Aut_ind_'+str(POPname)+'.txt', np.round(Aut_matrix,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_rfmix_maxanc_X_pop_'+str(POPname)+'.txt', np.round(X_matrix_pop,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_rfmix_maxanc_X_ind_'+str(POPname)+'.txt', np.round(X_matrix,2), delimiter='\t', fmt='%f')


for POPname in ["ACB","ASW","CLM","MXL","PEL","PUR"]:
    all_length_anc =[]
    for i in range(22):
        print(i)
        with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_rfmix_"+str(POPname)+"_maxanc_chr_"+str(i+1), 'rb') as f:
            length_anc=pickle.load(f)
        all_length_anc.append(length_anc)
    
    all_length_anc_sum=np.apply_along_axis(sum, 0,all_length_anc)
       
    #all_length_anc_prop=np.apply_along_axis(rel_prop,1,all_length_anc_sum)
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_anc_prop_ind_rfmix_maxanc_Aut_'+str(POPname)+'_bp.txt', np.round(all_length_anc_sum), delimiter='\t', fmt='%f')
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_rfmix_"+str(POPname)+"_maxanc_chr_X", 'rb') as f:
        all_length_anc_X=pickle.load(f)
    
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_anc_prop_ind_rfmix_maxanc_X_'+str(POPname)+'_bp.txt', np.round(all_length_anc_X), delimiter='\t', fmt='%f')
   
for POPname in ["CLM","MXL","PEL","PUR","ACB","ASW"]:
    all_lengths_counts =[]
    for i in range(22):
        print(i)
        with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_gnomix_maxanc_"+str(POPname)+"_chr_"+str(i+1), 'rb') as f:
            lengths_counts=pickle.load(f)
        all_lengths_counts.append(lengths_counts)
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_gnomix_maxanc_"+str(POPname)+"_chr_X", 'rb') as f:
        all_lengths_counts_X=pickle.load(f)
    Aut_matrix=np.apply_along_axis(sum,0,all_lengths_counts)
    X_matrix=all_lengths_counts_X
    Aut_matrix_pop=np.apply_along_axis(np.mean,1,Aut_matrix)
    X_matrix_pop=np.apply_along_axis(np.mean,1,X_matrix)
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_gnomix_maxanc_Aut_pop_'+str(POPname)+'.txt', np.round(Aut_matrix_pop,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_gnomix_maxanc_Aut_ind_'+str(POPname)+'.txt', np.round(Aut_matrix,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_gnomix_maxanc_X_pop_'+str(POPname)+'.txt', np.round(X_matrix_pop,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_gnomix_maxanc_X_ind_'+str(POPname)+'.txt', np.round(X_matrix,2), delimiter='\t', fmt='%f')


for POPname in ["ACB","ASW","CLM","MXL","PEL","PUR"]:
    all_length_anc =[]
    for i in range(22):
        print(i)
        with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_gnomix_"+str(POPname)+"_maxanc_chr_"+str(i+1), 'rb') as f:
            length_anc=pickle.load(f)
        all_length_anc.append(length_anc)
    
    all_length_anc_sum=np.apply_along_axis(sum, 0,all_length_anc)
       
    #all_length_anc_prop=np.apply_along_axis(rel_prop,1,all_length_anc_sum)
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_anc_prop_ind_gnomix_maxanc_Aut_'+str(POPname)+'_bp.txt', np.round(all_length_anc_sum), delimiter='\t', fmt='%f')
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_gnomix_maxanc_"+str(POPname)+"_chr_X", 'rb') as f:
        all_length_anc_X=pickle.load(f)
    
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_anc_prop_ind_gnomix_maxanc_X_'+str(POPname)+'_bp.txt', np.round(all_length_anc_X), delimiter='\t', fmt='%f')
   
