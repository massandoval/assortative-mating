
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
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_rfmix_maxanc_"+str(POPname)+"_chr_X_onlyfemales", 'rb') as f:
        all_lengths_counts_X=pickle.load(f)
    X_matrix=all_lengths_counts_X
    X_matrix_pop=np.apply_along_axis(np.mean,1,X_matrix)
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_rfmix_maxanc_X_pop_'+str(POPname)+'_onlyfemales.txt', np.round(X_matrix_pop,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_rfmix_maxanc_X_ind_'+str(POPname)+'_onlyfemales.txt', np.round(X_matrix,2), delimiter='\t', fmt='%f')


for POPname in ["ACB","ASW","CLM","MXL","PEL","PUR"]:
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_rfmix_"+str(POPname)+"_maxanc_chr_X_onlyfemales", 'rb') as f:
        all_length_anc_X=pickle.load(f)
    
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_anc_prop_ind_rfmix_maxanc_X_'+str(POPname)+'_onlyfemales_bp.txt', np.round(all_length_anc_X), delimiter='\t', fmt='%f')
   
for POPname in ["CLM","MXL","PEL","PUR","ACB","ASW"]:
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/lengths_counts_gnomix_maxanc_"+str(POPname)+"_chr_X_onlyfemales", 'rb') as f:
        all_lengths_counts_X=pickle.load(f)
    X_matrix=all_lengths_counts_X
    X_matrix_pop=np.apply_along_axis(np.mean,1,X_matrix)
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_gnomix_maxanc_X_pop_'+str(POPname)+'_onlyfemales.txt', np.round(X_matrix_pop,2), delimiter='\t', fmt='%f')
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_profile_gnomix_maxanc_X_ind_'+str(POPname)+'_onlyfemales.txt', np.round(X_matrix,2), delimiter='\t', fmt='%f')


for POPname in ["ACB","ASW","CLM","MXL","PEL","PUR"]:
    with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/length_anc_gnomix_"+str(POPname)+"_maxanc_chr_X_onlyfemales", 'rb') as f:
        all_length_anc_X=pickle.load(f)
    
    np.savetxt('/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/lengths_anc_prop_ind_gnomix_maxanc_X_'+str(POPname)+'_onlyfemales_bp.txt', np.round(all_length_anc_X), delimiter='\t', fmt='%f')
   
