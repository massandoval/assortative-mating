#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import scipy
#import sklearn
#from sklearn import metrics
#from sklearn.model_selection import train_test_split
#from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import os
#import plotnine as pn
#import pydot
#import pydotplus
import random
from math import log
#from sklearn.preprocessing import MinMaxScaler
#from plotnine import ggplot, stat_boxplot, ggtitle, aes, stat_function, theme, theme_minimal, scale_x_discrete, element_text, position_dodge, xlab, ylab, geom_bar, geom_col,geom_point, geom_line, stat_summary, scale_y_continuous,scale_x_continuous, scale_color_manual,scale_fill_manual,geom_boxplot,facet_wrap, facet_grid,geom_point,xlim,ylim,geom_bin2d,scale_fill_continuous,geom_density,ggsave


# In[2]:


POPnum = sys.argv[1]
falsetrue=sys.argv[2]
#POPnum = 1
#falsetrue="false"
bins=len(np.arange(1,11.5,0.5))+1	 #fragment windows
ancestries=3 # 
samplePOP=100
TOTALPOP=1000
GEN=19
scale=1000
SCALE=scale
POPname_list=["CLM","MXL","PEL","PUR","ACB","ASW"]
POP1_list=[0.082,0.043,0.032,0.138,0.881,0.765]
POP2_list=[0.277,0.5,0.761,0.15,0.005,0.041]
POPname=POPname_list[int(POPnum)-1]
POP1=POP1_list[int(POPnum)-1]
POP2=POP2_list[int(POPnum)-1]




# In[3]:

runspredictedfile="/rds/general/project/human-popgen-datasets/live/CNN_AM/log_pulses_"+str(falsetrue)+"_newgenmap/runs_predicted_"+str(POPname)+"_oct23.txt"
if (falsetrue=="false"):
    command_runs='ls -lth /rds/general/project/human-popgen-datasets/live/CNN_AM/log_pulses_'+str(falsetrue)+'_newgenmap/'+str(POPname)+'/fragments/Frags_AutX_log_pulses_'+str(falsetrue)+'_newgenmap_predicted_'+str(POPname)+'_AM1_*.txt | grep -v oct2023 | grep "K Oct" | sed "s/.*_AM1_\(.*\)_AM2_\(.*\)_AM3_\(.*\)_SB1_\(.*\)_SB2_\(.*\)_GEN.*run_\(.*\).txt/\\1 \\2 \\3 \\4 \\5 \\6/g" | sort -n -k6,6 > '+runspredictedfile

if (falsetrue=="true"):
    command_runs='ls -lth /rds/general/project/human-popgen-datasets/live/CNN_AM/log_pulses_'+str(falsetrue)+'_newgenmap/'+str(POPname)+'/fragments/Frags_AutX_log_pulses_'+str(falsetrue)+'_newgenmap_predicted_'+str(POPname)+'_AM1_*.txt | grep -v oct2023 | grep "K Oct" | sed "s/.*_AM1_\(.*\)_AM2_\(.*\)_AM3_\(.*\)_SB1_\(.*\)_SB2_\(.*\)_GEN.*run_\(.*\).txt/\\1 \\2 \\3 \\4 \\5 \\6/g" | sort -n -k6,6 > '+runspredictedfile

os.system(command_runs)
AMcomb=np.loadtxt(runspredictedfile, dtype='str', delimiter=' ')
n=len(AMcomb)
POP3=1-(POP1+POP2)


# In[4]:


x= ["NA" for i in range(n)]
count=0
count2=0
runvec = []
AM1vec = []
AM2vec = []
AM3vec = []
SB1vec = []
SB2vec = []
#prop_pop1_t1_vec = []
#prop_pop2_t1_vec = []
#prop_pop3_t1_vec = []
#prop_pop1_t2_vec = []
#prop_pop2_t2_vec = []
#prop_pop3_t2_vec = []
for AMFMSB in AMcomb:    
    AM1 = AMFMSB[0]
    AM2 = AMFMSB[1]
    AM3 = AMFMSB[2]
    SB1 = AMFMSB[3]
    SB2 = AMFMSB[4]
    run = AMFMSB[5]  
   #print(run)
    input_frag ='/rds/general/project/human-popgen-datasets/live/CNN_AM/log_pulses_'+str(falsetrue)+'_newgenmap/'+str(POPname)+'/fragments/Frags_AutX_log_pulses_'+str(falsetrue)+'_newgenmap_predicted_'+str(POPname)+'_AM1_'+str(AM1)+'_AM2_'+str(AM2)+'_AM3_'+str(AM3)+'_SB1_'+str(SB1)+'_SB2_'+str(SB2)+'_GEN_19_POP1_'+str(round(POP1*1000))+'_POP2_'+str(round(POP2*1000))+'_POP3_'+str(round(POP3*1000))+'_scale_1000_run_'+str(run)+'.txt'
    frag_lengths = np.loadtxt(input_frag, dtype='str', delimiter='; ')
    frag_lengths_2 = [[int(cell.replace(';', '')) for cell in row] for row in frag_lengths]
    frag_lengths_3 = np.reshape(frag_lengths_2,(samplePOP,2,ancestries,bins))
     ##sort individuals. first by max ancestry, then by largest mean length fragment of this ancestry.
    maxanc_vec = []
    meanlen_maxanc_vec = []
    prop_maxanc_vec =[]
    for ind in range(samplePOP):
        #n_frags=np.sum(frag_lengths_3[ind][0],axis=1)
        total_lengths=np.sum(frag_lengths_3[ind][0]*range(bins),axis=1)
        maxanc=np.argmax(total_lengths)
        prop_maxanc=(total_lengths)[maxanc]
        #meanlen_maxanc=(total_lengths/n_frags)[maxanc]
        maxanc_vec.append(maxanc)
        prop_maxanc_vec.append(prop_maxanc)
        #meanlen_maxanc_vec.append(meanlen_maxanc)
    maxanc_meanlen_array = np.array([range(samplePOP),maxanc_vec,prop_maxanc_vec])
    frag_lengths_4=frag_lengths_3[np.lexsort((maxanc_meanlen_array[1,:],maxanc_meanlen_array[2,:]))]
    #print(np.shape(frag_lengths_4))
    frag_lengths_5 = np.swapaxes(frag_lengths_4,1,3)    
    #frag_lengths_5 = frag_lengths_4[random.sample(range(0,len(frag_lengths_4)),len(frag_lengths_4))]    #Shuffle individuals
    runvec.append(float(run))
    AM1vec.append(float(AM1))
    AM2vec.append(float(AM2))     
    AM3vec.append(float(AM3)) 
    SB1vec.append(float(SB1)) 
    SB2vec.append(float(SB2))
    x[count]=frag_lengths_5
    count=count+1
    if (count%100==0):
        print(count)

x=np.array(x)
           


# In[5]:


list_indices=np.array(range(n))
random.shuffle(list_indices)
list_indices=list_indices.astype(int)

import pickle
with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/Bootstrap_"+str(POPname)+"_log_rfmix_maxanc", 'rb') as f:
    bootstrap_Aut= pickle.load(f)

import pickle
with open("/rds/general/user/amassand/home/SLIM_models/MODELS/pops/Bootstrap_"+str(POPname)+"_log_rfmix_maxanc_X_onlyfemales", 'rb') as f:
    bootstrap_X= pickle.load(f)

target_panel=np.loadtxt("/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/target_panel.txt",delimiter='\t',dtype="str")
targetAFRAME_panel=np.loadtxt("/rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/targetAFRAME_panel.txt",delimiter='\t',dtype="str")
all_inds=np.concatenate([target_panel,targetAFRAME_panel],axis=0)
inds_1=np.where(all_inds[:,1]==str(POPname))[0]*6
inds_2=np.where(all_inds[:,1]==str(POPname))[0]*6+1
inds_3=np.where(all_inds[:,1]==str(POPname))[0]*6+2
inds_4=np.where(all_inds[:,1]==str(POPname))[0]*6+3
inds_5=np.where(all_inds[:,1]==str(POPname))[0]*6+4
inds_6=np.where(all_inds[:,1]==str(POPname))[0]*6+5
indices=np.sort(np.concatenate([inds_1,inds_2,inds_3,inds_4,inds_5,inds_6]))

females_1=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F"))[0]*6
females_2=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F"))[0]*6+1
females_3=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F"))[0]*6+2
females_4=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F"))[0]*6+3
females_5=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F"))[0]*6+4
females_6=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F"))[0]*6+5
males_1=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="M"))[0]*6
males_2=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="M"))[0]*6+1
males_3=np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="M"))[0]*6+2
indices_X=np.sort(np.concatenate([females_1,females_2,females_3,females_4,females_5,females_6]))

n_females=np.shape(np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="F")))[1]
n_males=np.shape(np.where((all_inds[:,1]==str(POPname))&(all_inds[:,2]=="M")))[1]

bootstrap_X=np.array(bootstrap_X[0])
bootstrap_Aut=np.array(bootstrap_Aut[0])

x_real_pop_bootstrap=[]
for i in range(1000):
    bootstrap_pop_X=np.mean(bootstrap_X[i].reshape(22,int(np.shape(indices_X)[0]/3),3),axis=1)
    bootstrap_pop_Aut=np.mean(bootstrap_Aut[i].reshape(22,int(np.shape(indices)[0]/3),3),axis=1)
    x_real_pop=np.stack((bootstrap_pop_Aut,bootstrap_pop_X),axis=-1)
    x_real_pop_bootstrap.append(x_real_pop)    

x_real_pop_bootstrap=np.array(x_real_pop_bootstrap)
x_real_pop_bootstrap=x_real_pop_bootstrap[:,:,np.array((2,0,1)),:]
x_real_pop_bootstrap_mean=np.mean(x_real_pop_bootstrap[:,:,:,:],axis=0)


x=x[list_indices]
AMcomb=AMcomb[list_indices]
x_pop = np.mean(x.tolist(),axis=1)




# In[130]:


###Data frame to plot real fragments
ancestry=np.repeat([[1,2,3],],1000*22,axis=0).reshape(22*3*1000)
frag_length=np.repeat([np.repeat([np.array(range(22)),],3,axis=1).reshape(22*3)],1000,axis=0).reshape(22*3*1000)
counts_real_Aut_bootstrap=np.reshape([bootstrap_run.reshape(22*3) for bootstrap_run in x_real_pop_bootstrap[:,:,:,0]],1000*22*3)
counts_real_X_bootstrap=np.reshape([bootstrap_run.reshape(22*3) for bootstrap_run in x_real_pop_bootstrap[:,:,:,1]],1000*22*3)  
Aut_label=np.repeat("Aut",1000*22*3)
X_label=np.repeat("X",1000*22*3)
Bootstrap_label=np.repeat(range(1,1001),22*3)
dataframe_Aut_real=pd.DataFrame(data=[counts_real_Aut_bootstrap, ancestry,frag_length,Aut_label,Bootstrap_label]).T
dataframe_X_real=pd.DataFrame(data=[counts_real_X_bootstrap, ancestry,frag_length,X_label,Bootstrap_label]).T
df_real_AutX=pd.concat([dataframe_Aut_real,dataframe_X_real])
df_real_AutX.columns=['counts', 'ancestry', 'frag_length','AutX','Bootstrap']
df_real_AutX = df_real_AutX.astype({"counts":float, "frag_length": int, "AutX": "category","Bootstrap": int})

ancestry=np.repeat([[1,2,3],],22*np.shape(x_pop)[0],axis=0).reshape(22*3*np.shape(x_pop)[0])
frag_length=np.repeat([np.repeat([np.array(range(22)),],3,axis=1).reshape(22*3),],np.shape(x_pop)[0],axis=0).reshape(22*3*np.shape(x_pop)[0])
counts_sim_Aut=x_pop[:,:,:,0].reshape(22*3*np.shape(x_pop)[0])
counts_sim_X=x_pop[:,:,:,1].reshape(22*3*np.shape(x_pop)[0])
Aut_label=np.repeat("Aut",22*3*np.shape(x_pop)[0])
X_label=np.repeat("X",22*3*np.shape(x_pop)[0])
dataframe_Aut_sim=pd.DataFrame(data=[counts_sim_Aut, ancestry,frag_length,Aut_label]).T
dataframe_X_sim=pd.DataFrame(data=[counts_sim_X, ancestry,frag_length,X_label]).T
df_sim_AutX=pd.concat([dataframe_Aut_sim,dataframe_X_sim])
df_sim_AutX.columns=['counts', 'ancestry', 'frag_length','AutX']
df_sim_AutX = df_sim_AutX.astype({"counts":float, "frag_length": "category", "AutX": "category"})

df_sim_Aut=df_sim_AutX[df_sim_AutX["AutX"]=="Aut"]
df_sim_X=df_sim_AutX[df_sim_AutX["AutX"]=="X"]
df_real_Aut=df_real_AutX[df_real_AutX["AutX"]=="Aut"]
df_real_X=df_real_AutX[df_real_AutX["AutX"]=="X"]


# In[132]:


import scipy.stats
df_sim_X_summary=df_sim_X.groupby(['frag_length','ancestry'])["counts"].describe()
log_probs_X_cell=[np.array(np.log(scipy.stats.poisson.pmf(np.round(np.array(df_real_X["counts"][df_real_X["Bootstrap"]==bootstrap])),np.array(df_sim_X_summary["mean"])))) 
 for bootstrap in np.array(range(1,1001))]
log_probs_X=np.apply_along_axis(sum,1,log_probs_X_cell)
df_sim_Aut_summary=df_sim_Aut.groupby(['frag_length','ancestry'])["counts"].describe()
log_probs_Aut_cell=[np.array(np.log(scipy.stats.poisson.pmf(np.round(np.array(df_real_Aut["counts"][df_real_Aut["Bootstrap"]==bootstrap])),np.array(df_sim_Aut_summary["mean"])))) 
 for bootstrap in np.array(range(1,1001))]
log_probs_Aut=np.apply_along_axis(sum,1,log_probs_Aut_cell)


# In[169]:


df_log=pd.DataFrame({"frag_length":np.repeat(np.array([[item[0] for item in np.array(df_sim_X_summary.index)],]),1000,axis=0).reshape(1000*22*3),
                     "ancestry":np.repeat(np.array([[item[1] for item in np.array(df_sim_X_summary.index)],]),1000,axis=0).reshape(1000*22*3),                     
                     "Aut":np.reshape(log_probs_Aut_cell,22*3*1000),
                     "X":np.reshape(log_probs_X_cell,22*3*1000),
                     "pop":np.repeat(POPname,22*3*1000),
                     "bootstrap":np.repeat(range(1,1001),22*3),
                     "pulses":np.repeat(falsetrue,22*3*1000)})
log_prob_file="/rds/general/project/human-popgen-datasets/live/CNN_AM/log_prob_poisson_file_oct23.csv"
df_log.to_csv(log_prob_file, mode='a', header=False,  index=False, sep=",")





