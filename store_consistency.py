

import sys
import os
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/simplifiedCFD/final_files')
import file_management
import time as tm
import mne
from mne import create_info, EpochsArray
import scipy
import pandas as pd
tick = tm.time()
import random
#Calculate the avg modulation inxed over all_input2. Input3 and 4 are used in order to load the file.
from signal_analysis import continufySpikes2
from sklearn.feature_selection import mutual_info_regression

def get_autocorr(t1,t2,N):
    t2 = np.tile(np.reshape(t2, (-1, 1)),(1,2*N))
    shifts = np.arange(-N,N)
    shifted_array = np.empty_like(t2)
    for i, shift in enumerate(shifts):
        shifted_array[:,i] = np.roll(t2[:,i], shift)
    t1=t1[N:-N]
    shifted_array=shifted_array[N:-N]
    shifted_array=(shifted_array-shifted_array.mean())/shifted_array.std()
    t1=(t1-t1.mean())/t1.std()
    return shifts,np.dot(t1,shifted_array/len(t1))

def getConsistence(volt,label,N,type="corr"):
    cors=[]
    cors_control=[]
    inds=np.arange(-N,N)
    print(np.unique(volt["iseed"].values))
    for i in np.unique(volt["iseed"].values):
        for j in range(i):
            v1 = volt[volt["iseed"]==i][label].values
            v1 = (v1-v1.mean())/np.linalg.norm(v1)
            v2 = volt[volt["iseed"]==j][label].values
            v2 = (v2-v2.mean())/np.linalg.norm(v2)
            #ind_control = random.randint(1, len(v2)-1)
            ind_control = random.randint(int(len(v2)/3), int(2*len(v2)/3))
            print(ind_control)
            v1control = np.array( list(v1[ind_control:])+list(v1[:ind_control]) )
            if(type=="corr"):
                cors.append(get_autocorr(v1,v2,N)[1])
                cors_control.append(get_autocorr(v1control,v2,1)[1][-1])
            elif(type=="MI"):
                all_mi=[]
                for ind in inds:
                    v1_moved = np.array( list(v1[ind:])+list(v1[:ind]) )
                    all_mi.append(mutual_info_regression(v2.reshape(-1, 1),v1_moved))
                cors.append(np.array(all_mi).flatten())
                cors_control.append(mutual_info_regression(v2.reshape(-1, 1),v1control)[0])
    return inds,cors,cors_control

current_folder = os.getcwd()

input1 = int(sys.argv[1])
input2 = int(sys.argv[2])
input3 = sys.argv[3]
folder = os.path.join(current_folder, input3)

method = ["corr","MI"][input1]
N=[5,50][input2]
fs = 500*2
#all_inds = np.arange(-5,5,1)

#input_folder = folder+"/volt_pyr_ca3.lzma"
#volt = file_management.load_lzma(input_folder)
#inds,corrs,corrs_err = getConsistence(volt,"soma_volt_mean",N,type=method)
input_folder = folder+"/external_inputs/external_inputs_pyr_ca3.lzma"
spikes_pyr = file_management.load_lzma(input_folder)

w=5
pyrs = continufySpikes2(spikes_pyr,fs,"rate",w=w,mode="exp",dt=1000/fs,tmax = np.max(spikes_pyr["tvec"].values))
inds,cors,cors_control = getConsistence(pyrs,"rate",N,type=method)

title = "Consistence2_method"+method+"_N"+str(N)
file_management.save_lzma([inds,cors,cors_control],title,parent_dir=folder)

os.remove("Consistence_out"+str(input1)+"_"+str(input2)+".dat")
os.remove("Consistence_err"+str(input1)+"_"+str(input2)+".log")