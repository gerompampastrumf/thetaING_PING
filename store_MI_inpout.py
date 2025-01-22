

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
#Calculate the avg modulation inxed over all_input2. Input3 and 4 are used in order to load the file.
from signal_analysis import continufySpikes2
from sklearn.feature_selection import mutual_info_regression

def get_MI_lags(noise,spikes_pyr,all_inds,verbose=False):
    all_iseeds=np.unique(spikes_pyr["iseed"].values)
    all_mi=np.zeros((len(all_iseeds),len(all_inds)))
    for ind_seed,iseed in enumerate(all_iseeds):
        if(verbose):
            print("iseed",int(iseed))
        noise_temp = noise[noise["iseed"]==iseed]["rate"].values
        spikes_pyr_temp = spikes_pyr[(spikes_pyr["iseed"]==iseed)]["rate"].values
        #print(noise_temp[100:110],spikes_pyr_temp[100:110])
        for ind_ir,ir in enumerate(all_inds):
            noise_mov = np.array( list(noise_temp[ir:])+list(noise_temp[:ir]) )
            print(len(noise_mov),len(spikes_pyr_temp))
            all_mi[ind_seed,ind_ir]=mutual_info_regression(spikes_pyr_temp.reshape(-1, 1),noise_mov)       
    return all_mi

def get_autocorr(t1,t2,all_inds):
    t2 = np.tile(np.reshape(t2, (-1, 1)),(1,len(all_inds)))
    shifts = all_inds
    shifted_array = np.empty_like(t2)
    for i, shift in enumerate(shifts):
        shifted_array[:,i] = np.roll(t2[:,i], shift)
    t1=t1[len(all_inds)//2:-len(all_inds)//2]
    shifted_array=shifted_array[len(all_inds)//2:-len(all_inds)//2]
    shifted_array=(shifted_array-shifted_array.mean())/shifted_array.std()
    t1=(t1-t1.mean())/t1.std()
    return shifts,np.dot(t1,shifted_array/len(t1))

def get_corr_lags(noise,spikes_pyr,all_inds,verbose=False):
    all_iseeds=np.unique(spikes_pyr["iseed"].values)
    all_mi=np.zeros((len(all_iseeds),len(all_inds)))
    for ind_seed,iseed in enumerate(all_iseeds):
        if(verbose):
            print("iseed",int(iseed))
        noise_temp = noise[noise["iseed"]==iseed]["rate"].values
        spikes_pyr_temp = spikes_pyr[(spikes_pyr["iseed"]==iseed)]["rate"].values
        #print(noise_temp[100:110],spikes_pyr_temp[100:110])
        all_mi[ind_seed,:] = get_autocorr(noise_temp,spikes_pyr_temp,all_inds)[1]
    return all_mi

current_folder = os.getcwd()

input1 = int(sys.argv[1])
input2 = sys.argv[2]
isInputEC = bool(int(sys.argv[3]))
isMI = bool(int(sys.argv[4]))
folder = os.path.join(current_folder, input2)
print(isInputEC,sys.argv[3])
cells="/external_inputs/external_inputs_pyr"
label_cells = "pyr"
area="ca3"
w = [1,2,5,53][input1]
fs = 500*2
all_inds = np.arange(int(-fs/1000*80),int(fs/1000*80),1)
mode = ["gaussian","exp"][1]
#all_inds = np.arange(-5,5,1)

input_folder = folder+cells+"_ca3.lzma"
spikes_pyr = file_management.load_lzma(input_folder)
tmax = np.max(spikes_pyr["tvec"].values)
if(not isInputEC):
    input_folder = folder+"/external_inputs/external_inputs_DGnoise.lzma"
    print("HERE")
else:
    input_folder = "sigma0_2b/external_inputs4/external_inputs_theta_gen.lzma"
    print("HERE2")
noise_spikes = file_management.load_lzma(input_folder)


noise = continufySpikes2(noise_spikes,fs,"rate",w=w,mode=mode,dt=1000/fs,tmax=tmax)
spikes_pyr = continufySpikes2(spikes_pyr,fs,"rate",w=w,mode=mode,dt=1000/fs,tmax=tmax)
if(isMI):
    all_mi= get_MI_lags(noise,spikes_pyr,all_inds,verbose=True)
    label="MI"
else:
    all_mi= get_corr_lags(noise,spikes_pyr,all_inds,verbose=True)
    label="Corr"

if(not isInputEC):
    title = label+"_inds_w"+str(w)
else:
    title = label+"_inds_w"+str(w)+"_Ext"
if(mode!="gaussian"):
    title=title+"_"+mode
file_management.save_lzma([all_inds,all_mi],title,parent_dir=folder)

os.remove("MI_inds_out"+str(input1)+".dat")
os.remove("MI_inds_err"+str(input1)+".log")