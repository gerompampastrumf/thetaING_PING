

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
from signal_analysis import bandpass_filter_and_hilbert_transform2

def split_consecutive_series(ts):
    if len(ts) == 0:
        return []
    ts = np.array(ts)
    # Compute the differences between consecutive elements
    diffs = np.diff(ts)
    # Identify the indices where the difference is greater than 1
    gap_indices = np.where(diffs > 1)[0] + 1
    # Split the array at the gap indices
    sublists = np.split(ts, gap_indices)
    # Convert sublists from numpy arrays to lists
    result = [sublist.tolist() for sublist in sublists]
    return result

def get_MI_lags(noise,spikes_pyr,all_inds,theta_ts,theta_bins,verbose=False,fs=1000, f0=8.0, df0=2, norder=4):
    all_iseeds=np.unique(spikes_pyr["iseed"].values)
    all_mi=np.zeros((len(all_iseeds),len(theta_bins)-1,len(all_inds)))
    for ind_seed,iseed in enumerate(all_iseeds):
        if(verbose):
            print("iseed",int(iseed))
        noise_temp = noise[noise["iseed"]==iseed]["rate"].values
        spikes_pyr_temp = spikes_pyr[(spikes_pyr["iseed"]==iseed)]["rate"].values
        ref = theta_ts[(theta_ts["iseed"]==iseed)]["soma_volt_mean"].values
        _, _,  phase = bandpass_filter_and_hilbert_transform2(ref-ref.mean(), fs, f0, df0, norder)
        
        for ind_theta in range(len(theta_bins[:-1])):
            inds_phase = (phase>=theta_bins[ind_theta]) & (phase<theta_bins[ind_theta+1])
            inds_phase = split_consecutive_series(inds_phase)
            #spikes_pyr_temp2 = spikes_pyr_temp[inds_phase]
            #spikes_pyr_temp2 = np.copy(spikes_pyr_temp)
            #spikes_pyr_temp2[~inds_phase]=0.0
            
            for ind_ir,ir in enumerate(all_inds):
                noise_mov2 = np.array( list(noise_temp[ir:])+list(noise_temp[:ir]) )
                
                #noise_mov = noise_mov2[inds_phase]
                #noise_mov = np.copy(noise_mov2)
                #noise_mov[~inds_phase]=0.0
                all_mi_temp = []
                for ind_phase in inds_phase:
                    noise_mov=noise_mov2[ind_phase]
                    spikes_pyr_temp2=spikes_pyr_temp[ind_phase]
                    all_mi_temp.append(mutual_info_regression(spikes_pyr_temp2.reshape(-1, 1),noise_mov))
                all_mi[ind_seed,ind_theta,ind_ir]=np.mean(all_mi_temp)
                #all_mi[ind_seed,ind_theta,ind_ir]=mutual_info_regression(spikes_pyr_temp2.reshape(-1, 1),noise_mov)       
    return all_mi

current_folder = os.getcwd()

input1 = int(sys.argv[1])
input2 = sys.argv[2]
isInputEC = bool(int(sys.argv[3]))
folder = os.path.join(current_folder, input2)
print(isInputEC,sys.argv[3])
cells="/external_inputs/external_inputs_pyr"
label_cells = "pyr"
area="ca3"
w = [1,2,5][input1]
fs = 500*2
all_inds = np.arange(int(-fs/1000*80),int(fs/1000*80),1)
theta_bins = np.linspace(0,2*np.pi,5)

input_folder = folder+cells+"_ca3.lzma"
spikes_pyr = file_management.load_lzma(input_folder)
if(not isInputEC):
    input_folder = folder+"/external_inputs/external_inputs_DGnoise.lzma"
    print("HERE")
else:
    input_folder = "sigma0_2b/external_inputs4/external_inputs_theta_gen.lzma"
    print("HERE2")
noise_spikes = file_management.load_lzma(input_folder)

volt = file_management.load_lzma(folder+"/volt_pyr_ca3.lzma")
tmax = len(volt[volt["iseed"]==0]["soma_volt_mean"].values)/fs*1000+1000/fs
print(tmax)
noise = continufySpikes2(noise_spikes,fs,"rate",w=w,mode="gaussian",dt=1000/fs,tmax=tmax)
spikes_pyr = continufySpikes2(spikes_pyr,fs,"rate",w=w,mode="gaussian",dt=1000/fs,tmax=tmax)
all_mi= get_MI_lags(noise,spikes_pyr,all_inds,volt,theta_bins,verbose=True)
print(len(spikes_pyr))
if(not isInputEC):
    title = "MIb_theta_inds_w"+str(w)
else:
    title = "MIb_theta_inds_w"+str(w)+"_Ext"
file_management.save_lzma([all_inds,theta_bins,all_mi],title,parent_dir=folder)

os.remove("MI_inds_out"+str(input1)+".dat")
os.remove("MI_inds_err"+str(input1)+".log")