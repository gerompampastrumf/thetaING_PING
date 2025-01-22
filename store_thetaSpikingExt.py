

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
from signal_analysis import bandpass_filter_and_hilbert_transform2
from signal_analysis import bound_timenphases
from signal_analysis import continufySpikes2

def get_theta_spiking(spikes,timeseries,label,fs=500,f0=8.0,df0=2.0,norder=4,offset=0):
    all_phase=[]
    for iseed in np.unique(timeseries["iseed"].values):
        ts = timeseries[(timeseries["iseed"]==iseed)][label].values
        ts = ts -ts.mean()
        off = int(np.ceil(offset*fs))
        ts = ts[off:]
        spikes_temp = spikes[(spikes["iseed"]==iseed)]["tvec"].values/1e3
        _, _,  phase = bandpass_filter_and_hilbert_transform2(ts-ts.mean(), fs, f0, df0, norder)
        time=np.linspace(0,len(ts)*1/fs,len(ts)) 
        time,phase=bound_timenphases(time,phase)      
        all_phase.append(np.interp(spikes_temp, time, phase, left=None, right=None, period=None))
    return all_phase

current_folder = os.getcwd()

input1 = int(sys.argv[1])
input2 = int(sys.argv[2])
input3 = int(sys.argv[3])
input4 = sys.argv[4]
folder = os.path.join(current_folder, input4)
df0 = [2,4,6][input1]
#method = ["Hilbert","bycycle"][input2]
w=[1,5][input2]
cells=["/external_inputs/external_inputs_pyr","/spikes_bas"][input3]
label_cells = ["pyr","bas"][input3]
area="ca3"


fs = 500*2
input_folder = folder+cells+"_ca3.lzma"
spikes_pyr = file_management.load_lzma(input_folder)
input_folder = "sigma0_2b/external_inputs4/external_inputs_theta_gen.lzma"
theta_spikes = file_management.load_lzma(input_folder)
mode = "exp"
theta_inp = continufySpikes2(theta_spikes,fs,"theta_inp",w=w,mode=mode,dt=1000/fs,tmax=np.max(spikes_pyr["tvec"].values))

all_phase=get_theta_spiking(spikes_pyr,theta_inp,"theta_inp",fs=fs,f0=8.0,df0=df0)

title = "thetaSpikingExt_"+label_cells+"_w"+str(w)+"_mode"+mode+"_df"+str(df0)

file_management.save_lzma([all_phase],title,parent_dir=folder)

os.remove("thetaSpikingExt_out"+str(input1)+"_"+str(input2)+"_"+str(input3)+".dat")
os.remove("thetaSpikingExt_err"+str(input1)+"_"+str(input2)+"_"+str(input3)+".log")