

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
cells=["/external_inputs/external_inputs_pyr","/spikes_bas"][input3]
label_cells = ["pyr","bas"][input3]
area="ca3"


fs = 500*2
input_folder = folder+cells+"_ca3.lzma"
spikes_pyr = file_management.load_lzma(input_folder)
input_folder = folder+"/synapses_pyr_ca3.lzma"
syn_pyr = file_management.load_lzma(input_folder)
syn = pd.concat([-syn_pyr.filter(like="ec2180_mean").filter(like="iAdend3").sum(axis=1),syn_pyr["iseed"]],
                  axis=1,keys=["syn_ec2","iseed"])

all_phase=get_theta_spiking(spikes_pyr,syn,"syn_ec2",fs=fs,f0=8.0,df0=df0)

title = "thetaSpikingNew_"+label_cells+"_df"+str(df0)

file_management.save_lzma([all_phase],title,parent_dir=folder)

os.remove("thetaSpikingNew_out"+str(input1)+"_"+str(input2)+"_"+str(input3)+".dat")
os.remove("thetaSpikingNew_err"+str(input1)+"_"+str(input2)+"_"+str(input3)+".log")