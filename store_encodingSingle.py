from sklearn.feature_selection import mutual_info_regression
import sys
import os
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/simplifiedCFD/final_files')
import file_management
import time as tm
from signal_analysis import continufySpikes2

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

def get_encoding(rates,rates_n,noise,all_inds,isMI,Nseg,Tmax):
    mis_per_spike={}
    spikes_counted={}
    for iseed in np.unique(rates["iseed"].values):
        mis_per_spike[iseed]=[]
        spikes_counted[iseed]=[]
        all_t = noise[noise["iseed"]==iseed]["tvec"].values
        all_t = all_t[(all_t<Tmax-Nseg*all_inds[-1])&(all_t>-Nseg*all_inds[0])]
        for t in all_t:
            inds = (int(t)+Nseg*all_inds[0]-1,int(t)+Nseg*all_inds[-1]+1)
            rat_n = rates_n[rates_n["iseed"]==iseed]["rate"].values[inds[0]:inds[-1]]
            rat = rates[rates["iseed"]==iseed]["rate"].values[inds[0]:inds[-1]]
            all_mi=[]
            print(t,inds,len(rat_n),len(rat))
            if(isMI):
                for _,ir in enumerate(all_inds):
                    noise_mov = np.array( list(rat_n[ir:])+list(rat_n[:ir]) )
                    print(len(rat),len(noise_mov))
                    all_mi.append(mutual_info_regression(rat.reshape(-1, 1),noise_mov))
            else:
                _,all_mi=get_autocorr(rat_n,rat,all_inds)
            mis_per_spike[iseed].append(all_mi)
            spikes_counted[iseed].append(t)
        mis_per_spike[iseed]=np.array(mis_per_spike[iseed])
    return spikes_counted,mis_per_spike

input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = int(sys.argv[3])
input4 = int(sys.argv[4])
input5 = bool(int(sys.argv[5]))
input6 = int(sys.argv[6])

w=input3
mode=["exp","gaussian"][input4]
flag_MI = input5
Nseg = [1,2,4,8,16,32,64][input6]

fs=1000
all_inds = np.arange(-62,63,1)


spikes = file_management.load_lzma(input1+"/external_inputs/external_inputs_pyr_ca3.lzma")
noise = file_management.load_lzma(input2)
noise=noise[noise["tvec"].values<np.max(spikes["tvec"].values)]
print(noise)
rates_n = continufySpikes2(noise,fs,"rate",w=w,mode=mode,dt=1,tmax = np.max(spikes["tvec"].values))
rates = continufySpikes2(spikes,fs,"rate",w=w,mode=mode,dt=1,tmax = np.max(spikes["tvec"].values))
spikes_counted,mis_per_spike=get_encoding(rates,rates_n,noise,all_inds,flag_MI,Nseg,Tmax=np.max(spikes["tvec"].values))

title = "Encoding_MI"+str(flag_MI)+"_mode"+mode+"_w"+str(w)+"_N"+str(Nseg)
file_management.save_lzma([all_inds,spikes_counted,mis_per_spike],title,parent_dir=input1)

#os.remove("Encode_out"+str(input3)+"_"+str(input4)+"_"+str(input5)+"_"+str(input6)+".dat")
#os.remove("Encode_err"+str(input3)+"_"+str(input4)+"_"+str(input5)+"_"+str(input6)+".log")