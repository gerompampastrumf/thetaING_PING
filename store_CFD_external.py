import sys
import os
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/simplifiedCFD/final_files')
import file_management
import time as tm
from signal_analysis import compute_coherence_measurements
from signal_analysis import compute_coherence_measurements_diff_high_low
from signal_analysis import get_significant
from signal_analysis import continufySpikes2
import pandas as pd
tick = tm.time()
sys.setrecursionlimit(10000)
#Calculate the avg modulation inxed over all_input2. Input3 and 4 are used in order to load the file.
def get_avg_cfd(fgamma_lims,ftheta_lims,fs,df1,label1,df2,label2,surrogate_test=False,nperseg=2048,testing="cluster",flag_z=False,thres=0.005):  
    all_psi = []
    all_nu = []
    ftheta_min,ftheta_max=ftheta_lims
    fgamma_min,fgamma_max,dfgamma=fgamma_lims
    all_inputs1=np.unique(df1["iseed"].values)
    n_temps=0
    for input1 in all_inputs1:
        print("sim: "+str(input1))
        n_temps+=1
        ts1 = df1[label1][(df1['iseed'] == input1)].values
        ts1 = np.array(ts1).astype(np.float32)
        ts1=ts1[len(ts1)//15:]
        ts2 = df2[label2][(df2['iseed'] == input1)].values
        ts2 = np.array(ts2).astype(np.float32)
        ts2=ts2[len(ts2)//15:]
        if(flag_z):
            ts1-=np.mean(ts1)
            ts1/=np.std(ts1)
            ts2-=np.mean(ts2)
            ts2/=np.std(ts2)
        ftheta, fgamma, psi_temp, nu_temp, psi_sur_temp,nu_sur_temp = compute_coherence_measurements_diff_high_low(ts1,ts2, fs, nperseg, 0.5*nperseg, ftheta_min=ftheta_min, 
                               ftheta_max=ftheta_max, fgamma_min=fgamma_min, fgamma_max=fgamma_max, dfgamma=dfgamma,
                               norder=4, nfft=nperseg, surrogate_test=surrogate_test, nsurrogates=1000,choiseSurr = 2,ibeta=4)
        
        psi_temp= np.array(psi_temp)
        nu_temp=np.array(nu_temp)
        
        if(surrogate_test):
            psi_sur_temp=np.array(psi_sur_temp)
            mask = get_significant(psi_temp,psi_sur_temp,thres1 = thres,thres2 = thres,testing=testing)
            psi_temp = psi_temp*mask
        
        all_psi.append(psi_temp)
        all_nu.append(nu_temp)
    return ftheta,fgamma,all_psi,all_nu


input_folder = "sigma0_2b/external_inputs4/external_inputs_theta_gen.lzma"
theta_spikes = file_management.load_lzma(input_folder)
fs = 500*2
current_folder = os.getcwd()

input2 = sys.argv[1]
input3 = int(sys.argv[2])
folder = os.path.join(current_folder, input2)
flag_z = True
    
area="ca3"
surrogate_test=True
testing=["cluster","point"][0]
thres = [0.01,0.005,0.05][input3]
w=5

for input1 in range(6):
    ts_choise = ["ica","volt_pyr","external_inputs/external_inputs_pyr","volt_bas","lfp","synapses_pyr"][input1]
    
    df= file_management.load_lzma(os.path.join(folder,ts_choise+"_"+area+".lzma"))
    ised = np.unique(df["iseed"].values)[0]

    if(ts_choise=="lfp"):
        elec_pos = np.unique(df["electrode"].values)[0]
        ts_ = df[ts_choise][(df["electrode"]==elec_pos)&(df["iseed"]==ised)].values
        pos_vals = np.unique(df["electrode"].values)
    elif(ts_choise=="ica"):
        elec_comp=np.unique(df["component"].values)[0]
        elec_pos = np.unique(df["electrode"].values)[0]
        ts_ = df[ts_choise][(df["electrode"]==elec_pos)&(df["component"]==elec_comp)&(df["iseed"]==ised)].values
        #pos_vals = ["Bdend","soma","Adend1","Adend2","Adend3"]
        pos_vals = ["soma"]
    elif(ts_choise=="synapses_pyr"):
        all_comp=["Bdend","soma","Adend1","Adend2","Adend3"]
        df = pd.concat([*[df.filter(like=comp).filter(like="mean").sum(axis=1) for comp in all_comp],*[df["iseed"]]],axis=1,keys=[*all_comp,"iseed"])
        ts_ = df[df["iseed"]==ised][all_comp[0]].values
        df.fillna(0,inplace=True)
        pos_vals = all_comp
    elif("volt" in ts_choise):
        pos_vals=[k for k in list(df.keys()) if "volt_mean" in k]
        all_comp=pos_vals
        df =pd.concat([*[df.filter(like=comp).filter(like="volt_mean").sum(axis=1) for comp in all_comp],*[df["iseed"]]],axis=1,keys=[*all_comp,"iseed"])
        ts_ = df[pos_vals[0]][(df["iseed"]==ised)].values
    elif("external_inputs_pyr" in ts_choise):
        df = continufySpikes2(df,fs,"rates_pyr",w=w,mode="exp",dt=1000/fs,tmax=np.max(df["tvec"].values))
        ts_ = df["rates_pyr"][(df["iseed"]==ised)].values
        pos_vals=["none"]

    theta_inp = continufySpikes2(theta_spikes,fs,"theta_inp",w=w,mode="exp",dt=1000/fs,tmax=len(ts_)/fs*1000)

    ftheta_lims = 2,15
    fgamma_lims = 20,120,1
    nperseg = int(fs/500*1024)

    all_psi= {}
    all_psi_std={}
    all_nu = {}
    all_nu_std={}

    for ind,pos in enumerate(pos_vals):
        print(pos)
        if(ts_choise=="lfp"):
            df_pos = df[df["electrode"]==pos]
        elif(ts_choise=="ica"):
            df_pos = df[(df["electrode"]==elec_pos)&(df["component"]==pos)]
        elif(ts_choise=="synapses_pyr"):
            df_pos = df[[pos,"iseed"]]
            df_pos.rename(columns={pos:ts_choise},inplace=True)
        elif("volt" in ts_choise):
            df_pos = df[[pos,"iseed"]]
            df_pos.rename(columns={pos:ts_choise},inplace=True)
        elif("external_inputs_pyr"in ts_choise):
            ts_choise = "rates_pyr"
            df_pos = df.copy()
        ftheta,fgamma,psi_tmp,nu_tmp= get_avg_cfd(fgamma_lims,ftheta_lims,fs,theta_inp,"theta_inp",df_pos,ts_choise,surrogate_test=surrogate_test,nperseg=nperseg,testing=testing,flag_z=flag_z,thres=thres)
        all_psi[pos]=np.mean(psi_tmp,axis=0)
        all_psi_std[pos]=np.std(psi_tmp,axis=0)
        all_nu[pos]=np.mean(nu_tmp,axis=0)
        all_nu_std[pos]=np.std(nu_tmp,axis=0)
        

    if(surrogate_test):
        title = "CFD"+str([2,4,5][input3])+"Ext_fft"+str(nperseg)+"_z"+str(flag_z)+"_sign"+testing+"_"+ts_choise+"_"+area
    else:
        title = "CFDExt_fft"+str(nperseg)+"_z"+str(flag_z)+"_"+ts_choise+"_"+area

    file_management.save_lzma([all_psi,all_psi_std,all_nu,all_nu_std,fgamma,ftheta],title,parent_dir=folder)

    file_name = __file__
    store = 'cp '+__file__.replace(current_folder+"/","") +" "+ folder.replace(current_folder+"/","") +__file__.replace(current_folder,"")

tock = tm.time()
diff = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

print("All data saved after ", horas,'h', mint, 'min', seg, 's')
os.popen(store)

os.remove("CFD_out"+".dat")
os.remove("CFD_err"+".log")