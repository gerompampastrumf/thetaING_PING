import sys
import os
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/simplifiedCFD/final_files')
import file_management
import time as tm
from signal_analysis import compute_coherence_measurements
from signal_analysis import get_significant
from signal_analysis import continufySpikes2
import pandas as pd
tick = tm.time()
sys.setrecursionlimit(10000)
#Calculate the avg modulation inxed over all_input2. Input3 and 4 are used in order to load the file.
def get_avg_cfd(fgamma_lims,ftheta_lims,fs,df,label,surrogate_test=False,nperseg=2048,testing="cluster",flag_z=False,thres = 0.005): 
    all_psi = []
    all_nu = []
    ftheta_min,ftheta_max=ftheta_lims
    fgamma_min,fgamma_max,dfgamma=fgamma_lims
    all_inputs1=np.unique(df["iseed"].values)
    n_temps=0
    for input1 in all_inputs1:
        print("sim: "+str(input1))
        n_temps+=1
        ts = df[label][(df['iseed'] == input1)].values
        ts=ts[len(ts)//15:]
        if(flag_z):
            ts-=np.mean(ts)
            ts/=np.std(ts)
        ftheta, fgamma, psi_temp, nu_temp, psi_sur_temp,nu_sur_temp = compute_coherence_measurements(ts, fs, nperseg, 0.5*nperseg, ftheta_min=ftheta_min, 
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

#input1 = int(sys.argv[1])
fs = 500*2
current_folder = os.getcwd()
for input1 in range(6):
    ts_choise = ["ica","volt_pyr","external_inputs/external_inputs_pyr","volt_bas","lfp","synapses_pyr"][input1]
    input2 = sys.argv[1]
    input3 = int(sys.argv[2])
    input4 = int(sys.argv[3])

    folder = os.path.join(current_folder, input2)
    flag_z = True
    area="ca3"
    surrogate_test=False
    testing=["cluster","point"][0]
    df= file_management.load_lzma(os.path.join(folder,ts_choise+"_"+area+".lzma"))
    ised = np.unique(df["iseed"].values)[0]

    if(ts_choise=="lfp"):
        elec_pos = np.unique(df["electrode"].values)[0]
        pos_vals = np.unique(df["electrode"].values)
    elif(ts_choise=="ica"):
        elec_comp=np.unique(df["component"].values)[0]
        elec_pos = np.unique(df["electrode"].values)[0]
        pos_vals = ["Bdend","soma","Adend1","Adend2","Adend3"]
        #pos_vals = ["soma"]
    elif(ts_choise=="synapses_pyr"):
        all_comp=["Bdend","soma","Adend1","Adend2","Adend3"]
        df = pd.concat([*[df.filter(like=comp).filter(like="mean").sum(axis=1) for comp in all_comp],*[df["iseed"]]],axis=1,keys=[*all_comp,"iseed"])
        df.fillna(0,inplace=True)
        pos_vals = all_comp
    elif("volt" in ts_choise):
        pos_vals=[k for k in list(df.keys()) if "volt_mean" in k]
        all_comp=pos_vals
        df =pd.concat([*[df.filter(like=comp).filter(like="volt_mean").sum(axis=1) for comp in all_comp],*[df["iseed"]]],axis=1,keys=[*all_comp,"iseed"])
    elif("external_inputs_pyr" in ts_choise):
        w=[1,5,150][input3]
        mode=["exp","gaussian"][input4]
        ratesnam = "rates_pyr_w"+str(w)+"_mode"+mode
        df = continufySpikes2(df,fs,ratesnam,w=w,mode=mode,dt=1000/fs,tmax=np.max(df["tvec"].values))
        pos_vals=["none"]

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
            ts_choise = ratesnam
            df_pos = df.copy()
        ftheta,fgamma,psi_tmp,nu_tmp= get_avg_cfd(fgamma_lims,ftheta_lims,fs,df_pos,ts_choise,surrogate_test=surrogate_test,nperseg=nperseg,testing=testing,flag_z=flag_z)
        all_psi[pos]=np.mean(psi_tmp,axis=0)
        all_psi_std[pos]=np.std(psi_tmp,axis=0)
        all_nu[pos]=np.mean(nu_tmp,axis=0)
        all_nu_std[pos]=np.std(nu_tmp,axis=0)
        

    if(surrogate_test):
        title = "CFD4_fft"+str(nperseg)+"_z"+str(flag_z)+"_sign"+testing+"_"+ts_choise+"_"+area
    else:
        title = "CFD_fft"+str(nperseg)+"_z"+str(flag_z)+"_"+ts_choise+"_"+area

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

#os.remove("CFD_out"+str(input1)+".dat")
#os.remove("CFD_err"+str(input1)+".log")
os.remove("CFD_out"+".dat")
os.remove("CFD_err"+".log")