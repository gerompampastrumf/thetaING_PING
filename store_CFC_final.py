import sys
import os
import pandas as pd
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/simplifiedCFD/final_files')
import file_management
import time as tm
from signal_analysis import compute_mi
from signal_analysis import compute_coherence_measurements
from signal_analysis import compute_mi_conelty
from signal_analysis import continufySpikes2
tick = tm.time()
#Calculate the avg modulation inxed over iseed and bseed
def get_avg_mi(fgamma,ftheta,fs,df,label,surrogate_test=False,flag_cyclebytheta=False,flag_canolty=False,flag_z=False): 
    all_mi = []
    all_inputs1=np.unique(df["iseed"].values)
    n_temps=0
    for input1 in all_inputs1:
        n_temps+=1
        ts = df[label][(df['iseed'] == input1)].values
        ts=ts[len(ts)//15:]
        if(flag_z):
            ts-=np.mean(ts)
            ts/=np.std(ts)
        if(flag_canolty):
            ftheta, fgamma, mi_temp, mi_sur_temp = compute_mi_conelty(ts, ftheta, fgamma, fs=fs,
                                                    surrogate_test=surrogate_test, nsurrogates=100)
        else:
            ftheta, fgamma, mi_temp, mi_sur_temp = compute_mi(ts, ftheta, fgamma, fs=fs, 
                                      nbins=20, surrogate_test=surrogate_test, nsurrogates=100
                                      ,choiseSurr = 2,flag_cyclebytheta=flag_cyclebytheta)
        
        mi_temp = np.swapaxes(mi_temp,0,1)
        mi_sur_temp= np.swapaxes(mi_sur_temp,1,2)
        if(surrogate_test):
            inds_insig = mi_temp<=np.quantile(np.array(mi_sur_temp),0.995,axis=0)
            mi_temp[inds_insig] = 0
            
        all_mi.append(mi_temp)
    print(n_temps)
    return all_mi

input1 = int(sys.argv[1])
input2 = sys.argv[2]
input3 = int(sys.argv[3])
input4 = int(sys.argv[4])
current_folder = os.getcwd()
folder = os.path.join(current_folder, input2)
flag_z = True
ts_choise = ["external_inputs/external_inputs_pyr","ica","volt_pyr","volt_bas","lfp","synapses_pyr"][input1]

fs = 500*2
area="ca3"
surrogate_test=False
flag_canolty=True
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
    pos_vals = ["Bdend","soma","Adend1","Adend2","Adend3"]
elif(ts_choise=="synapses_pyr"):
    all_comp=["Bdend","soma","Adend1","Adend2","Adend3"]
    df = pd.concat([*[df.filter(like=comp).filter(like="mean").sum(axis=1) for comp in all_comp],*[df["iseed"]]],axis=1,keys=[*all_comp,"iseed"])
    ts_ = df[all_comp[0]].values
    df.fillna(0,inplace=True)
    pos_vals = all_comp
elif("volt" in ts_choise):
    #pos_vals=["soma_volt_mean"]
    pos_vals=[k for k in list(df.keys()) if "volt_mean" in k]
    all_comp=pos_vals
    df =pd.concat([*[df.filter(like=comp).filter(like="volt_mean").sum(axis=1) for comp in all_comp],*[df["iseed"]]],axis=1,keys=[*all_comp,"iseed"])
    ts_ = df[pos_vals[0]][(df["iseed"]==ised)].values
    print("ts_ len: ",len(ts_))
elif("external_inputs_pyr" in ts_choise):
    w=[1,5,150][input3]
    mode=["exp","gaussian"][input4]
    ratesnam = "rates_pyr_w"+str(w)+"_mode"+mode
    df = continufySpikes2(df,fs,ratesnam,w=w,mode=mode,dt=1000/fs,tmax=np.max(df["tvec"].values))
    pos_vals=["none"]
    ts_ = df[ratesnam][(df["iseed"]==ised)].values
    print("ts_ len: ",len(ts_))

ftheta_min,ftheta_max = 2,15
fgamma_min,fgamma_max,dfgamma = 20,120,1
#nperseg = int(fs/500*1024)
nperseg = int(fs/500*1024)
all_mi = {}
all_mi_std = {}
flag_cyclebytheta=False

ftheta, fgamma, _, _, _, _ = compute_coherence_measurements(ts_, fs, nperseg, 0.5*nperseg, ftheta_min=ftheta_min,
                            ftheta_max=ftheta_max, fgamma_min=fgamma_min, fgamma_max=fgamma_max, 
                            dfgamma=dfgamma,norder=4, nfft=nperseg,surrogate_test=False, nsurrogates=1,choiseSurr = None)
for ind,pos in enumerate(pos_vals):
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
    all_mi_temp =get_avg_mi(fgamma,ftheta,fs,df_pos,ts_choise,surrogate_test=surrogate_test,flag_cyclebytheta=flag_cyclebytheta,flag_canolty=flag_canolty,flag_z=flag_z) 
    all_mi[pos]= np.mean(all_mi_temp,axis=0)
    all_mi_std[pos]= np.std(all_mi_temp,axis=0)

if(surrogate_test):
    title = "CFC_Can"+str(flag_canolty)+"_fft"+str(nperseg)+"_z"+str(flag_z)+"_sign_test"+ts_choise+"_"+area
else:
    title = "CFC_Can"+str(flag_canolty)+"_fft"+str(nperseg)+"_z"+str(flag_z)+"_"+ts_choise+"_"+area

file_management.save_lzma([all_mi,all_mi_std,fgamma,ftheta],title,parent_dir=folder)

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

os.remove("CFC_out"+str(input1)+".dat")
os.remove("CFC_err"+str(input1)+".log")
