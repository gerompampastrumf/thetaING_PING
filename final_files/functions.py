import pandas as pd 
import numpy as np 

def process_spike_data(data_unified):
    ''' process spike data. Produces a dataframe with columns: tvec, idvec'''
    data = {"tvec":[],"idvec":[]}
    aux = data_unified[0]
    for i in aux: 
        if aux[i]:
            data["tvec"].extend(aux[i])
            data["idvec"].extend(i*np.ones(len(aux[i])))

    return pd.DataFrame(data,columns = ["tvec","idvec"])

def process_volt_data(data_unified):
    ''' process voltage data. Produces a dataframe with columns: key_volt_mean, key_vot_std'''    
    keys = data_unified.keys()
    data = {}
    for key in keys:
        aux = []
        for i in data_unified[key][0]:
            aux.append(data_unified[key][0][i])
        data[f"{key}_volt_mean"] = np.mean(aux,axis=0)
        data[f"{key}_volt_std"]  = np.std(aux,axis=0)
    return pd.DataFrame(data,columns = data.keys())    

def process_syn_data( data_unified ): #syncurrents
    ''' process synaptic current data. Produces a dataframe with columns: ikey_mean, ikey_std'''
    keys = data_unified.keys()
    data = {}
    for key in keys:
        aux = []
        for i in data_unified[key][0]:
            aux.append(data_unified[key][0][i])
        data[f"i{key}_mean"] = np.mean(aux,axis=0)
        data[f"i{key}_std"]  = np.std(aux,axis=0)
    return pd.DataFrame(data,columns = data.keys())
