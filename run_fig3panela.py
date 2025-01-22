#RUNS the theta-ING motif given a seed (input1=[1,20])
#and a multiplier for the bas to pyr strength.

import warnings
warnings.filterwarnings("ignore")
from neuron import h
import numpy as np
import sys
import os
from neuron.units import ms, mV
import time as tm
import shutil
sys.path.append('/home/dimitrios/Neurons/thetaING_PING/final_files/')

from network_hippocampus_ca3 import *
from LFPsimpy import LfpElectrode
import file_management
from functions import *

tick = tm.time()

# inputs values from bash
inputs_argvs = [] 
column_labels = []
for argv in sys.argv[1:]:
    inputs_argvs.append(int(argv))

if len(inputs_argvs) > 2:
    input1 = inputs_argvs[0] # background noise basal seed
    input2 = inputs_argvs[1] 
    input3 = inputs_argvs[2] 
    #input2 = inputs_argvs[3] # second variable to scan
    # must be added the other variables in case they are used
    # ... = inputs[2]
    # ... = inputs[3]
    # ...
    number_of_argvs = len(sys.argv[1:])
    argvs=""
    for i in range(number_of_argvs-1):
        argvs += sys.argv[1:][i]+'_'
    argvs+=sys.argv[-1]

    column_labels.append("iseed")
    column_labels.append("input2")
    # column_labels ...

elif len(inputs_argvs) == 2: 

    input1 = inputs_argvs[0] # background noise basal seed
    input2 = inputs_argvs[1]
    # control the seed for the background noise generation (online generation)
        # and would also select the trial of the background noise (offline generation, not implemented)
    number_of_argvs = len(sys.argv[1:])
    argvs=""
    for i in range(number_of_argvs-1):
        argvs += sys.argv[1:][i]+'_'
    argvs+=sys.argv[-1]

    column_labels.append("iseed")
    column_labels.append("input2")
else: 
    print("provide at least two inputs")
    sys.exit()

simulation_time = 60000.0 # total time of the simulation
time_resolution = 1.0 # time resolution of the measurements

DoMakeNoise = True # Control the external noise
MakeNetStim = True # Online generation of the background noise 
DoMakeExternalInputs = True # Control the external inputs
MakeCellStim         = False # Online generation of the external inputs

create_ca1_network  = False # Create the CA1 network as well
record_all_synapses = True # Record all the synapses
record_lfp          = True # Record the LFP and transmembrane currents
record_spikes_dendr = {"Bdend":None,"Adend1":None,"Adend2":None,"Adend3":-30}
synlist = [ "Adend3AMPA_ec2180","Adend3NMDA_ec2180","Adend1AMPA_dgreg", "Adend1NMDA_dgreg","somaGABA_bas", "somaAMPA_noise", "somaGABA_noise"]

record_mp = {"Bdend": False, "soma": True, "Adend1": False, "Adend2": False, "Adend3": False}

# what data to save (be consistent to what you record)
save_data_volt   = True
save_data_spikes = True 
# these are the ones that require consistency what previously declared
save_data_syn    = True
save_data_lfp    = True
save_data_ica    = True

current_folder = os.getcwd()
inputs_folder = os.path.join(current_folder, "sigma0_2b/external_inputs4")#0_2b stands for width=0.2*125ms, 1e4 neurons

file_name = __file__
file_name = file_name.replace(current_folder,"") 
file_name = file_name.split('_')[-1][:-3] 
save_folder = os.path.join(current_folder, file_name) +"-"+str(input2)
if not os.path.exists(save_folder+"/external_inputs"):
    os.makedirs(save_folder+"/external_inputs")
#shutil.copy(inputs_folder+"/external_inputs_theta_gen.lzma",save_folder+"/external_inputs/external_inputs_theta_gen.lzma")
iseed=input1

def unify_data(variable, cells, nr=1):
    f''' unify data from all cells in the network: membrane voltage, synaptic currents, spikes, etc '''
    local_data = { cell.id: list( np.round( cell.__dict__[variable].to_python(),nr)) for cell in cells }
    all_data = pc.py_alltoall( [local_data] + [None] * (pc.nhost() - 1) )
    return all_data

"""
#==============PARAMETERS====================
"""

''' ############################################################################
            Parameters: neurons -> neurons
#############################################################################'''

i_bas = [-1e-2,-1e-3,0.0,1e-3,1e-2][2]
#i_pyr = 0.05

wTheta_pyr = np.linspace(0,2,11)[6]
nTheta_pyr = np.linspace(0,2,11)[3]
bas_bas = [0.0,0.5,1.2,2.0,3.0][2]
bas_pyr = [0.0,0.1,0.2,0.3][input2]
pyr_bas = [0.5,1.0,2.0,5.0,10.0][1]*0.0

n_ec_neurons=500
noise_bas = 4.0
noise_bas_neg = [1.2,1.5,2.0,5.0,10.0][-1]
ec_bas = [0.2,0.3,0.4,0.5,0.6][3]
ec_pyr=[0.0,0.03,0.04,0.05,0.06][2]
pyr_gaba_noise=[2.0,5.0,10.0,50.0,100.0][2]
pyr_noise = [1.75,2.1,2.2,2.3,2.4][2]
i_pyr = 0.0
isi_noise = 1.0
bas_pyr2=[1.0,2.0,5.0,10.0,20.0][2]
theta_delay = [0,10,20,30,40][2]

weights_neurons_neurons = {}
nsyns_neurons_neurons = {}
syn_neurons_neurons  = {}
delay_neurons_neurons = {}

weights_neurons_neurons["bas_ca3_to_pyr_ca3"] = [[0.576e-3*bas_pyr*bas_pyr2/2]] # [0.8*0.72e-3]
weights_neurons_neurons["pyr_ca3_to_bas_ca3"] = [[1.0e-3*pyr_bas/2,1.0e-3*0.1*pyr_bas/2]] #[0.36e-3+scale*1.38e-3, 1.38e-3]
weights_neurons_neurons["bas_ca3_to_bas_ca3"] = [[4.5e-3*bas_bas*1.5*2.0/2]]

nsyns_neurons_neurons["bas_ca3_to_pyr_ca3"] = [15*2]
nsyns_neurons_neurons["pyr_ca3_to_bas_ca3"] = [40*2]
nsyns_neurons_neurons["bas_ca3_to_bas_ca3"] = [15*2]

syn_neurons_neurons["bas_ca3_to_pyr_ca3"] = [["somaGABA_bas"]]                   # [["somaGABA"]]
syn_neurons_neurons["pyr_ca3_to_bas_ca3"] = [["somaAMPA_pyr","somaNMDA_pyr"]]    # [["somaAMPAf", "somaNMDA"]]
syn_neurons_neurons["bas_ca3_to_bas_ca3"] = [["somaGABA_bas"]]                   # [["somaGABAf"]]

delay_neurons_neurons["bas_ca3_to_pyr_ca3"] = [1.5,0.2]
delay_neurons_neurons["pyr_ca3_to_bas_ca3"] = [1.5,0.2]
delay_neurons_neurons["bas_ca3_to_bas_ca3"] = [1.5,0.2]


''' ############################################################################
            Parameters: inputs -> neurons
#############################################################################'''
weights_inputs_neurons = {}
nsyns_inputs_neurons   = {}
syn_inputs_neurons     = {}
delay_inputs_neurons   = {}

weights_inputs_neurons["theta_gen_to_pyr_ca3"]= [[1.0e-3*ec_pyr,1.0e-3*0.1*ec_pyr]]
weights_inputs_neurons["theta_gen_to_bas_ca3"]= [[1.0e-3*ec_bas,1.0e-3*0.1*ec_bas]]
delay_inputs_neurons["theta_gen_to_pyr_ca3"] = [theta_delay,0.2]
delay_inputs_neurons["theta_gen_to_bas_ca3"] = [20.0+1.0,0.2]
syn_inputs_neurons["theta_gen_to_pyr_ca3"]    = [["Adend3AMPA_ec2180","Adend3NMDA_ec2180"]]
syn_inputs_neurons["theta_gen_to_bas_ca3"]    = [["somaAMPA_ec2180","somaNMDA_ec2180"]]
nsyns_inputs_neurons["theta_gen_to_pyr_ca3"]  = [n_ec_neurons]
nsyns_inputs_neurons["theta_gen_to_bas_ca3"]  = [n_ec_neurons]

keys = ["pyr_ca3","bas_ca3"]
weights_noise_neurons = dict.fromkeys(keys)
for key in keys:
    weights_noise_neurons[key] = {}

weights_noise_neurons["pyr_ca3"]["Adend1NMDA_dgreg"]   = 1e-05*isi_noise/10.0
weights_noise_neurons["pyr_ca3"]["Adend1AMPA_dgreg"]   = 1e-05*isi_noise
weights_noise_neurons["pyr_ca3"]["somaAMPA_noise"]   = 0.05e-3*isi_noise*2.0*pyr_noise
weights_noise_neurons["pyr_ca3"]["somaGABA_noise"]   = 0.012e-3*pyr_gaba_noise*isi_noise*2.0
weights_noise_neurons["olm_ca3"]={}
weights_noise_neurons["bas_ca3"]["somaAMPA_noise"]   = 0.02e-3*noise_bas*isi_noise
weights_noise_neurons["bas_ca3"]["somaGABA_noise"]   = 0.2e-3*noise_bas_neg*isi_noise
"""
isi_pyrNoise=1.05 in original here I put the default value 1.0
isi_noise_neurons["pyr_ca3"]["Adend1NMDA_pyrCA3"] = isi_pyrNoise
isi_noise_neurons["pyr_ca3"]["Adend1AMPA_pyrCA3"] = isi_pyrNoise
isi_noise_neurons["pyr_ca3"]["somaAMPA_noise"] = 1.0
isi_noise_neurons["pyr_ca3"]["somaGABA_noise"] = 1.0
isi_noise_neurons["bas_ca3"]["somaAMPA_noise"] = 1.0
isi_noise_neurons["bas_ca3"]["somaGABA_noise"] = 1.0
"""
"""
#==============CREATION====================
"""
net = Network( weights_inputs_neurons  = weights_inputs_neurons,
               nsyns_inputs_neurons    = nsyns_inputs_neurons,
               syn_inputs_neurons      = syn_inputs_neurons,
               weights_neurons_neurons = weights_neurons_neurons,
               nsyns_neurons_neurons   = nsyns_neurons_neurons,
               syn_neurons_neurons     = syn_neurons_neurons,
               weights_noise_neurons   = weights_noise_neurons,
               delay_neurons_neurons   = delay_neurons_neurons,
               delay_inputs_neurons    = delay_inputs_neurons,
               iseed=3421*(2*iseed+1),  # seed for external inputs, this parameter will not be used (only for "online" generation)
               bseed=100003*(2*iseed+1),   # seed for background inputs    
               create_ca1_network = create_ca1_network,
               DoMakeNoise = DoMakeNoise,
               DoMakeExternalInputs = DoMakeExternalInputs,
               MakeCellStim = MakeCellStim,
               MakeNetStim  = MakeNetStim,
               inputs_folder = inputs_folder,
               external_inputs_iseed = iseed, # !!!
               n_pyr_ca3 = 200,
               n_bas_ca3 = 40,
               n_olm_ca3 = 1,
               n_pyr_ca1 = 0,
               n_bas_ca1 = 0,
               n_olm_ca1 = 0,
               n_cck_ca1 = 0,
               nseg_pyr  = 3,
               record_mp = record_mp,
               resolution = time_resolution,
               i_bas=i_bas,
               i_pyr=i_pyr*40e-3,
               ISI_noise=isi_noise,
               record_spikes_dendr=record_spikes_dendr)

net.set_noise_inputs(simulation_time) # set background noise and external input
if net.MakeNetStim:
    net.init_NetStims()  # init rngs of background

'''###########################################################################
                     external inputs vecstims
############################################################################'''

inputs, vecstims = [],[]
net.inputs_list = []
net.vsl_ = []

if net.DoMakeExternalInputs:
    print("external inputs with vecstims")
    print(net.external_inputs_data.keys())
    for key in net.external_inputs_data.keys():
        if key.endswith("ca3"):
            print(key)
            idv      = net.external_inputs_data[key]["idv"]
            spikes   = net.external_inputs_data[key]["spikes"]
            trg      = net.external_inputs_data[key]["population"]
            syn_list = net.external_inputs_data[key]["synapses"]
            delays   = net.external_inputs_data[key]["delays"]
            w_list   = net.external_inputs_data[key]["weights"]
            conn     = net.external_inputs_data[key]["connectivity"]

            for k,conn_ in enumerate(conn):
                for post_id, all_pre in enumerate(conn_):
                    net.inputs_list.append([ ])
                    for j, pre_id in enumerate(all_pre):
                        idvs = np.where(idv==pre_id)[0]
                        net.inputs_list[-1].append( np.sort(spikes[idvs]))
                        inputs.append(h.Vector( np.sort(spikes[idvs]) ))
                        vecstims.append(h.VecStim())
                        vecstims[-1].play(inputs[-1])

                        spike_list = vecstims[-1]
                        net.vsl_.append(spike_list)
                        for syn,w in zip(syn_list[k], w_list[k]):
                            net.ncl_.append(h.NetCon(spike_list, trg.cell[post_id].__dict__[syn].syn, 0, delays[k][j,post_id], np.random.normal(loc=w,scale=w*0.2)))

    tock = tm.time()
    diff = tock-tick
    horas = int(diff/3600)
    diff  = diff-horas*3600
    mint  = int(diff/60)
    diff  = diff-mint*60
    seg   = int(diff)

    print('Vecstim done after:', horas,'h', mint, 'min', seg, 's')
''' #######################################################################'''

pc = h.ParallelContext()
pc.set_maxstep(100*ms)
t = h.Vector().record(h._ref_t)

# h.steps_per_ms = 1/h.dt


if record_all_synapses:
    # this record only the specific synapses in the ca3 pyramidal cells
    # there is a function called record_all_synapses() in the net class that record all synapses in the network
    # but the synlist should be provided first to the neurons.
    for cell in net.pyr_ca3.cell:
        cell.syn_list = synlist
        cell.record_synapses()

#Record LFP
if record_lfp:
    # electrode_y_coords = np.arange(0,850,50)
    zcoords = [10,385]
    electrode = {}
    electrode["ca3"] = []

    for i, z in enumerate(zcoords): 
        electrode["ca3"].append( LfpElectrode(x=25.0, y=25.0, z=z, sampling_period=time_resolution, neuron_type = "Pyramidal CA3"))
      
h.init()
h.finitialize()
h.fcurrent()
h.frecord_init()
h.tstop=simulation_time
h.dt = 0.1
h.celsius = 34
pc.psolve(simulation_time*ms) # simulation running starts

tock = tm.time()
diff = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

tock = tm.time()
diff = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

print("Simulation done after", horas,'h', mint, 'min', seg, 's')
print("-------------------------------------------")
print("Preparing the data to be saved") 

t = np.array(t.to_python())
all_spikes, all_volt, lfp, ica = {},{},{},{}

# only saving the data if the simulation was run in the master node

if pc.id() == 0:

    if save_data_volt:
        #############################################
        print("Unifying the membrane potential data")
        #############################################
        data_volt = {}
        all_volt  = dict.fromkeys(["pyr_ca3","bas_ca3","olm_ca3"])
        keys = [] 
        for key in record_mp.keys():
            if record_mp[key] == True:
                keys.append(key)
        if keys: 
            all_volt["pyr_ca3"] = dict.fromkeys(keys)
            for key in keys:        
                all_volt["pyr_ca3"][key] = unify_data(f"{key}_volt", cells=net.__dict__["pyr_ca3"].cell,nr=4)
            
            data_volt["pyr_ca3"] = process_volt_data(all_volt["pyr_ca3"])

            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                data_volt["pyr_ca3"][column_labels[i]] = [value]*len(data_volt["pyr_ca3"])
            
            title = f"volt_pyr_ca3_{argvs}.lzma"
            file_management.save_lzma( data_volt["pyr_ca3"], title, parent_dir = save_folder)

            if record_mp["soma"] == True:
                all_volt["bas_ca3"] = dict.fromkeys(["soma"])
                all_volt["olm_ca3"] = dict.fromkeys(["soma"])
                all_volt["bas_ca3"]["soma"] = unify_data("soma_volt", cells=net.__dict__["bas_ca3"].cell,nr=4)
                all_volt["olm_ca3"]["soma"] = unify_data("soma_volt", cells=net.__dict__["olm_ca3"].cell,nr=4)

                data_volt["bas_ca3"] = process_volt_data(all_volt["bas_ca3"])
                data_volt["olm_ca3"] = process_volt_data(all_volt["olm_ca3"])

                for i in range(number_of_argvs): # adding the input values
                    value = inputs_argvs[i]
                    data_volt["bas_ca3"][column_labels[i]] = [value]*len(data_volt["bas_ca3"])
                    data_volt["olm_ca3"][column_labels[i]] = [value]*len(data_volt["olm_ca3"])

                title = f"volt_bas_ca3_{argvs}.lzma"
                file_management.save_lzma( data_volt["bas_ca3"], title, parent_dir = save_folder)
                title = f"volt_olm_ca3_{argvs}.lzma"
                file_management.save_lzma( data_volt["olm_ca3"], title, parent_dir = save_folder)

    if save_data_spikes:
        ################################
        print("Unifying the spike data")
        ################################
        all_spikes = dict.fromkeys(["pyr_ca3","bas_ca3","olm_ca3"])
        data_spikes = {} 

        all_spikes["pyr_ca3"] = unify_data("spike_times", net.__dict__["pyr_ca3"].cell)
        all_spikes["bas_ca3"] = unify_data("spike_times", net.__dict__["bas_ca3"].cell)
        all_spikes["olm_ca3"] = unify_data("spike_times", net.__dict__["olm_ca3"].cell)

        data_spikes["pyr_ca3"] = process_spike_data( all_spikes["pyr_ca3"] )
        data_spikes["bas_ca3"] = process_spike_data( all_spikes["bas_ca3"] )
        data_spikes["olm_ca3"] = process_spike_data( all_spikes["olm_ca3"] )

        for i in range(number_of_argvs): # adding the input values
            value = inputs_argvs[i]
            data_spikes["pyr_ca3"][column_labels[i]] = [value]*len(data_spikes["pyr_ca3"])
            data_spikes["bas_ca3"][column_labels[i]] = [value]*len(data_spikes["bas_ca3"])
            data_spikes["olm_ca3"][column_labels[i]] = [value]*len(data_spikes["olm_ca3"])
        
        title = f"external_inputs_pyr_ca3_{argvs}.lzma"
        file_management.save_lzma( data_spikes["pyr_ca3"], title, parent_dir = save_folder)
        title = f"spikes_bas_ca3_{argvs}.lzma"
        file_management.save_lzma( data_spikes["bas_ca3"], title, parent_dir = save_folder)
        title = f"spikes_olm_ca3_{argvs}.lzma"
        file_management.save_lzma( data_spikes["olm_ca3"], title, parent_dir = save_folder)

    dendr_spikes={"tvec":[],"idvec":[],"component":[]}
    recorded_dend_sp=False
    data_dendr_spikes={}
    for key in record_spikes_dendr:  
        if(record_spikes_dendr[key] is not None):
            print("DENDR recorded "+key)
            recorded_dend_sp = True
            dendr_spikes_temp = unify_data(key+"_spike_times", net.__dict__["pyr_ca3"].cell)
            dendr_spikes_temp = process_spike_data(dendr_spikes_temp)
            dendr_spikes["tvec"].extend(np.array(dendr_spikes_temp["tvec"]))
            dendr_spikes["idvec"].extend(np.array(dendr_spikes_temp["idvec"]))
            dendr_spikes["component"].extend(np.array([key]*len(dendr_spikes_temp["idvec"])))
            dendr_spikes=pd.DataFrame(dendr_spikes)
            for i in range(number_of_argvs): # adding the input values
                value = inputs_argvs[i]
                dendr_spikes[column_labels[i]] = [value]*len(dendr_spikes_temp)
    if(recorded_dend_sp):
        title = f"dendr_spikes_pyr_ca3_{argvs}.lzma"
        file_management.save_lzma( dendr_spikes, title, parent_dir = save_folder)

    if save_data_syn:
        if record_all_synapses:
            ################################
            print("Unifying the synaptic currents")
            ################################

            all_syn = {}
            all_syn["ca3"] = {}
            for key in synlist:
                all_syn["ca3"][key] = unify_data(f"i{key}", cells=net.__dict__["pyr_ca3"].cell,nr=4)
            data_syn = {} 
            data_syn["ca3"] = {}
            data_syn["ca3"] = process_syn_data( all_syn["ca3"] )

            for i in range(number_of_argvs): # adding the inpOut values
                value = inputs_argvs[i]
                data_syn["ca3"][column_labels[i]] = [value]*len(data_syn["ca3"])
            
            title = f"synapses_pyr_ca3_{argvs}.lzma"
            file_management.save_lzma( data_syn["ca3"], title, parent_dir = save_folder)

    if save_data_lfp:

        # LFPsimpy is already implemented with the parallelization framework
        if record_lfp:
            data_lfp = {}
            data_lfp["ca3"] =  {"lfp":[],"electrode":[]}
            
            for i, lfp_ in enumerate(electrode["ca3"]):
                data_lfp["ca3"]["lfp"].extend( np.array(lfp_.values) )
                data_lfp["ca3"]["electrode"].extend( [zcoords[i]]*len(lfp_.values) ) # z-component inseatd of names, in case we put more
            
            data_lfp["ca3"] = pd.DataFrame(data_lfp["ca3"])
            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                # data_lfp["ca3"][f"input{i+1}"] = []
                data_lfp["ca3"][column_labels[i]] = [value]*len(data_lfp["ca3"]) 
            
            title = f"lfp_ca3_{argvs}.lzma"
            file_management.save_lzma(data_lfp["ca3"], title, parent_dir=save_folder)
    
    if save_data_ica:
        if record_lfp:
            # LFPsimpy is already implemented with the parallelization framework
            data_ica = {}
            data_ica["ca3"] = {"ica":[],"electrode":[],"component":[]}
            
            for i, lfp_ in enumerate(electrode["ca3"]):
                for j,sublfp_ in enumerate(np.array(lfp_.values_per_section).T):
                    data_ica["ca3"]["ica"].extend( np.array(sublfp_) )
                    data_ica["ca3"]["electrode"].extend( [zcoords[i]]*len(sublfp_) )
                    data_ica["ca3"]["component"].extend( [["Bdend","soma","Adend1","Adend2","Adend3"][j]]*len(sublfp_) )

            data_ica["ca3"] = pd.DataFrame(data_ica["ca3"])
            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                data_ica["ca3"][column_labels[i]] = [value]*len(data_ica["ca3"]) 
            
            title = f"ica_ca3_{argvs}.lzma"
            file_management.save_lzma(data_ica["ca3"], title, parent_dir=save_folder)



error_file = "error_"+str(input1)+"_"+str(input2)+".log"
output_file = "output_"+str(input1)+"_"+str(input2)+".dat"
os.remove(error_file)
os.remove(output_file)

pc.barrier()
pc.done()
h.quit()