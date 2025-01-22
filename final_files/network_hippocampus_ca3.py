"""

Ultima modificacion 03/11/2022 
nseg: nuevo argumento cuando se definen las neuronas  

Ultima modifiación 08/09/2022 
Se modifica las carpetas donde están almacenados los external inputs. 
Ahora hay dos carpetas: external_inputs_neuron y external_inputs_python. 
Dentro de las cuales estan los inputs de todas aquellas poblaciones excepto la del dg burst.
Para esta hay un total de 21 subcarpetas indicando el nivel de burst.  

Mejor, vamos a hacer que el file_folder, lo vamos a llamar inputs_folder y paso 
directamente como entrada la carpeta donde estan los inputs.

Ultima modificacion 07/06/2022

Esta version empieza como una copia de la version optimiazada de hippocampo ca3 netork_test3.py
(07/09/2021)

Me dispongo aqui a hacer la red completa del hippocampo con todas sus conexiones

Update 19/02/2022
De momento la neurona CCK está omitida.

Falta por implementar el factor de escala. Este factor sirve para determinar los
autenticos valores de las conductancias de AMPA puesto que en verdad el NMDA
que utiliza este modelo es una combinacion de AMPA y NMDA, cosa que no he llegado
a entender. Por tanto lo que hice en su momento (y de recuperar/recordar) para
implementarlo aqui

Update 11/03/2022
(spanish) el codigo ha sido modificado de nuevo de tal manera que la dinamica de ca3
se reproducira de la misma manera independiente de si se simula ca1. Por otro lado,
he creado una nueva seed para el wiring de los inputs externos(EC,DG,..) para satisfacer
lo anteriormente comentado. Ademas se crea la opcion para guardar los inputs
generados como background, que alimentan a cada neurona.

Update 23/03/2022
(spanish) añado otra semilla mas para la genereacion de una distribucion de delays
para todas las conexiones

Last update 24/03/2022

Last update 28/03/2022 inputs de giro dentado a las basket, no estaba especificado.

Last update 31/03/2022 Asilo los gains de ec2 y dg en sus correspondientes contribuciones...
de cara a poder aislar cada uno analizar la correlacion entre los inputs.

#Unifico los inputs de AMPA/NMDA en las mismas conexiones. Por tanto, necesito obtener primero las conexiones

# Implemento el VecStim a ver si me va...

# 10/05/2022
# Cleaner version

"""

from neuron import h
import numpy as np
import random
import sys
import os
# import warnings
# sys.path.append('/home/jaime/Desktop/hippocampus/files')
from neurons import *
import file_management
# from external_inputs import *

#warnings.simplefilter(action='ignore', category=FutureWarning)
#warnings.simplefilter(action='ignore', category=RuntimeWarning
pc = h.ParallelContext()

class RandomStream():
    # this is also obselete, but I left it here just in case
    def __init__(self, sead, generator):
        self.generator = generator
        self.r = h.Random()
        self.sead = sead
        self.start()

    def set_generator(self):
        if self.generator == 'gaussian':
            g = self.r.normal(0,1)
        if self.generator == 'exponential':
            g = self.r.negexp(1)
        return g

    def start(self):
        return self.r.MCellRan4(self.sead,self.sead)

    def repick(self):
        return self.r.repick()

class Population:
    "Population of cells"
    # cell_type -- pyr, bas, olm
    # n -- number of cells in the population
    # x, y, z -- initial position for the first Cell
    # dx -- an increment of the x-position for the cell location
    # amp, dur, delay -- parameters for the IClamp in the soma
    def __init__(self, cell_type, pgidlist, n , x, y, z, dx, amp, dur, delay, nseg, record_mp, resolution,record_spikes_dendr, pyr=False, naux=0,
    x_position = None, y_position = None):
        self.cell = [] # List of cells in the population
        self.nc   = [] # NetCon list for recording spikes
        self.pgidlist = pgidlist
        self.ng   = len(self.pgidlist) # number of cells in host
        self.n    = n  # number of cells
        self.x    = x
        self.y    = y
        self.z    = z
        self.record_spikes_dendr = record_spikes_dendr
        #self.ltimevec = h.List() # list of Vectors for recording spikes, one per cell
        #self.lidvec = h.List()
        #self.tvec = h.Vector()
        #self.idvec = h.Vector()
        self.nssidx = {}
        self.nseidx = {}
        self.ncsidx = {}
        self.nceidx = {}
        self.ncd  = [] # only for pyramidal

        self.naux = int(naux) # to distinguish gid from i-th position of the neruon in the network
        self.nc_conn = []
        if pyr:
            for i in self.pgidlist: #range(n):
                #self.cell.append(cell_type(x+i*dx,y,z,gGID))
                self.cell.append(cell_type( x_position[i],y_position[i],z,i, nseg, record_mp, resolution,record_spikes_dendr=self.record_spikes_dendr))
                self.cell[-1].somaInj.amp   = amp
                self.cell[-1].somaInj.dur   = dur
                self.cell[-1].somaInj.delay = delay
                pc.cell(i, self.cell[-1].spike_detector)
        else:
            for i in self.pgidlist: #range(n):
                #self.cell.append(cell_type(x+i*dx,y,z,gGID))
                self.cell.append(cell_type(x+(i-self.naux)*dx,y,z,i, nseg, record_mp, resolution))
                self.cell[-1].somaInj.amp   = amp
                self.cell[-1].somaInj.dur   = dur
                self.cell[-1].somaInj.delay = delay
                pc.cell(i, self.cell[-1].spike_detector)
        for section in self.cell[0].sec_list: # define a vector per each comparment/secion to later save the synaptic currents
            self.__dict__[section+"_isyn"] = {}

        for syn in self.cell[0].syn_list:
            self.__dict__[syn+"_mean"] = []
            self.__dict__[syn+"_std"]  = []

    def set_r(self, syn, r):
        for c in self.cell:
            c.__dict__[syn].syn.r = r

''' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'''
class Network:

    def __init__(self, weights_inputs_neurons, nsyns_inputs_neurons, syn_inputs_neurons,
               weights_neurons_neurons, nsyns_neurons_neurons, syn_neurons_neurons,
               weights_noise_neurons, nseg_pyr, delay_inputs_neurons, delay_neurons_neurons, 
               connections=True,
               DoMakeNoise=True,
               DoMakeExternalInputs=True,
               MakeNetStim=True,
               MakeCellStim=True,
               bseed=1234,
               wseed=4321,
               iseed=3421,
               jseed=2134,
               SaveConn=False,
               noise_burst=0.5,
               create_ca1_network=False,
               record_noise=False,
               ISI_noise = 1.0,
               jitter=0.2,
               delays_in_noise=True,
               delay_mean_noise=2.0,
               delay_mean_inputs=2.0,
               delay_mean_connections=2.0,
               save_input_connection=False,
               topology="rectangular",
               inputs_folder = '/home/jaime/Desktop/hippocampus/external_inputs/baseline/',
            #    burst_level_label = '10',
               external_inputs_iseed=0,
               n_pyr_ca3 = 800, 
               n_bas_ca3 = 100, 
               n_olm_ca3 = 30, 
               n_pyr_ca1 = 800, 
               n_bas_ca1 = 100, 
               n_olm_ca1 = 30, 
               n_cck_ca1 = 70, 
               record_mp = {"soma":True, "Bdend":False, "Adend1":False, "Adend2":False, "Adend3":False},
               record_spikes_dendr = {"Adend1":None,"Adend2":None,"Adend3":None,"Bdend":None},
               resolution = 1,
               i_bas=0,
               i_pyr=40e-3):

        self.nseg_pyr = nseg_pyr
        print("Setting Cells")
        self.record_spikes_dendr = record_spikes_dendr
        self.weights_inputs_neurons  = weights_inputs_neurons
        self.nsyns_inputs_neurons    = nsyns_inputs_neurons
        self.syn_inputs_neurons      = syn_inputs_neurons
        self.weights_neurons_neurons = weights_neurons_neurons
        self.nsyns_neurons_neurons   = nsyns_neurons_neurons
        self.syn_neurons_neurons     = syn_neurons_neurons
        self.weights_noise_neurons   = weights_noise_neurons
        self.delay_inputs_neurons    = delay_inputs_neurons
        self.delay_neurons_neurons   = delay_neurons_neurons
        self.record_mp = record_mp #new 
        self.resolution = resolution 

        self.ISI_noise = ISI_noise
        self.create_ca1_network = create_ca1_network
        self.record_noise = record_noise
        self.save_input_connection = save_input_connection
        self.i_bas = i_bas
        self.i_pyr=i_pyr
        self.delays_in_noise        = delays_in_noise
        self.delay_mean_noise       = delay_mean_noise
        self.delay_mean_inputs      = delay_mean_inputs
        self.delay_mean_connections = delay_mean_connections

        self.topology  = topology
        self.n_pyr_ca3 = n_pyr_ca3
        self.n_bas_ca3 = n_bas_ca3
        self.n_olm_ca3 = n_olm_ca3
        self.n_pyr_ca1 = n_pyr_ca1
        self.n_bas_ca1 = n_bas_ca1
        self.n_olm_ca1 = n_olm_ca1
        self.n_cck_ca1 = n_cck_ca1
        self.set_pyramidal_topology(type=self.topology) #this will change the number of neurons

        self.n_pop = [self.n_pyr_ca3, self.n_bas_ca3, self.n_olm_ca3]
        self.cell_labels = ["pyr_ca3", "bas_ca3", "olm_ca3"]
        if self.create_ca1_network:
            self.n_pop += [self.n_pyr_ca1, self.n_bas_ca1, self.n_olm_ca1, self.n_cck_ca1]
            self.cell_labels += ["pyr_ca1", "bas_ca1", "olm_ca1", "cck_ca1"]

        self.n_pop = np.cumsum(self.n_pop)
        self.n = self.n_pop[-1]

        self.set_gids()
        self.create_populations()

        self.bseed = bseed # seed for background noise should be ordened in order to keep same realization of when running ca3 and ca1 together
        self.wseed = wseed # seed for wiring
        self.iseed = iseed # seed for external inputs (EC2, EC3, MS and DG)
        self.jseed = jseed # wiring for external inputs
        self.dseed = 1993  # for distribution of delays on connections

        self.jitter = jitter # sigma of the distribution of delays
        self.DoMakeNoise = DoMakeNoise
        self.DoMakeExternalInputs = DoMakeExternalInputs
        self.MakeNetStim  = MakeNetStim
        self.MakeCellStim = MakeCellStim # added on 14/03/2022

        self.SCGain   = 0.5
        self.SaveConn = SaveConn

        self.time_initial = 50.0 # time in which external inputs are generated

        self.noise_burst = noise_burst

        self.inputs_folder         = inputs_folder
        # self.burst_level_label    = burst_level_label esta sobra
        self.external_inputs_iseed = external_inputs_iseed

        self.set_delays_distribution()
        if connections:
            print("Setting Connections")
            self.set_all_conns()

    def set_pyramidal_topology(self, type):
        dx = 50.0
        if type == "circular":
            npyr = np.array([6*i for i in range(1,17)])[2:]
            radius = np.array([dx*i for i in range(1,17)])[2:]
            self.pyr_ca3_x, self.pyr_ca3_y = [],[]
            for i,n in enumerate(npyr):
                theta = np.linspace(0,2*np.pi, n+1)[:-1]
                self.pyr_ca3_x.append(radius[i]*np.cos(theta))
                self.pyr_ca3_y.append(radius[i]*np.sin(theta))
            self.pyr_ca3_x = np.concatenate(self.pyr_ca3_x)
            self.pyr_ca3_y = np.concatenate(self.pyr_ca3_y)
            self.n_pyr_ca3 = np.sum(npyr)

            if self.create_ca1_network:
                self.pyr_ca1_x, self.pyr_ca1_y = [],[]
                for i,n in enumerate(npyr):
                    theta = np.linspace(0,2*np.pi, n+1)[:-1]
                    self.pyr_ca1_x.append(1e9+radius[i]*np.cos(theta))
                    self.pyr_ca1_y.append(radius[i]*np.sin(theta))
                self.pyr_ca1_x = np.concatenate(self.pyr_ca1_x)
                self.pyr_ca1_y = np.concatenate(self.pyr_ca1_y)
                self.n_pyr_ca1 = np.sum(npyr)

        elif type == "rectangular":
            n_y = int(np.round(np.sqrt(self.n_pyr_ca3/2)))
            n_x=2*n_y
            if(n_x*n_y!=self.n_pyr_ca3):
                raise Exception("@Network.set_pyramidal_topology(), number of pyrs does not fit in a n_x,2*n_x grid")
            
            self.pyr_ca3_x = np.concatenate([np.arange(0,n_x,1)*dx for i in range(n_y)])-n_x*dx/2
            self.pyr_ca3_y = np.concatenate([np.ones(n_x)*dx*i for i in range(n_y)])-n_y*dx/2
            
            if self.create_ca1_network:
                self.pyr_ca1_x = np.concatenate([np.arange(0,n_x,1)*dx for i in range(n_y)])-n_x*dx/2
                self.pyr_ca1_y = np.concatenate([np.ones(n_x)*dx*i for i in range(n_y)])-n_y*dx/2
                
        else:
            raise Exception("@Network.set_pyramidal_topology(), Not implemented row topology yet")

    def set_gids(self):
        """ set the gidlist on this host. """
        self.gidlist=[]
        for i in range(pc.id(), self.n ,pc.nhost()):
            self.gidlist.append(i)
            pc.set_gid2node(i, int(pc.id()) )

    def create_populations(self):
        gid_pyr_ca3, gid_bas_ca3, gid_olm_ca3 = [],[],[]
        gid_pyr_ca1, gid_bas_ca1, gid_olm_ca1,gid_cck_ca1 = [],[],[],[] #,[],[],[]

        for i in self.gidlist:
            if np.logical_and(i>= 0,i< self.n_pop[0]): gid_pyr_ca3.append(i)
            if np.logical_and(i>=self.n_pop[0],i<self.n_pop[1]): gid_bas_ca3.append(i)
            if np.logical_and(i>=self.n_pop[1],i<self.n_pop[2]): gid_olm_ca3.append(i)
            if self.create_ca1_network:
                if np.logical_and(i>=self.n_pop[2],i<self.n_pop[3]): gid_pyr_ca1.append(i)
                if np.logical_and(i>=self.n_pop[3],i<self.n_pop[4]): gid_bas_ca1.append(i)
                if np.logical_and(i>=self.n_pop[4],i<self.n_pop[5]): gid_olm_ca1.append(i)
                if np.logical_and(i>=self.n_pop[5],i<self.n_pop[6]): gid_cck_ca1.append(i) 

        # esto hay que retocarlo un poco.
        self.pyr_ca3 = Population(cell_type=PyrAdr_CA3, pgidlist=gid_pyr_ca3, n=self.n_pyr_ca3, x= 0, y=0, z=0, dx=50.0, amp= self.i_pyr, dur=1e9, delay=2*h.dt, nseg = self.nseg_pyr, naux = 0, pyr=True, x_position=self.pyr_ca3_x,y_position=self.pyr_ca3_y, record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # Estas corrientes se han anadido
        self.bas_ca3 = Population(cell_type=Bwb,        pgidlist=gid_bas_ca3, n=self.n_bas_ca3, x=10, y=0, z=0, dx=50.0, amp= self.i_bas,     dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[0],record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # para obtener el baseline state
        self.olm_ca3 = Population(cell_type=Ow,         pgidlist=gid_olm_ca3, n=self.n_olm_ca3, x=20, y=0, z=0, dx=50.0, amp=-25e-3, dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[1],record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # dur tiene que se mayor o igual que el tiempo de simulacion
        self.cell_ca3 = self.pyr_ca3.cell+self.bas_ca3.cell+self.olm_ca3.cell

        if self.create_ca1_network:
            self.pyr_ca1 = Population(cell_type=PyrAdr_CA1, pgidlist=gid_pyr_ca1, n=self.n_pyr_ca1, x=1e3+0,  y=0, z=0, dx=50, amp= 40e-3, dur=1e9, delay=2*h.dt, nseg = self.nseg_pyr, naux=self.n_pop[2],record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # Estas corrientes se han anadido
            self.bas_ca1 = Population(cell_type=Bwb,        pgidlist=gid_bas_ca1, n=self.n_bas_ca1, x=1e3+10, y=0, z=0, dx=50, amp= 0,     dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[3],record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # para obtener el baseline state
            self.olm_ca1 = Population(cell_type=Ow,         pgidlist=gid_olm_ca1, n=self.n_olm_ca1, x=1e3+20, y=0, z=0, dx=50, amp=-25e-3, dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[4],record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # dur tiene que se mayor o igual que el tiempo de simulacion
            self.cck_ca1 = Population(cell_type=Cck_cell,   pgidlist=gid_cck_ca1, n=self.n_cck_ca1, x=1e3+30, y=0, z=0, dx=50, amp= 10e-3, dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[5],record_mp = self.record_mp, resolution=self.resolution,record_spikes_dendr=self.record_spikes_dendr) # dur tiene que se mayor o igual que el tiempo de simulacion
            self.cell_ca1 = self.pyr_ca1.cell+self.bas_ca1.cell+self.olm_ca1.cell+self.cck_ca1.cell

    def set_noise_inputs(self,simdur):
        if self.DoMakeNoise:
            if self.MakeNetStim:
                self.make_all_NetStims(simdur)
            else:
                print("Simulation without external noide. It will be fixed in the future.")
        else:
            print("Simulation without external noise")

        if self.DoMakeExternalInputs:
            if self.MakeCellStim:
                print("Simulation without external inputs. It will be fixed in the future.")
                self.make_all_CellStims(simdur) # generated by Python
            else:
                self.load_all_CellStims(simdur) # generated by Python
        else:
            print("Simulation without external inputs")

        print("Done!")

    #both should be called @ beginning of each sim
    def init_NetStims(self):
        # h.mcell_ran4_init(self.iseed)
        for i in range(len(self.nrl)):
            rds = self.nrl[i]
            sead = self.nrlsead[i]
            rds.MCellRan4(sead,sead)
            rds.negexp(1)

    def init_CellStims(self):
        print('---------------------')
        print('INIT CELLSTIMS')
        print('---------------------')
        for i in range(len(self.nrl_)):
            rds = self.nrl_[i]
            rds.start()
            rds.set_generator()

    def make_NetStims(self, po, syn, delay, w, ISI, time_limit, add_to_sead): #added delay in 24th
        ''' noise for background '''
        po.nssidx[syn] = len(self.nsl) #index into net.nsl
        po.ncsidx[syn] = len(self.ncl) #index into net.
        self.netstims_tvec[po.cell[0].label+"_"+syn] = h.Vector()
        self.netstims_ivec[po.cell[0].label+"_"+syn] = h.Vector()
        for i in range(po.ng):
            cell = po.cell[i]

            ns = h.NetStim()
            ns.interval = ISI
            ns.noise = 1
            ns.number = (1e3/ISI)*time_limit
            ns.start = 0

            nc = h.NetCon(ns,cell.__dict__[syn].syn)
            #print(po.pgidlist[i]-po.naux,po.naux)
            nc.delay = delay[po.pgidlist[i]-po.naux]
            nc.weight[0] = w
            if pc.id() == 0:
                if self.record_noise:
                    nc.record(self.netstims_tvec[po.cell[0].label+"_"+syn], self.netstims_ivec[po.cell[0].label+"_"+syn],i)

            rds = h.Random()
            rds.negexp(1)            # set random # generator using negexp(1) - avg interval in NetStim
            sead = po.pgidlist[i] + self.bseed + add_to_sead #to avoid sead = 0 and having an element to control
            rds.MCellRan4(sead,sead) # seeds are in order, shouldn't matter
            ns.noiseFromRandom(rds)  # use random # generator for this NetStim

            #ns.start = rds.discunif(0,1e3) # start inputs random time btwn 0-1e3 ms to avoid artificial sync
            #rds.MCellRan4(sead,sead) # reinit rand # generator
            #nc_list.append(nc)
            self.nsl.append(ns)
            self.ncl.append(nc)
            cell.ncl.append(nc)
            self.nrl.append(rds) #esta lista de random es la que no entiendo muy bien
            self.nrlsead.append(sead)

        add_to_sead = add_to_sead + po.n

        po.nseidx[syn] = len(self.nsl)-1
        po.nceidx[syn] = len(self.ncl)-1
        return add_to_sead

    def set_delays_distribution(self):
        # I generate delays for every input and noise source
        self.sigma_noise  = 0.0
        if self.delays_in_noise:
            self.sigma_noise = 1.0

        rdtmp = self.dseed
        rdtmp = self.set_delays_to_synaptic_noise(delay_mean=self.delay_mean_noise, sigma=self.sigma_noise,seed=rdtmp)
        rdtmp = self.set_delays_to_synaptic_inputs(seed=rdtmp)
        rdtmp = self.set_delays_to_synpatic_connections(seed=rdtmp)

    def set_delays_to_synpatic_connections(self, seed):
        ''' delays neurons -> neurons connections '''
        print("Setting the delays: neurons to neurons")

        self.delays_neurons_neurons = dict.fromkeys(self.weights_neurons_neurons.keys() )
        for key in self.delays_neurons_neurons.keys():
            if key.endswith("ca3"):
                np.random.seed(seed)
                nsyn = self.nsyns_neurons_neurons[key]
                po = self.__dict__[key.split("_to_")[1]]
                self.delays_neurons_neurons[key] = [] 
                delay_mean = self.delay_neurons_neurons[key][0]
                delay_std = self.delay_neurons_neurons[key][1]
                for n in nsyn:
                    self.delays_neurons_neurons[key].append( np.abs(delay_mean+np.random.normal(0,delay_std,(n,po.n))) )
                    seed += 1

            if self.create_ca1_network & key.endswith("ca1"):
                np.random.seed(seed)
                nsyn = self.nsyns_neurons_neurons[key]
                po = self.__dict__[key.split("_to_")[1]]
                self.delays_neurons_neurons[key] = []
                delay_mean = self.delay_neurons_neurons[key][0]
                delay_std = self.delay_neurons_neurons[key][1]
                for n in nsyn:
                    self.delays_neurons_neurons[key].append( np.abs(delay_mean+np.random.normal(0,delay_std,(nsyn,po.n))) ) 
                    seed += 1
        return seed

    def set_delays_to_synaptic_inputs(self, seed):
        ''' delays inputs -> neurons connections '''
        print("Setting the delays: inputs to neurons")
        self.delays_inputs_neurons = dict.fromkeys(self.weights_inputs_neurons.keys())
        for key in self.delays_inputs_neurons.keys():
            if key.endswith("ca3"):
                np.random.seed(seed)
                nsyn = self.nsyns_inputs_neurons[key]
                po = self.__dict__[key.split("_to_")[1]]
                self.delays_inputs_neurons[key] = [] 
                delay_mean = self.delay_inputs_neurons[key][0]
                delay_std = self.delay_inputs_neurons[key][1]
                for n in nsyn:
                    self.delays_inputs_neurons[key].append( np.abs(delay_mean+np.random.normal(0,delay_std,(n,po.n))) )
                    seed += 1
            
            if self.create_ca1_network & key.endswith("ca1"):
                np.random.seed(seed)
                nsyn = self.nsyns_inputs_neurons[key]
                po = self.__dict__[key.split("_to_")[1]]
                self.delays_inputs_neurons[key] = []
                delay_mean = self.delay_inputs_neurons[key][0]
                delay_std = self.delay_inputs_neurons[key][1]
                for n in nsyn:
                    self.delays_inputs_neurons[key].append( np.abs(delay_mean+np.random.normal(0,delay_std,(n,po.n))) )
                    seed += 1

        return seed

    def set_delays_to_synaptic_noise(self,delay_mean,sigma,seed):
        ''' delays inputs noise background -> neurons '''
        print("Setting the delays: noise to neurons")
        self.delays_noise_neurons = dict.fromkeys(self.cell_labels)
        for key in self.cell_labels:
            np.random.seed(seed)
            po = self.__dict__[key]
            if key.split("_")[0] == "pyr":
                keys = ["somaAMPA_noise","Adend1NMDA_dgreg","Adend1AMPA_dgreg","somaGABA_noise","Adend3AMPA_noise","Adend3GABA_noise","Adend3NMDA_noise"]
            else:
                keys = self.weights_noise_neurons[key].keys()
            self.delays_noise_neurons[key] = dict.fromkeys(keys)
            for k in keys:
                self.delays_noise_neurons[key][k] = np.abs(delay_mean+sigma*np.random.normal(0.0,self.jitter,po.n))
            seed+=1
        return seed

    def make_all_NetStims(self,simdur): #,rdmseed):
        print("Making NetStims")
        self.netstims_tvec = {} # saving all
        self.netstims_ivec = {}
        self.nsl = [] # NetStim List
        self.ncl = [] # NetCon List
        self.nrl = [] # Random List for NetStims
        self.nrlsead = [] #List of seeds for NetStim random

        print("Making Background Noise")
        rdtmp = 0 # add_to_sead value - incremented in make_NetStims
        for cell in self.cell_labels:
            if cell.endswith("ca3"):
                print("to "+cell.split("_")[0]+" "+cell.split("_")[1])
                po = self.__dict__[cell]
                for syn in self.weights_noise_neurons[cell].keys():
                    delays = self.delays_noise_neurons[cell][syn]
                    weight = self.weights_noise_neurons[cell][syn]
                    rdtmp  = self.make_NetStims(po=po, syn=syn, delay=delays, w=weight,  ISI=self.ISI_noise,  time_limit=simdur, add_to_sead=rdtmp)#, sead=rdtmp)

            if self.create_ca1_network & cell.endswith("ca1"):
                print("to "+cell.split("_")[0]+" "+cell.split("_")[1])
                po = self.__dict__[cell]
                for syn in self.weights_noise_neurons[cell].keys():
                    delays = self.delays_noise_neurons[cell][syn]
                    weight = self.weights_noise_neurons[cell][syn]
                    rdtmp  = self.make_NetStims(po=po, syn=syn, delay=delays, w=weight,  ISI=self.ISI_noise,  time_limit=simdur, add_to_sead=rdtmp)#, sead=rdtmp)
    
    def set_CellStims_connections(self, ncellstims):
        print("set_CellStims_connections executed")
        keys = [w.split("_to_")[0] for w in self.weights_inputs_neurons]
        #keys = ["sep_180","sep_360","ec2_180","ec2_360","ec3_180","ec3_360","dg_regular","dg_burst"]
        inputs = dict.fromkeys(keys)
        for i, key in enumerate(keys):
            inputs[key] = ncellstims[i]

        random.seed(self.jseed) # for wiring
        keys = self.weights_inputs_neurons.keys()
        self.connectivity_inputs_neurons = dict.fromkeys(keys)
        for key in keys:
            if key.endswith("ca3"):
                ninputs = inputs[key.split("_to_")[0]]
                po = self.__dict__[key.split("_to_")[1]]
                self.connectivity_inputs_neurons[key] = []
                for n in self.nsyns_inputs_neurons[key]:
                    self.connectivity_inputs_neurons[key].append( self.set_connectivity(ninputs, po.n, n) ) 
                
            if self.create_ca1_network & key.endswith("ca1"):
                ninputs = inputs[key.split("_to_")[0]]
                po = self.__dict__[key.split("_to_")[1]]
                self.connectivity_inputs_neurons[key] = []
                for n in self.nsyns_inputs_neurons[key]:
                    self.connectivity_inputs_neurons[key].append( self.set_connectivity(ninputs, po.n, n) )
    
    def load_all_CellStims(self,simdur):
        self.ncl_ = []
        print("Loading CellsStims")
        keys = [w.split("_to_")[0] for w in self.weights_inputs_neurons]
        #keys = ['sep_180','sep_360','ec2_180','ec2_360','ec3_180','ec3_360','dg_regular','dg_burst']

        self.tvec = dict.fromkeys(keys)
        self.idvec = dict.fromkeys(keys)
        t0 = 0.0
        tf = simdur+t0
        ncellstims = []
        for key in keys:
            file = f"external_inputs_{key}.lzma"
            print(os.path.join(self.inputs_folder,file))
            data = file_management.load_lzma(os.path.join(self.inputs_folder,file))
            x = np.array( data["tvec"][data["iseed"] == self.external_inputs_iseed].values )
            y = np.array( data["idvec"][data["iseed"]== self.external_inputs_iseed].values )
            window = np.logical_and(x>=t0,x<=tf)
            self.tvec[key]   = x[window]-t0
            self.idvec[key]  = y[window]
            # ncellstims.append(np.max(self.idvec[key])+1)
            ncellstims.append( np.max(y)+1 )
            
        print("Setting the connections")
        self.set_CellStims_connections(ncellstims)

        print("Creating dictionary with all parameters")
        keys1 = self.weights_inputs_neurons.keys()
        keys2 = ["idv", "spikes", "population", "synapses", "delays", "weights", "connectivity"]
        self.external_inputs_data = dict.fromkeys(keys1)
        for key in keys1:
            if key.endswith("ca3"):
                self.external_inputs_data[key] = dict.fromkeys(keys2)
                inputs = key.split("_to_")[0]
                po = self.__dict__[key.split("_to_")[1]]
                self.external_inputs_data[key]["idv"]          = self.idvec[inputs]
                self.external_inputs_data[key]["spikes"]       = self.tvec[inputs]
                self.external_inputs_data[key]["population"]   = po
                self.external_inputs_data[key]["synapses"]     = self.syn_inputs_neurons[key]
                self.external_inputs_data[key]["delays"]       = self.delays_inputs_neurons[key]
                self.external_inputs_data[key]["weights"]      = self.weights_inputs_neurons[key]
                self.external_inputs_data[key]["connectivity"] = self.connectivity_inputs_neurons[key]

            if self.create_ca1_network & key.endswith("ca1"):
                self.external_inputs_data[key] = dict.fromkeys(keys2)
                inputs = key.split("_to_")[0]
                po = self.__dict__[key.split("_to_")[1]]
                self.external_inputs_data[key]["idv"]          = self.idvec[inputs]
                self.external_inputs_data[key]["spikes"]       = self.tvec[inputs]
                self.external_inputs_data[key]["population"]   = po
                self.external_inputs_data[key]["synapses"]     = self.syn_inputs_neurons[key]
                self.external_inputs_data[key]["delays"]       = self.delays_inputs_neurons[key]
                self.external_inputs_data[key]["weights"]      = self.weights_inputs_neurons[key]
                self.external_inputs_data[key]["connectivity"] = self.connectivity_inputs_neurons[key]
    
    def record_external_inputs(self, CellStim, tvec, idvec):
        nil = h.Section(name='nil')
        nil.insert('hh')
        syn_ = h.ExpSyn(nil(0.5))
        for i in range(len(CellStim)):
            nc = h.NetCon(CellStim[i],syn_)
            nc.record(tvec,idvec,i)

    def set_input_connection(self,ns_src,trg,syn_list,delay,w_list,conn, save=False, name='noname'):
        '''Important: syn_list and w_list must be lists.'''
        for k, conn_ in enumerate(conn):
            for post_id, all_pre in enumerate(conn_):
                for j, pre_id in enumerate(all_pre):
                    for syn,w in zip(syn_list[k],w_list[k]):
                        self.ncl_.append(h.NetCon(ns_src[pre_id], trg.cell[post_id].__dict__[syn].syn, 0, delay[j,post_id], w))

    def make_conn(self, preN, postN, conv):
        conn = np.zeros((postN,conv),dtype=np.int16)
        for i in range(postN):
            conn[i,:]=random.sample(range(preN),conv)
        return conn

    def set_connectivity(self, nsrc, ntrg, conv):
        conn = self.make_conn(nsrc,ntrg,conv)
        conn_ = conn[ np.arange(pc.id(), ntrg, pc.nhost()),:]
        return conn_

    def set_cells_connectivity(self):
        random.seed(self.wseed)
        ''' Set the connectivity '''
        keys = self.weights_neurons_neurons.keys()
        self.connectivity_neurons_neurons = dict.fromkeys(keys)
        for key in keys:
            if key.endswith("ca3"):
                print(key)
                from_po = self.__dict__[key.split("_to_")[0]]
                to_po   = self.__dict__[key.split("_to_")[1]]
                self.connectivity_neurons_neurons[key] = []
                for n in self.nsyns_neurons_neurons[key]:
                    self.connectivity_neurons_neurons[key].append( self.set_connectivity(from_po.n, to_po.n, n) )
            
            if self.create_ca1_network & key.endswith("ca1"):
                from_po = self.__dict__[key.split("_to_")[0]]
                to_po   = self.__dict__[key.split("_to_")[1]]
                self.connectivity_neurons_neurons[key] = [] 
                for n in self.nsyns_neurons_neurons[key]:
                    self.connectivity_neurons_neurons[key].append( self.set_connectivity(from_po.n, to_po.n, n) )

    def set_all_conns(self):
        self.set_cells_connectivity()
        print("Making connections among populations")
        print("------------------------------------")
        for key in self.weights_neurons_neurons.keys():
            if key.endswith("ca3"):
                # print(key)
                from_po = self.__dict__[key.split("_to_")[0]]
                to_po   = self.__dict__[key.split("_to_")[1]]
                delays  = self.delays_neurons_neurons[key]
                syn     = self.syn_neurons_neurons[key]
                weights = self.weights_neurons_neurons[key]
                conn    = self.connectivity_neurons_neurons[key]
                self.set_connections(from_po, to_po, syn, delays, weights, conn)
            if self.create_ca1_network & key.endswith("ca1"):
                from_po = self.__dict__[key.split("_to_")[0]]
                to_po   = self.__dict__[key.split("_to_")[1]]
                delays  = self.delays_neurons_neurons[key]
                syn     = self.syn_neurons_neurons[key]
                weights = self.weights_neurons_neurons[key]
                conn    = self.connectivity_neurons_neurons[key]
                self.set_connections(from_po, to_po, syn, delays, weights, conn)

    def set_conn_weight(self, conn, weight):
        for nc in conn:
            nc.weight[0] = weight

    def record_synaptic_currents(self):
        ''' save synaptic currents that every neuron receives '''
        for cell in self.cell_ca3: 
            cell.record_synapses()
        if self.create_ca1_network:
            for cell in self.cell_ca1:
                cell.record_synapses()

    def set_connections(self,src,trg,syn_list,delay,w_list,conn, print_=False): # mirar esto con detalle
        #conn = self.make_conn(src.n,trg.ng,conv, trg)
        #conn = self.make_conn(src.n,trg.n,conv, trg)
        #conn = conn[ np.arange(pc.id(), trg.n, pc.nhost()),:]
        for k, conn_ in enumerate(conn):
            for post_id, all_pre in enumerate(conn_):
                #trg_id = trg.pgidlist[post_id]
                #delays = delay[:,post_id]
                #print(post_id,delays)
                for j, pre_id in enumerate(all_pre):
                    src_id = pre_id+src.naux
                    for syn,w in zip(syn_list[k],w_list[k]):
                        #trg.nc.append(pc.gid_connect(src_id, pc.gid2cell(trg_id).__dict__[syn].syn))
                        trg.nc.append(pc.gid_connect(pre_id+src.naux, trg.cell[post_id].__dict__[syn].syn))
                        trg.nc[-1].weight[0] = w
                        trg.nc[-1].delay = delay[k][j,post_id]
                        trg.nc[-1].threshold = 0.0

        if self.SaveConn: # we do not use it, but I left it here just in case
            try:
                print(self.nqcon.size())
            except:
                self.nqcon = h.NQS("id1","id2","w","syn")
                self.nqcon.strdec("syn")
            for post_id, all_pre in enumerate(conn):
                for j, pre_id in enumerate(all_pre):
                    self.nqcon.append(src.cell[pre_id].id,trg.cell[post_id].id,w,syn)

    def calc_lfp(self): # don't use it, but again, I left it here just in case
        # lfp is modeled as a difference between voltages in distal apical and basal compartemnts
        self.vlfp_ca3 = h.Vector(self.pyr_ca3.cell[0].Adend3_volt.size()) #lfp in neuron Vector
        for cell in self.pyr_ca3.cell:
            self.vlfp_ca3.add(cell.Adend3_volt)
            self.vlfp_ca3.sub(cell.Bdend_volt)
        #self.vlfp_ca3.div(len(self.pyr_ca3.cell)) # normalize lfp by amount of pyr cells
        self.vlfp_ca3=np.array(self.vlfp_ca3.to_python())

        if self.create_ca1_network:
            self.vlfp_ca1 = h.Vector(self.pyr_ca1.cell[0].Adend3_volt.size())
            for cell in self.pyr_ca1.cell:
                self.vlfp_ca1.add(cell.Adend3_volt)
                self.vlfp_ca1.sub(cell.Bdend_volt)
            #self.vlfp_ca1.div(len(self.pyr_ca1.cell)) # normalize lfp by amount of pyr cells
            self.vlfp_ca1=np.array(self.vlfp_ca1.to_python())

    def make_all_CellStims(self,simdur):
        # to tired to transcribe this function 
        # it has to replicate what it is done int external_input.py
        pass 

        # '''
        # This function load the external inputs created by make_external_inputs_python
        # '''
        # self.ncl_ =[]
        # print("Making external inputs")
        # self.idvec, self.tvec, self.nsl_n, labels, self.input_time = make_external_inputs_python(simdur, self.iseed, pburst=0.5)

        # print("Making cellstims connectionss")
        # self.set_CellStims_connections(self.nsl_n)

        # print("Creating dictionary with all parameters")
        # keys1 = self.weights_inputs_neurons.keys()
        # keys2 = ["idv", "spikes", "population", "synapses", "delays", "weights", "connectivity"]
        # self.external_inputs_data = dict.fromkeys(keys1)
        # for key in keys1:
        #     if key.endswith("ca3"):
        #         self.external_inputs_data[key] = dict.fromkeys(keys2)
        #         inputs = key.split("_to_")[0]
        #         po = self.__dict__[key.split("_to_")[1]]
        #         self.external_inputs_data[key]["idv"]          = self.idvec[inputs]
        #         self.external_inputs_data[key]["spikes"]       = self.tvec[inputs]
        #         self.external_inputs_data[key]["population"]   = po
        #         self.external_inputs_data[key]["synapses"]     = self.syn_inputs_neurons[key]
        #         self.external_inputs_data[key]["delays"]       = self.delays_inputs_neurons[key]
        #         self.external_inputs_data[key]["weights"]      = self.weights_inputs_neurons[key]
        #         self.external_inputs_data[key]["connectivity"] = self.connectivity_inputs_neurons[key]

        #     if self.create_ca1_network & key.endswith("ca1"):
        #         self.external_inputs_data[key] = dict.fromkeys(keys2)
        #         inputs = key.split("_to_")[0]
        #         po = self.__dict__[key.split("_to_")[1]]
        #         self.external_inputs_data[key]["idv"]          = self.idvec[inputs]
        #         self.external_inputs_data[key]["spikes"]       = self.tvec[inputs]
        #         self.external_inputs_data[key]["population"]   = po
        #         self.external_inputs_data[key]["synapses"]     = self.syn_inputs_neurons[key]
        #         self.external_inputs_data[key]["delays"]       = self.delays_inputs_neurons[key]
        #         self.external_inputs_data[key]["weights"]      = self.weights_inputs_neurons[key]
        #         self.external_inputs_data[key]["connectivity"] = self.connectivity_inputs_neurons[key]