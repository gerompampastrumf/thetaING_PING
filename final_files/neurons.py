"""
03/09/2021

ultimo intento para hacer que funcione el pc = h.ParallelContext

Last modification: creation of a synpase in Adend3 in pyr CA3 and CA1
"""
###############################################################################
#
# Synapases and general cell
#
###############################################################################
from neuron import h,gui
import numpy as np
h.load_file('stdrun.hoc')
#os.system('nrnivmodl')

class Synapse:
    def __init__(self, sect, loc, tau1, tau2, e):
        self.syn      = h.MyExp2SynBB(loc, sec=sect)
        self.syn.tau1 = tau1
        self.syn.tau2 = tau2
        self.syn.e    = e
        '''
        self.syn= h.ExpSyn(loc, sec=sect)#h.MyExp2SynBB(loc, sec=sect)
        self.syn.tau= 2.0+
        '''
class SynapseNMDA:
    def __init__(self, sect, loc, tau1, tau2, tau1NMDA, tau2NMDA, r, e):
        self.syn            = h.MyExp2SynNMDABB_b(loc, sec=sect)
        self.syn.tau1       = tau1
        self.syn.tau2       = tau2
        self.syn.tau1NMDA   = tau1NMDA
        self.syn.tau2NMDA   = tau2NMDA
        self.syn.r          = r
        self.syn.e          = e

class Cell:
    "General cell"
    def __init__(self,x,y,z,id, nseg, record, resolution=1.0,record_spikes_dendr={"Adend1":None,"Adend2":None,"Adend3":None,"Bdend":None}): # new variable, record
        self.x=x
        self.y=y
        self.z=z
        self.id=id
        self.all_sec = []
        self.syn_list = []
        self.sec_list = []
        self.nseg = nseg
        self.record = record
        self.resolution = resolution
        self.record_spikes_dendr=record_spikes_dendr
        self.set_label()
        self.add_comp('soma', self.record["soma"])
        self.dendr_spike_detector= []
        self.set_morphology()
        self.set_conductances()
        self.set_synapses()
        self.set_inj()
        self.calc_area()

        self.spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_detector.threshold = 0.0
        self.spike_times = h.Vector()
        self.spike_detector.record(self.spike_times)
        self.nc  = [] # connection between population src-> this cell
        self.ncl = [] # external noise netcons
        self.ncs = [] # artificial neurons netcons
    def __repr__(self):
        return '{}[{}]'.format(self.name, self.id)

    def set_label(self): # added 11/03/2022
        pass

    def set_morphology(self):
        pass

    def set_conductances(self):
        pass

    def set_synapses(self):
        pass

    def record_synapses(self): # added 06/07/2021

        for name in self.syn_list:
            if "NMDA" not in name:
                self.__dict__["i"+name] = h.Vector()
                self.__dict__["i"+name].record( self.__dict__[name].syn._ref_i, self.resolution)
            else: # iNMDA synapses consits of a sum of an AMPA and NMDA contributions. Why? don't know
                self.__dict__["i"+name] = h.Vector()
                self.__dict__["i"+name].record( self.__dict__[name].syn._ref_iNMDA, self.resolution)

    def set_inj(self):
        self.somaInj = h.IClamp(0.5, sec=self.soma)

    def add_comp(self, name, rec,rec_spike_thres=None):
        #self.__dict__[name] = h.Section()
        self.__dict__[name] = h.Section(name=name, cell=self) # modify 31/08/2021 important to make pc.gid2cell() work
        # this issue has been solve via:  https://github.com/ahwillia/PyNeuron-Toolbox/issues/10
        # https://github.com/ahwillia/PyNeuron-Toolbox/commit/36ada2da6b2d508a3e7f054c40239044103a2234

        self.all_sec.append(self.__dict__[name])
        self.sec_list.append(name)
        # Record voltage
        if rec:
            self.__dict__[name+"_volt"] = h.Vector()
            self.__dict__[name+"_volt"].record(self.__dict__[name](0.5)._ref_v,self.resolution)
        if rec_spike_thres is not None:
            self.dendr_spike_detector.append(h.NetCon(self.__dict__[name](0.5)._ref_v, None, sec=self.__dict__[name]))
            self.dendr_spike_detector[-1].threshold = rec_spike_thres
            self.__dict__[name+"_spike_times"] = h.Vector()
            self.dendr_spike_detector[-1].record(self.__dict__[name+"_spike_times"])            

    # def plot_volt(self, name, fig=1):
    #     figure(fig)
    #     volt = self.__dict__[name+"_volt"].to_python()
    #     plot(arange(len(volt))*h.dt, volt)

    def calc_area(self):
        self.total_area = 0
        self.n = 0
        for sect in self.all_sec:
            self.total_area += h.area(0.5,sec=sect)
            self.n+=1

###############################################################################
#
# Basket Cell -- Bwb
#
###############################################################################

class Bwb(Cell):
    name = "Basket cell"
    def set_label(self):
        self.label='bas'

    def set_morphology(self):
        total_area = 10000 # um2
        self.soma.nseg  = self.nseg
        self.soma.cm    = 1      # uF/cm2
        diam = np.sqrt(total_area) # um
        L  = diam/np.pi  # um
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)

    def set_conductances(self):
        self.soma.insert('pas')   # leak current
        self.soma.e_pas = -65     # mV
        self.soma.g_pas = 0.1e-3  # S/cm2

        self.soma.insert('Nafbwb')
        self.soma.insert('Kdrbwb')

    def set_synapses(self):
        # external noise 
        self.somaAMPA_noise   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)   
        self.somaGABA_noise   = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
        # external inputs CA3        
        self.somaAMPA_ec2180  = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0) 
        self.somaNMDA_ec2180  = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaAMPA_ec2360  = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)                                         # not initially used
        self.somaNMDA_ec2360  = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)  # not initially used
        self.somaAMPA_dgreg   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.somaNMDA_dgreg   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaAMPA_dgburst = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0) 
        self.somaNMDA_dgburst = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        # external inputs in CA1
        self.somaAMPA_ec3180  = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.somaNMDA_ec3180  = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaAMPA_ec3360  = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.somaNMDA_ec3360  = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaAMPA_pyrCA3  = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.somaNMDA_pyrCA3  = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        # external inputs common 
        self.somaGABA_olm  = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        self.somaGABA_sep360  = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        # connections 
        self.somaGABA_bas     = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        self.somaAMPA_pyr     = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
        self.somaNMDA_pyr     = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0) 
       
        # self.somaGABAss = Synapse(sect=self.soma, loc=0.5, tau1=20.0, tau2=40.0, e=-80)# Synapse(sect=self.soma, loc=0.5, tau1=20.0, tau2=40.0, e=-80) # only for septal input # ojo!
        # self.somaGABAss = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        # self.somaGABAss = Synapse(sect=self.soma, loc=0.5, tau1=0,07, tau2=9.1, e=-80) # only for septal inpu
        # self.somaNMDA   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        # self.syn_list   = ['somaAMPAf','somaGABAf','somaGABAss','somaNMDA']

###############################################################################
#
# OLM cell -- Ow
#
###############################################################################
class Ow(Cell):
    name = "OLM cell"
    def set_label(self):
        self.label = "olm"

    def set_morphology(self):
        total_area = 10000 # um2
        self.soma.nseg  = self.nseg
        self.soma.cm    = 1      # uF/cm2
        diam = np.sqrt(total_area) # um
        L    = diam/np.pi  # um
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)

    def set_conductances(self):
        self.soma.insert('pas')
        self.soma.e_pas = -65     # mV
        self.soma.g_pas = 0.1e-3  # S/cm2

        self.soma.insert('Nafbwb')
        self.soma.insert('Kdrbwb')

        self.soma.insert('Iholmw')
        self.soma.insert('Caolmw')
        self.soma.insert('ICaolmw')
        self.soma.insert('KCaolmw')

    def set_synapses(self):
        # external noise 
        self.somaGABA_noise  = Synapse( sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80.0)
        self.somaAMPA_noise  = Synapse( sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # external inputs (both CA3 and CA1) 
        self.somaGABA_bas  = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        self.somaGABA_sep360  = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        # connections
        self.somaAMPA_pyr = Synapse( sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        self.somaNMDA_pyr = SynapseNMDA( sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)

        # self.somaGABAss = Synapse(    sect=self.soma, loc=0.5, tau1=0.07,   tau2=9.1, e=-80.0) # ojo! #Synapse(    sect=self.soma, loc=0.5, tau1=20.0,   tau2=40.0, e=-80.0) # only for septal input
        # self.somaGABAss = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        # self.somaNMDA   = SynapseNMDA( sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.syn_list   = ['somaGABAf','somaAMPAf','somaGABAss','somaNMDA']

###############################################################################
#
# CA3 Pyramidal Cell -- PyrAdr_CA3
#
###############################################################################

class PyrAdr_CA3(Cell):
    name = "Pyramidal CA3 cell"
    def set_label(self):
        self.label="pyr"

    def set_morphology(self):
        self.add_comp('Bdend', self.record["Bdend"],self.record_spikes_dendr["Bdend"])
        self.add_comp('Adend1',self.record["Adend1"],self.record_spikes_dendr["Adend1"])
        self.add_comp('Adend2',self.record["Adend2"],self.record_spikes_dendr["Adend2"])
        self.add_comp('Adend3',self.record["Adend3"],self.record_spikes_dendr["Adend3"])

        self.soma.nseg = self.nseg 
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,          20, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+20,       20, sec=self.soma)
        
        self.Bdend.nseg = self.nseg
        h.pt3dclear(sec=self.Bdend)
        h.pt3dadd(self.x, self.y, self.z,          2, sec=self.Bdend)
        h.pt3dadd(self.x, self.y, self.z-200,      2, sec=self.Bdend)
        
        self.Adend1.nseg = self.nseg
        h.pt3dclear(sec=self.Adend1)
        h.pt3dadd(self.x, self.y, self.z+20,       2, sec=self.Adend1)
        h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend1)
        
        self.Adend2.nseg = self.nseg
        h.pt3dclear(sec=self.Adend2)
        h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend2)
        h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend2)
        
        self.Adend3.nseg = self.nseg
        h.pt3dclear(sec=self.Adend3)
        h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend3)
        h.pt3dadd(self.x, self.y, self.z+20+150*3, 2, sec=self.Adend3)

        self.Bdend.connect(self.soma(0),      1)
        # self.soma.connect(self.Bdend,      1, 0)
        self.Adend1.connect(self.soma(1),     0)
        self.Adend2.connect(self.Adend1(1),   0)
        self.Adend3.connect(self.Adend2(1),   0)

    def set_conductances(self):
        for sect in self.all_sec:
            sect.insert('pas')
            sect.insert('nacurrent')
            sect.insert('kacurrent')
            sect.insert('kdrcurrent')
            sect.insert('hcurrent')
            sect.cm = 1 # Capacitancia
            sect.Ra = 150 #kOhm/cm2 # Citoplasmatic resistivity
            for seg in sect:
                seg.pas.g = 0.0000357 # 1/Rm (Rm=28 kOhm cm2)
                seg.pas.e = -70 # mV

        for seg in self.Adend1:
            seg.nacurrent.ki = 0.5
            seg.kacurrent.g  = 0.072
            # seg.hcurrent.v50 = -82
            seg.hcurrent.g   = 0.0002

        for seg in self.Adend2:
            seg.nacurrent.ki = 0.5
            seg.kacurrent.g  = 0
            seg.kacurrent.gd = 0.120
            seg.hcurrent.v50 = -90
            seg.hcurrent.g   = 0.0004

        for seg in self.Adend3:
            # to account for spunes and smaller dendritic branches, Cm was increased
            # and Rm decreased by a factor of 2 in the distal apical dendritic compartment
            seg.cm           = 2 # 2 * capacitance
            seg.pas.g        = 0.0000714 #(1/(Rm/2))
            seg.nacurrent.ki = 0.5
            seg.kacurrent.g  = 0
            seg.kacurrent.gd = 0.200
            seg.hcurrent.v50 = -90
            seg.hcurrent.g   = 0.0007
            seg.nacurrent.ki = 1

    def set_synapses(self):
        # external noise 
        self.somaAMPA_noise    = Synapse( sect=self.soma,   loc=0.3, tau1=0.05, tau2=5.3,  e=0 )   
        self.somaGABA_noise    = Synapse( sect=self.soma,   loc=0.7, tau1=0.07, tau2=9.1,  e=-80 ) 
        self.Adend3AMPA_noise  = Synapse( sect=self.Adend3, loc=0.3, tau1=0.05, tau2=5.3,  e=0 )   
        self.Adend3GABA_noise  = Synapse( sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80 ) 
        self.Adend3NMDA_noise  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # external inputs CA3        
        self.Adend3AMPA_ec2180  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 ) 
        self.Adend3NMDA_ec2180  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        self.Adend3AMPA_ec2360  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )                                         # not initially used
        self.Adend3NMDA_ec2360  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )  # not initially used
        self.Adend1AMPA_dgreg   = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )
        self.Adend1NMDA_dgreg   = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        self.Adend1AMPA_dgburst = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0 ) 
        self.Adend1NMDA_dgburst = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # external inputs in CA1 (assiming this as the pyr CA1) s
        self.Adend3AMPA_ec3180  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )                                         # not initially used
        self.Adend3NMDA_ec3180  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )  # not initially used
        self.Adend3AMPA_ec3360  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )
        self.Adend3NMDA_ec3360  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        self.Adend1AMPA_pyrCA3  = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )
        self.Adend1NMDA_pyrCA3  = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # connections in CA3
        self.somaGABA_bas      = Synapse(      sect=self.soma,   loc=0.7, tau1=0.07, tau2=9.1, e=-80 )
        self.Adend3GABA_olm    = Synapse(      sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1, e=-80 )
        self.BdendAMPA_pyr     = Synapse(      sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, e=0   )
        self.BdendNMDA_pyr     = SynapseNMDA(  sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # connections in CA1 (same as previous)
        self.somaGABA_cck      = Synapse(      sect=self.soma,   loc=0.7, tau1=0.07, tau2=9.1, e=-80 )
        self.Adend2GABA_cck    = Synapse(      sect=self.Adend2, loc=0.5, tau1=0.07, tau2=9.1, e=-80 )
      
        # self.somaGABAf   = Synapse(    sect=self.soma,   loc=0.5, tau1=0.07, tau2=9.1, e=-80.0) 
        # self.somaGABAfb  = Synapse(    sect=self.soma,   loc=0.3, tau1=0.07, tau2=9.1, e=-80.0) 
        # self.somaAMPAf   = Synapse(    sect=self.soma,   loc=0.7, tau1=0.05, tau2=5.3, e=0.0)
        
        # self.BdendAMPAf  = Synapse(     sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # self.BdendNMDA   = SynapseNMDA( sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.Adend1AMPAf = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # self.Adend1NMDA  = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.Adend2GABAf = Synapse(     sect=self.Adend2, loc=0.5, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3GABAf = Synapse(     sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3GABAfb = Synapse(    sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3AMPAf = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0.0)
        # self.Adend3NMDA  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
      
        ### Synapses onto Adend1
        # self.syn_list = ['Adend3AMPAf','Adend3NMDA','Adend3GABAf', 'Adend2GABAf','Adend1AMPAf','Adend1NMDA',
        #                 'somaAMPAf','somaGABAf','somaGABAfb','BdendAMPAf','BdendNMDA']


###############################################################################
#
# CA3 Pyramidal Cell but with less axial resistance and A-type K+ currents -- PyrAdr_CA3
#
###############################################################################

class PyrAdr_CA3_moreExc(Cell):
    name = "Pyramidal CA3 cell"
    def set_label(self):
        self.label="pyr"

    def set_morphology(self):
        self.add_comp('Bdend', self.record["Bdend"])
        self.add_comp('Adend1',self.record["Adend1"])
        self.add_comp('Adend2',self.record["Adend2"])
        self.add_comp('Adend3',self.record["Adend3"])

        self.soma.nseg = self.nseg 
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,          20, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+20,       20, sec=self.soma)
        
        self.Bdend.nseg = self.nseg
        h.pt3dclear(sec=self.Bdend)
        h.pt3dadd(self.x, self.y, self.z,          2, sec=self.Bdend)
        h.pt3dadd(self.x, self.y, self.z-200,      2, sec=self.Bdend)
        
        self.Adend1.nseg = self.nseg
        h.pt3dclear(sec=self.Adend1)
        h.pt3dadd(self.x, self.y, self.z+20,       2, sec=self.Adend1)
        h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend1)
        
        self.Adend2.nseg = self.nseg
        h.pt3dclear(sec=self.Adend2)
        h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend2)
        h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend2)
        
        self.Adend3.nseg = self.nseg
        h.pt3dclear(sec=self.Adend3)
        h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend3)
        h.pt3dadd(self.x, self.y, self.z+20+150*3, 2, sec=self.Adend3)

        self.Bdend.connect(self.soma(0),      1)
        # self.soma.connect(self.Bdend,      1, 0)
        self.Adend1.connect(self.soma(1),     0)
        self.Adend2.connect(self.Adend1(1),   0)
        self.Adend3.connect(self.Adend2(1),   0)

    def set_conductances(self):
        for sect in self.all_sec:
            sect.insert('pas')
            sect.insert('nacurrent')
            sect.insert('kacurrent')
            sect.insert('kdrcurrent')
            sect.insert('hcurrent')
            sect.cm = 1 # Capacitancia
            sect.Ra = 150.0/3.0 #kOhm/cm2 # Citoplasmatic resistivity
            for seg in sect:
                seg.pas.g = 0.0000357 # 1/Rm (Rm=28 kOhm cm2)
                seg.pas.e = -70 # mV
        for seg in self.soma:
            seg.kacurrent.g = 0.048/3.0

        for seg in self.Adend1:
            seg.nacurrent.ki = 0.5
            seg.kacurrent.g  = 0.072
            # seg.hcurrent.v50 = -82
            seg.hcurrent.g   = 0.0002

        for seg in self.Adend2:
            seg.nacurrent.ki = 0.5
            seg.kacurrent.g  = 0
            seg.kacurrent.gd = 0.120
            seg.hcurrent.v50 = -90
            seg.hcurrent.g   = 0.0004

        for seg in self.Adend3:
            # to account for spunes and smaller dendritic branches, Cm was increased
            # and Rm decreased by a factor of 2 in the distal apical dendritic compartment
            seg.cm           = 2 # 2 * capacitance
            seg.pas.g        = 0.0000714 #(1/(Rm/2))
            seg.nacurrent.ki = 0.5
            seg.kacurrent.g  = 0
            seg.kacurrent.gd = 0.200
            seg.hcurrent.v50 = -90
            seg.hcurrent.g   = 0.0007
            seg.nacurrent.ki = 1

    def set_synapses(self):
        # external noise 
        self.somaAMPA_noise    = Synapse( sect=self.soma,   loc=0.3, tau1=0.05, tau2=5.3,  e=0 )   
        self.somaGABA_noise    = Synapse( sect=self.soma,   loc=0.7, tau1=0.07, tau2=9.1,  e=-80 ) 
        self.Adend3AMPA_noise  = Synapse( sect=self.Adend3, loc=0.3, tau1=0.05, tau2=5.3,  e=0 )   
        self.Adend3GABA_noise  = Synapse( sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80 ) 
        
        # external inputs CA3        
        self.Adend3AMPA_ec2180  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 ) 
        self.Adend3NMDA_ec2180  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        self.Adend3AMPA_ec2360  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )                                         # not initially used
        self.Adend3NMDA_ec2360  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )  # not initially used
        self.Adend1AMPA_dgreg   = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )
        self.Adend1NMDA_dgreg   = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        self.Adend1AMPA_dgburst = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0 ) 
        self.Adend1NMDA_dgburst = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # external inputs in CA1 (assiming this as the pyr CA1) s
        self.Adend3AMPA_ec3180  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )                                         # not initially used
        self.Adend3NMDA_ec3180  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )  # not initially used
        self.Adend3AMPA_ec3360  = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )
        self.Adend3NMDA_ec3360  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        self.Adend1AMPA_pyrCA3  = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0 )
        self.Adend1NMDA_pyrCA3  = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # connections in CA3
        self.somaGABA_bas      = Synapse(      sect=self.soma,   loc=0.7, tau1=0.07, tau2=9.1, e=-80 )
        self.Adend3GABA_olm    = Synapse(      sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1, e=-80 )
        self.BdendAMPA_pyr     = Synapse(      sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, e=0   )
        self.BdendNMDA_pyr     = SynapseNMDA(  sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0 )
        
        # connections in CA1 (same as previous)
        self.somaGABA_cck      = Synapse(      sect=self.soma,   loc=0.7, tau1=0.07, tau2=9.1, e=-80 )
        self.Adend2GABA_cck    = Synapse(      sect=self.Adend2, loc=0.5, tau1=0.07, tau2=9.1, e=-80 )
      
        # self.somaGABAf   = Synapse(    sect=self.soma,   loc=0.5, tau1=0.07, tau2=9.1, e=-80.0) 
        # self.somaGABAfb  = Synapse(    sect=self.soma,   loc=0.3, tau1=0.07, tau2=9.1, e=-80.0) 
        # self.somaAMPAf   = Synapse(    sect=self.soma,   loc=0.7, tau1=0.05, tau2=5.3, e=0.0)
        
        # self.BdendAMPAf  = Synapse(     sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # self.BdendNMDA   = SynapseNMDA( sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.Adend1AMPAf = Synapse(     sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # self.Adend1NMDA  = SynapseNMDA( sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.Adend2GABAf = Synapse(     sect=self.Adend2, loc=0.5, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3GABAf = Synapse(     sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3GABAfb = Synapse(    sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3AMPAf = Synapse(     sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0.0)
        # self.Adend3NMDA  = SynapseNMDA( sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
      
        ### Synapses onto Adend1
        # self.syn_list = ['Adend3AMPAf','Adend3NMDA','Adend3GABAf', 'Adend2GABAf','Adend1AMPAf','Adend1NMDA',
        #                 'somaAMPAf','somaGABAf','somaGABAfb','BdendAMPAf','BdendNMDA']


###############################################################################
#
# CA1 Pyramidal Cell -- PyrAdr_CA1
#
###############################################################################

class PyrAdr_CA1(Cell):
    name = "Pyramidal CA1 cell"
    def set_label(self):
        self.label = "pyr"

    def set_morphology(self):
        self.add_comp('Bdend', self.record["Bdend"])
        self.add_comp('Adend1',self.record["Adend1"])
        self.add_comp('Adend2',self.record["Adend2"])
        self.add_comp('Adend3',self.record["Adend3"])

        self.soma.nseg = self.nseg
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,          20, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+20,       20, sec=self.soma)

        self.Bdend.nseg = self.nseg
        h.pt3dclear(sec=self.Bdend)
        h.pt3dadd(self.x, self.y, self.z,          2, sec=self.Bdend)
        h.pt3dadd(self.x, self.y, self.z-200,      2, sec=self.Bdend)
        
        self.Adend1.nseg = self.nseg
        h.pt3dclear(sec=self.Adend1)
        h.pt3dadd(self.x, self.y, self.z+20,       2, sec=self.Adend1)
        h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend1)
        
        self.Adend2.nseg = self.nseg
        h.pt3dclear(sec=self.Adend2)
        h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend2)
        h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend2)
        
        self.Adend3.nseg = self.nseg
        h.pt3dclear(sec=self.Adend3)
        h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend3)
        h.pt3dadd(self.x, self.y, self.z+20+150*3, 2, sec=self.Adend3)

        #self.Bdend.connect(self.soma,      0, 0)
        #self.Adend1.connect(self.soma,   0.5, 0)
        self.Bdend.connect(self.soma,      0, 1)
        self.Adend1.connect(self.soma,     1, 0)
        self.Adend2.connect(self.Adend1,   1, 0)
        self.Adend3.connect(self.Adend2,   1, 0)

    def set_conductances(self):
        for sect in self.all_sec:
            sect.insert('pas')
            sect.insert('nacurrent')
            sect.insert('kacurrent')
            sect.insert('kdrcurrent')
            sect.insert('hcurrent')
            sect.cm = 1 #Capacitancia
            sect.Ra = 150 #citoplasmatic resistivity
            for seg in sect: 
                seg.pas.g = 0.0000357 #1/Rm (Rm=28 kOhm cm2)
                seg.pas.e = -70 # mV
               
        for sec in self.Adend1:
            sec.nacurrent.ki = 0.5
            sec.kacurrent.g  = 0.072

        for sec in self.Adend2:
            sec.nacurrent.ki = 0.5
            sec.kacurrent.g  = 0
            sec.kacurrent.gd = 0.120

        for sec in self.Adend3:
            #to account for spunes and smaller dendritic branches, Cm was increased
            #and Rm decreased by a factor of 2 in the distal apical dendritic compartment
            sec.cm           = 2 # 2 * capacitance
            sec.pas.g        = 0.0000714 #(1/(Rm/2))
            sec.nacurrent.ki = 0.5
            sec.kacurrent.g  = 0
            sec.kacurrent.gd = 0.200
        
        for sec in self.Bdend:
            sec.nacurrent.ki  = 1
        
        # What makes the difference between CA3 y CA1, the expressing of protein HCN2/HCN1
        for sec in self.Adend1:
            sec.hcurrent.g = 0
            sec.hcurrent.gd = 0.0002
        for sec in self.Adend2: 
            sec.hcurrent.g = 0
            sec.hcurrent.gd = 0.0004
        for sec in self.Adend3:
            sec.hcurrent.g = 0
            sec.hcurrent.gd = 0.0007
        for sec in self.Bdend:
            sec.hcurrent.g  = 0
            sec.hcurrent.gd  = 0.0001
        for sec in self.soma:
            sec.hcurrent.g   = 0
            sec.hcurrent.gd   = 0.0001

    def set_synapses(self):
        # external noise 
        self.somaAMPA_noise    = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)   
        self.somaGABA_noise    = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
        self.Adend3AMPA_noise  = Synapse(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0)   
        self.Adend3GABA_noise  = Synapse(sect=self.Adend3, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
        # external inputs
        self.Adend3AMPA_ec3180  = Synapse(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0)                                         # not initially used
        self.Adend3NMDA_ec3180  = SynapseNMDA(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)  # not initially used
        self.Adend3AMPA_ec3360  = Synapse(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.Adend3NMDA_ec3360  = SynapseNMDA(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.Adend1AMPA_pyrCA3  = Synapse(sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.Adend1NMDA_pyrCA3  = SynapseNMDA(sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        # connections
        self.somaGABA_bas      = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        self.Adend3GABA_olm    = Synapse(sect=self.Adend3, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        self.BdendAMPA_pyr     = Synapse(sect=self.Bdend, loc=0.5, tau1=0.05, tau2=5.3, e=0)
        self.BdendNMDA_pyr     = SynapseNMDA(sect=self.Bdend, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaGABA_cck      = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
        self.Adend2GABA_cck    = Synapse(sect=self.Adend2, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
      
        # self.somaGABAf   = Synapse(    sect=self.soma,   loc=0.5, tau1=0.07, tau2=9.1, e=-80.0) 
        # self.somaGABAfb  = Synapse(    sect=self.soma,   loc=0.3, tau1=0.07, tau2=9.1, e=-80.0) 
        # self.somaAMPAf   = Synapse(    sect=self.soma,   loc=0.7, tau1=0.05, tau2=5.3, e=0.0)
        # self.BdendAMPAf  = Synapse(    sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # self.BdendNMDA   = SynapseNMDA(sect=self.Bdend,  loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.Adend1AMPAf = Synapse(    sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)
        # self.Adend1NMDA  = SynapseNMDA(sect=self.Adend1, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        # self.Adend2GABAf = Synapse(    sect=self.Adend2, loc=0.5, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3GABAf = Synapse(    sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80.0)
        # self.Adend3GABAfb = Synapse(    sect=self.Adend3, loc=0.7, tau1=0.07, tau2=9.1,  e=-80.0)

        # self.Adend3AMPAf = Synapse(    sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,  e=0.0)
        # self.Adend3NMDA  = SynapseNMDA(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
        ### Synapses onto Adend1
        # self.syn_list = ['Adend3AMPAf','Adend3NMDA','Adend3GABAf', 'Adend2GABAf','Adend1AMPAf','Adend1NMDA',
        #                 'somaAMPAf','somaGABAf','somaGABAfb','BdendAMPAf','BdendNMDA']

###############################################################################
# class Cck_cell(Cell):  # the previous one !!
#     name = "CCK-expressing cell"
#     def set_label(self):
#         self.label = "cck"

#     def set_morphology(self):
#         # punctual neuron
#         total_area = 10000
#         self.soma.nseg = self.nseg
#         self.soma.cm = 1.1
#         diam = np.sqrt(total_area) # um
#         L  = diam/np.pi  # um

#         #diam = 10
#         #L = 20
#         h.pt3dclear(sec=self.soma)
#         h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
#         h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)

#     def set_conductances(self):
#         self.soma.Ra = 100 # 150

#         self.soma.insert("ichan2cck")
#         self.soma.insert("ccanl")    # Ca2+ pump
#         self.soma.insert("borgka")
#         self.soma.insert("nca")      # N-type Ca2+ conductance
#         self.soma.insert("lca")      # L-type Ca2+ conducante
#         self.soma.insert("gskch")    # Ca2+ dependent K (SK) conducance
#         self.soma.insert("mykca")    # Ca2+ and voltage dependent k+ (BK) conductance
#         self.soma.insert("Ih")
#         self.soma.insert("constant")
#         self.soma.cao = 2
#         self.soma.cai = 1e-5

#         self.soma.ek = -85
#         self.soma.ena = 55
#         self.soma.eca = 130
#         self.soma.eh = -40

#         self.soma(0.5).ichan2cck.gnabar = 0.188*0.9*1.2
#         self.soma(0.5).ichan2cck.gkbar  = 0.013*0.85*1.2
#         self.soma(0.5).ichan2cck.gl     = 0.00006*0.6*1.027
#         self.soma(0.5).ichan2cck.el     = -70.0
#         self.soma(0.5).borgka.gkabar    = 0.0006
#         #self.soma(0.5).borgka.ek       = -85.0
#         self.soma(0.5).nca.gncabar      = 0.0000016 # check to modify- original 0.004
#         self.soma(0.5).lca.glcabar      = 0.000025
#         self.soma(0.5).gskch.gskbar     = 0.0004
#         self.soma(0.5).mykca.gkbar      = 0.072
#         self.soma(0.5).Ih.gkhbar        = 0.000025
#         self.soma(0.5).Ih.alpha         = 100
#         self.soma(0.5).Ih.slope         = 10
#         self.soma(0.5).Ih.amp           = 0.01
#         self.soma(0.5).Ih.taumin        = 20
#         self.soma(0.5).constant.ic      = 0#-0.00455

#     def set_synapses(self):
#         # external noise 
#         self.somaAMPA_noise    = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)   
#         self.somaGABA_noise    = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
#         # external inputs 
#         self.somaAMPA_ec3180   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)                                         # not initially used
#         self.somaNMDA_ec3180   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)  # not initially used
#         self.somaAMPA_ec3360   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
#         self.somaNMDA_ec3360   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
#         self.somaAMPA_pyrCA3   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
#         self.somaNMDA_pyrCA3   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
#         self.somaGABA_sep180   = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
#         self.somaGABA_sep360   = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
#         # connections
#         self.somaAMPA_pyr      = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
#         self.somaNMDA_pyr      = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
#         self.somaGABA_cck      = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 

class Cck_cell(Cell): # the correct one!!
    name = "CCK-expressing cell"
    def set_label(self):
        self.label = "cck"

    def set_morphology(self):
        # punctual neuron
        total_area = 10000
        self.soma.nseg = self.nseg
        self.soma.cm = 1.1
        diam = np.sqrt(total_area) # um
        L  = diam/np.pi  # um

        #diam = 10
        #L = 20
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)

    def set_conductances(self):
        self.soma.Ra = 100#150

        self.soma.insert("ichan2cck")
        self.soma.insert("ccanl")    # Ca2+ pump
        self.soma.insert("borgka")
        self.soma.insert("nca")      # N-type Ca2+ conductance
        self.soma.insert("lca")      # L-type Ca2+ conducante
        self.soma.insert("gskch")    # Ca2+ dependent K (SK) conducance
        self.soma.insert("mykca")    # Ca2+ and voltage dependent k+ (BK) conductance
        self.soma.insert("Ih")
        self.soma.insert("constant")
        self.soma.cao = 2
        self.soma.cai = 1e-5

        self.soma.ek = -85
        self.soma.ena = 55
        self.soma.eca = 130
        self.soma.eh = -40
        
        # new parameters in the new reduction
        # if you want to reproduce previous results, please use the previous values 
        self.soma(0.5).ichan2cck.gnabar = 5.82722712*1e-2 #0.188*0.9*1.2 
        self.soma(0.5).ichan2cck.gkbar  = 4.84406334*1e-02 #0.013*0.85*1.2
        self.soma(0.5).ichan2cck.gl     = 0.00006*0.6*1.027
        self.soma(0.5).ichan2cck.el     = -70.0
        self.soma(0.5).borgka.gkabar    = 0.0006
        #self.soma(0.5).borgka.ek       = -85.0
        self.soma(0.5).nca.gncabar      = 0.0000016 # check to modify- original 0.004
        self.soma(0.5).lca.glcabar      = 0.000025
        self.soma(0.5).gskch.gskbar     = 0.0004
        self.soma(0.5).mykca.gkbar      = 0.072
        self.soma(0.5).Ih.gkhbar        = 0.000025
        self.soma(0.5).Ih.alpha         = 100
        self.soma(0.5).Ih.slope         = 10
        self.soma(0.5).Ih.amp           = 0.01
        self.soma(0.5).Ih.taumin        = 20
        self.soma(0.5).constant.ic      = 0#-0.00455

    def set_synapses(self):
        # external noise 
        self.somaAMPA_noise    = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)   
        self.somaGABA_noise    = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
        # external inputs 
        self.somaAMPA_ec3180   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)                                         # not initially used
        self.somaNMDA_ec3180   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)  # not initially used
        self.somaAMPA_ec3360   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.somaNMDA_ec3360   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaAMPA_pyrCA3   = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3,  e=0)
        self.somaNMDA_pyrCA3   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaGABA_sep180   = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
        self.somaGABA_sep360   = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
        # connections
        self.somaAMPA_pyr      = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
        self.somaNMDA_pyr      = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0)
        self.somaGABA_cck      = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1,  e=-80) 
          

class Cck_spyros_modified(Cell):
    name = "CCK-expressing cell"
    def set_label(self):
        self.label = "cck"

    def set_morphology(self):
        # punctual neuron
        total_area = 10000
        self.soma.nseg = self.nseg
        self.soma.cm = 1.1
        diam = np.sqrt(total_area) # um
        L  = diam/np.pi  # um

        diam = 10
        L = 20
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
        h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)

    def set_conductances(self):
        self.soma.Ra = 100#150

        self.soma.insert("ichan2cck")
        self.soma.insert("ccanl")    # Ca2+ pump
        self.soma.insert("borgka")
        self.soma.insert("nca")      # N-type Ca2+ conductance
        self.soma.insert("lca")      # L-type Ca2+ conducante
        self.soma.insert("gskch")    # Ca2+ dependent K (SK) conducance
        self.soma.insert("mykca")    # Ca2+ and voltage dependent k+ (BK) conductance
        self.soma.insert("Ih")
        self.soma.insert("constant")
        self.soma.cao = 2
        self.soma.cai = 1e-5

        self.soma.ek = -85
        self.soma.ena = 55
        self.soma.eca = 130
        self.soma.eh = -40

        self.soma(0.5).ichan2cck.gnabar = 0.188*0.9*1.2
        self.soma(0.5).ichan2cck.gkbar  = 0.013*0.85*1.2
        self.soma(0.5).ichan2cck.gl    = 0.00006*0.6
        self.soma(0.5).ichan2cck.el    = -70.0
        self.soma(0.5).borgka.gkabar    = 0.0006
        #self.soma(0.5).borgka.ek       = -85.0
        self.soma(0.5).nca.gncabar     = 0.0000016 # check to modify- original 0.004
        self.soma(0.5).lca.glcabar     = 0.000025
        self.soma(0.5).gskch.gskbar    = 0.0004
        self.soma(0.5).mykca.gkbar     = 0.072
        self.soma(0.5).Ih.gkhbar       = 0.000025
        self.soma(0.5).Ih.alpha        = 100
        self.soma(0.5).Ih.slope        = 10
        self.soma(0.5).Ih.amp          = 0.01
        self.soma(0.5).Ih.taumin       = 20
        self.soma(0.5).constant.ic     = -0.00455

    def set_synapses(self):
        self.somaGABAss = Synapse(sect=self.soma, loc=0.5, tau1=20.0, tau2=40.0, e=-80.0) # only for septal input
        self.somaAMPAf  = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)    # background
        self.somaGABAf  = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80.0)  # backgorund
        self.somaNMDA   = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)


#####################################################################################################
''' whole cell to reduce '''
####################################################################################################
# class Cck_morphology(Cell):
#     name = "CCK-expressing cell"
#     def set_label(self):
#         self.label = "cck"

#     def set_morphology(self):

#         self.add_comp("radProx1", True)
#         self.add_comp("radMed1", True)
#         self.add_comp("radDist1",True)
#         self.add_comp("lmM1", True)
#         self.add_comp("lmt1", True)
#         self.add_comp("radProx2", True)
#         self.add_comp("radMed2", True)
#         self.add_comp("radDist2", True)
#         self.add_comp("lmM2", True)
#         self.add_comp("lmt2",True)
#         self.add_comp("oriProx1", True)
#         self.add_comp("oriMed1", True)
#         self.add_comp("oriDist1",True)
#         self.add_comp("oriProx2", True)
#         self.add_comp("oriMed2",True)
#         self.add_comp("oriDist2",True)

#         """
#         h.pt3dclear(sec=self.soma)
#         h.pt3dadd(self.x,    self.y, self.z, 10, sec=self.soma)
#         h.pt3dadd(self.x+15, self.y, self.z, 10, sec=self.soma)
#         h.pt3dclear(sec=self.radProx2 )
#         h.pt3dadd(self.x+15, self.y,    self.z, 4, sec=self.radProx2)
#         h.pt3dadd(self.x+15, self.y+30, self.z, 4, sec=self.radProx2)
#         h.pt3dclear(sec=self.radMed2)
#         h.pt3dadd(self.x+45, self.y+30, self.z, 3, sec=self.radMed2)
#         h.pt3dadd(self.x+75, self.y+60, self.z, 3, sec=self.radMed2)
#         h.pt3dclear(sec=self.radDist2 )
#         h.pt3dadd(self.x+75, self.y+60, self.z, 2, sec=self.radDist2)
#         h.pt3dadd(self.x+90, self.y+75, self.z, 2, sec=self.radDist2)
#         h.pt3dclear(sec=self.lmM2)
#         h.pt3dadd(self.x+90,  self.y+75, self.z, 1.5, sec=self.lmM2)
#         h.pt3dadd(self.x+105, self.y+90, self.z, 1.5, sec=self.lmM2)
#         h.pt3dclear(sec=self.lmt2)
#         h.pt3dadd(self.x+105, self.y+90,  self.z, 1, sec=self.lmt2)
#         h.pt3dadd(self.x+120, self.y+105, self.z, 1, sec=self.lmt2)
#         h.pt3dclear(sec=self.radProx1)
#         h.pt3dadd(self.x,    self.y,    self.z, 4, sec=self.radProx1)
#         h.pt3dadd(self.x-14, self.y+15, self.z, 4, sec=self.radProx1)
#         h.pt3dclear(sec=self.radMed1)
#         h.pt3dadd(self.x-14, self.y+15, self.z, 3, sec=self.radMed1)
#         h.pt3dadd(self.x-29, self.y+30, self.z, 3, sec=self.radMed1)
#         h.pt3dclear(sec=self.)
#         h.pt3dadd(self.x-29, self.y+30, self.z, 2, sec=self.radDist1)
#         h.pt3dadd(self.x-44, self.y+45, self.z, 2, sec=self.radDist1)
#         h.pt3dclear(sec=self.lmM1)
#         h.pt3dadd(self.x-44, self.y+45, self.z, 1.5, sec=self.lmM1)
#         h.pt3dadd(self.x-59, self.y+60, self.z, 1.5, sec=self.lmM1)
#         h.pt3dclear(sec=self.lmt1)
#         h.pt3dadd(self.x-59, self.y+60, self.z, 1, sec=self.lmt1)
#         h.pt3dadd(self.x-89, self.y+90, self.z, 1, sec=self.lmt1)
#         h.pt3dclear(sec=self.oriProx1)
#         h.pt3dadd(self.x-29, self.y, self.z, 2, sec=self.oriProx1)
#         h.pt3dadd(self.x-29, self.y, self.z, 2, sec=self.oriProx1)
#         h.pt3dclear(sec=self.oriMed1)
#         h.pt3dadd(self.x-29, self.y-29, self.z, 1.5, sec=self.oriMed1)
#         h.pt3dadd(self.x-59, self.y-59, self.z, 1.5, sec=self.oriMed1)
#         h.pt3dclear(sec=self.oriDist1)
#         h.pt3dadd(self.x-59, self.y-59, self.z, 1.0, sec=self.oriDist1)
#         h.pt3dadd(self.x-89, self.y-89, self.z, 1.0, sec=self.oriDist1)
#         h.pt3dclear(sec=self.oriProx2)
#         h.pt3dadd(self.x+15, self.y,    self.z, 2, sec=self.oriProx2)
#         h.pt3dadd(self.x+45, self.y-29, self.z, 2, sec=self.oriProx2)
#         h.pt3dclear(sec=self.oriMed2)
#         h.pt3dadd(self.x+45, self.y-29, self.z, 1.5, sec=self.oriMed2)
#         h.pt3dadd(self.x+75, self.y-59, self.z, 1.5, sec=self.oriMed2)
#         h.pt3dclear(sec=self.oriDist2)
#         h.pt3dadd(self.x+75,  self.y-59, self.z, 1.0, sec=self.oriDist2)
#         h.pt3dadd(self.x+105, self.y-89, self.z, 1.0, sec=self.oriDist2)
#         """
#         self.soma.L = 20
#         self.soma.diam = 10
#         self.radProx1.L = 100
#         self.radProx1.diam = 4
#         self.radMed1.L = 100
#         self.radMed1.diam = 3
#         self.radDist1.L = 200
#         self.radDist1.diam = 2
#         self.lmM1.L = 100
#         self.lmM1.diam = 1.5
#         self.lmt1.L = 100
#         self.lmt1.diam = 1
#         self.radProx2.L = 100
#         self.radProx2.diam = 4
#         self.radMed2.L = 100
#         self.radMed2.diam = 3
#         self.radDist2.L = 200
#         self.radDist2.diam = 2
#         self.lmM2.L = 100
#         self.lmM2.diam = 1.5
#         self.lmt2.L = 100
#         self.lmt2.diam = 1
#         self.oriProx1.L = 100
#         self.oriProx1.diam = 2
#         self.oriMed1.L = 100
#         self.oriMed1. diam = 1.5
#         self.oriDist1.L = 100
#         self.oriDist1.diam = 1
#         self.oriProx2.L = 100
#         self.oriProx2.diam = 2
#         self.oriMed2.L = 100
#         self.oriMed2.diam = 1.5
#         self.oriDist2.L = 100
#         self.oriDist2.diam = 1

#         for sect in self.all_sec:
#             sect.Ra = 100#150
#             sect.cm = 1.1
#             sect.nseg = int((sect.L/(0.1*h.lambda_f(100)) + .9)/2)*2 + 1

#         self.radProx1.connect(self.soma(0),0)
#         self.radMed1.connect(self.radProx1(1),0)
#         self.radDist1.connect(self.radMed1(1),0)
#         self.lmM1.connect(self.radDist1(1),0)
#         self.lmt1.connect(self.lmM1(1),0)
#         self.radProx2.connect(self.soma(1),0)
#         self.radMed2.connect(self.radProx2(1),0)
#         self.radDist2.connect(self.radMed2(1),0)
#         self.lmM2.connect(self.radDist2(1),0)
#         self.lmt2.connect(self.lmM2(1),0)
#         self.oriProx1.connect(self.soma(0),0)
#         self.oriMed1.connect(self.oriProx1(1),0)
#         self.oriDist1.connect(self.oriMed1(1),0)
#         self.oriProx2.connect(self.soma(1),0)
#         self.oriMed2.connect(self.oriProx2(1),0)
#         self.oriDist2.connect(self.oriMed2(1),0)

#     def set_conductances(self):
#         self.soma.insert("constant")
#         for seg in self.soma:
#             seg.constant.ic = -0.00455

#         for sect in self.all_sec:
#             #sect.Ra = 100#150
#             #sect.cm = 1.1
#             sect.insert("ichan2cck")
#             sect.insert("ccanl")    # Ca2+ pump
#             sect.insert("borgka")
#             sect.insert("nca")      # N-type Ca2+ conductance
#             sect.insert("lca")      # L-type Ca2+ conducante
#             sect.insert("gskch")    # Ca2+ dependent K (SK) conducance
#             sect.insert("mykca")    # Ca2+ and voltage dependent k+ (BK) conductance
#             sect.insert("Ih")
#             #sect.insert("constant")
#             sect.cao = 2
#             sect.cai = 1e-5

#             sect.ek = -85
#             sect.ena = 55
#             sect.eca = 130
#             sect.eh = -40

#             for seg in sect:
#                 seg.ichan2cck.gnabar = 0.188*0.9*1.2
#                 seg.ichan2cck.gkbar  = 0.013*0.85*1.2
#                 seg.ichan2cck.gl     = 0.00006*0.6*1.027 # to make the validation at 34 celsius
#                 seg.ichan2cck.el     = -70.0
#                 seg.borgka.gkabar    = 0.0006
#                 #self.soma(0.5).borgka.ek       = -85.0
#                 seg.nca.gncabar      = 0.0000016 # check to modify- original 0.004
#                 seg.lca.glcabar      = 0.000025
#                 seg.gskch.gskbar     = 0.0004
#                 seg.mykca.gkbar      = 0.072
#                 seg.Ih.gkhbar        = 0.000025
#                 seg.Ih.alpha         = 100
#                 seg.Ih.slope         = 10
#                 seg.Ih.amp           = 0.01
#                 seg.Ih.taumin        = 20
#                 #sect(0.5).constant.ic      = 0#-0.00455

#     def set_synapses(self):
#         self.somaAMPAf = Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0.0)   # background
#         self.somaGABAf = Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80.0) # backgorund
#         self.somaNMDA  = SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15.0, tau2NMDA=150.0, r=1, e=0.0)
