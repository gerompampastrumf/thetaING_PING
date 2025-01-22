import numpy as np
import pandas as pd
import scipy 
from mne import create_info, EpochsArray
import mne
from pactools import Comodulogram

# import file_management
from spectrum import *
from scipy.signal import hilbert
import pandas as pd
from skimage.measure import find_contours

def firingRate_fromSpikes(tvec,tmax=None,dt=0.1,sg=5,mode="gaussian"):  
    tvec = np.array(tvec)
    if(tmax is None):
        tmax = np.max(tvec)
    time_binned = np.arange(0, tmax, dt)
    rates,t = np.histogram(tvec,bins=time_binned)
    if(sg!=0):
        if(mode=="gaussian"):
            nlen = 10
            len_win = nlen*sg if nlen*sg%2==1 else nlen*sg+1
            pos    = np.arange(0,len_win,1)-len_win//2
            window = np.exp(-pos**2/(2*sg**2))/(sg*np.sqrt(2*np.pi))
        #elif(mode=="numpy"):
        #    return time_binned,scipy.ndimage.gaussian_filter1d(rates,sg)
        elif(mode=="exp"):
            nlen = 10
            len_win = nlen*sg if nlen*sg%2==1 else nlen*sg+1
            pos    = np.arange(0,len_win,1)-len_win//2
            window = np.exp(-pos/(sg))/sg
            window[pos<0]=0
        else:
            window = np.ones(sg)/sg
        rates  = np.convolve(rates,window,mode='same')
    return time_binned, rates


def continufySpikes2(spikes,fs,label,w=0,tmax=None,mode="mean",dt=0.1):
    timeseries={label:[],"iseed":[]}
    if(tmax is None):
        tmax=np.max(spikes["tvec"].values)
    for iseed in np.unique(spikes["iseed"].values):
        _,rates = firingRate_fromSpikes(spikes[spikes["iseed"]==iseed]["tvec"].values,tmax=tmax,mode=mode,sg=w,dt=dt)
        if(dt/1000<1/fs):
            rescale = int(1/((dt/1000)*fs))
            rates = np.sum(rates.reshape((int(len(rates)/rescale), int(rescale))),axis=1)
        timeseries[label].extend(rates)
        timeseries["iseed"].extend(np.zeros(len(rates),np.int32)+iseed)
    return pd.DataFrame(timeseries)


def continufySpikes(spikes,fs,label,w=0,tmax=None):
    timeseries={label:[],"iseed":[]}
    time_bins = np.arange(0,np.max(spikes["tvec"].values),1000/fs)
    for iseed in np.unique(spikes["iseed"].values):
        rates,t=np.histogram(spikes[spikes["iseed"]==iseed]["tvec"].values,bins=time_bins)
        t=t[:-1]
        if(tmax is not None):
            rates = rates[t<tmax]
            t = t[t<tmax]
        if(w):
            rates=np.convolve(rates, np.ones(w)/w, mode='same')
        timeseries[label].extend(rates)
        timeseries["iseed"].extend(np.zeros(len(rates))+iseed)
    return pd.DataFrame(timeseries)

def bound_timenphases(time,phase):
    ts_max = time[:-1][(phase[:-1]>3*np.pi/4) & (phase[1:]<np.pi/4)]
    ts_min = time[1:][(phase[:-1]>3*np.pi/4) & (phase[1:]<np.pi/4)]
    ts_twopi = ts_max+np.random.uniform(0.0+1e-6,(ts_min-ts_max)-1e-6)
    full_times = np.concatenate((time,ts_twopi,ts_twopi+1e-8))
    full_phase = np.concatenate((phase,np.ones(len(ts_twopi))*2*np.pi,np.zeros(len(ts_twopi))))
    arr1inds = full_times.argsort()
    return full_times[arr1inds],full_phase[arr1inds]

def get_activity(pyr_spikes,t0=1250,measure_type="firingRate"):
    activity=[]
    if(len(pyr_spikes["tvec"].values)==0):
        return 0.0,0.0
    Tmax = np.max(pyr_spikes["tvec"].values)
    N_neurons = len(np.unique(pyr_spikes["idvec"].values))
    for iseed in np.unique(pyr_spikes["iseed"].values):
        times = pyr_spikes[(pyr_spikes["iseed"]==iseed)]["tvec"].values
        
        times = times[times>t0]
        if(measure_type=="firingRate"):
            act_temp = len(times)/(Tmax-t0)*1000/N_neurons
        elif(measure_type=="perCycle"):
            act_temp = len(times)/(Tmax-t0)*125/N_neurons
        activity.append(act_temp)
    return np.mean(activity),np.std(activity)

def plot_contours(area,ftheta,fgamma,ax,thr=0.90):
    area = np.copy(area)
    thres = np.quantile(area,thr)
    area[area<thres]=0
    all_conts=find_contours(area)
    conts = all_conts[np.argmax([len(conts) for conts in all_conts])]
    ax.plot(ftheta[np.round(conts[:,1]).astype(int)],fgamma[np.round(conts[:,0]).astype(int)],c="k")

def hilbert_phase(x):
    sig=x-x.mean()
    std = sig.std()
    #print("Len std: ",len(sig),std)
    sig/=std
    analytic_sig = hilbert(sig)
    instantaneous_phase = np.angle(analytic_sig)
    amplitude_envelope  = np.abs(analytic_sig)
    return np.mod(instantaneous_phase,2*np.pi), (amplitude_envelope+x.mean())*std

def is_valid_move(matrix, visited, row, col):
    rows, cols = len(matrix), len(matrix[0])
    return (0 <= row < rows and 0 <= col < cols and
            matrix[row][col] == 1 and not visited[row][col])

def find_clusters(matrix):
    rows, cols = len(matrix), len(matrix[0])
    visited = [[False for _ in range(cols)] for _ in range(rows)]
    clusters = []

    for i in range(rows):
        for j in range(cols):
            if matrix[i][j] == 1 and not visited[i][j]:
                stack = [(i, j)]
                cluster = []

                while stack:
                    current_row, current_col = stack.pop()
                    if not visited[current_row][current_col]:
                        visited[current_row][current_col] = True
                        cluster.append((current_row, current_col))

                        # Explore neighbors
                        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]
                        for move in moves:
                            new_row, new_col = current_row + move[0], current_col + move[1]
                            if is_valid_move(matrix, visited, new_row, new_col):
                                stack.append((new_row, new_col))

                clusters.append(cluster)
    return clusters

import warnings
def compute_mi_conelty(timeseries, ftheta, fgamma, fs=1e3, surrogate_test=False, nsurrogates=100,Filter_theta=2.0,Filter_gamma=10,method="canolty"):
    estimator = Comodulogram(fs=fs, low_fq_range=ftheta,high_fq_range=fgamma,low_fq_width=Filter_theta,high_fq_width=Filter_gamma,method=method,progress_bar=False)
    estimator.fit(timeseries)
    mi = estimator.comod_
    if(surrogate_test):
        warnings.warn("CFC surrogate not yet made",UserWarning)
    return ftheta, fgamma, mi, np.zeros((nsurrogates, len(ftheta), len(fgamma)))

def compute_mi_conelty_diff_high_low(timeseries1,timeseries2, ftheta, fgamma, fs=1e3, surrogate_test=False, nsurrogates=100,Filter_theta=2.0):
    estimator = Comodulogram(fs=fs, low_fq_range=ftheta,high_fq_range=fgamma,low_fq_width=Filter_theta, method="canolty",progress_bar=False)
    estimator.fit(timeseries1,high_sig=timeseries2)
    mi = estimator.comod_
    if(surrogate_test):
        warnings.warn("CFC surrogate not yet made",UserWarning)
    return ftheta, fgamma, mi, np.zeros((nsurrogates, len(ftheta), len(fgamma)))

from scipy import signal

def round_up_to_odd(f):
    f = int(np.ceil(f))
    return f + 1 if f % 2 == 0 else f

def eegfilt(data,srate,locutoff,hicutoff,epochframes=0,filtorder=0,revfilt=0,firtype="firls",causal=0):
    
    frames = len(data)
    nyq            = srate*0.5;  # Nyquist frequency
    MINFREQ = 0;

    minfac         = 3
    min_filtorder  = 15
    trans          = 0.15

    if locutoff > 0 and hicutoff > 0 and locutoff > hicutoff:
        raise ValueError('locutoff > hicutoff ???')

    if locutoff < 0 or hicutoff < 0:
        raise ValueError('locutoff | hicutoff < 0 ???')

    if locutoff > nyq:
        raise ValueError('Low cutoff frequency cannot be > srate/2')

    if hicutoff > nyq:
        raise ValueError('High cutoff frequency cannot be > srate/2')

    if filtorder == 0:
        if locutoff > 0:
            filtorder = minfac * int(srate / locutoff)
            filtorder=round_up_to_odd(filtorder)
        elif hicutoff > 0:
            filtorder = minfac * int(srate / hicutoff)
            filtorder=round_up_to_odd(filtorder)
        min_filtorder = 15  # Define min_filtorder if it's not defined previously
        if filtorder < min_filtorder:
            filtorder = min_filtorder
            filtorder=round_up_to_odd(filtorder)

    if epochframes == 0:
        epochframes = frames  # default

    epochs = frames // epochframes

    if epochs * epochframes != frames:
        raise ValueError('epochframes does not divide frames.')

    if filtorder * 3 > epochframes:  # MATLAB filtfilt() restriction
        print(f'eegfilt(): filter order is {filtorder}.')
        raise ValueError('epochframes must be at least 3 times the filtorder.')

    if (1 + trans) * hicutoff / nyq > 1:
        raise ValueError('high cutoff frequency too close to Nyquist frequency')

    if locutoff > 0 and hicutoff > 0:  # bandpass filter
        if revfilt:
            # print("eegfilt() - performing {}-point notch filtering.".format(filtorder))
            pass
        else:
            # print("eegfilt() - performing {}-point bandpass filtering.".format(filtorder))
            pass
        # print("If a message, 'Matrix is close to singular or badly scaled,' appears,")
        # print("then MATLAB has failed to design a good filter. As a workaround,")
        # print("for band-pass filtering, first highpass the data, then lowpass it.")

        if firtype == 'firls':
            f = [MINFREQ, (1 - trans) * locutoff / nyq, locutoff / nyq, hicutoff / nyq, (1 + trans) * hicutoff / nyq, 1]
            # print("eegfilt() - low transition band width is {:.1f} Hz; high trans. band width, {:.1f} Hz.".format((f[2] - f[1]) * srate / 2, (f[4] - f[3]) * srate / 2))
            m = [0, 0, 1, 1, 0, 0]
        elif firtype == 'fir1':
            from scipy.signal import fir1
            filtwts = fir1(filtorder, [locutoff, hicutoff] / (srate / 2))
    elif locutoff > 0:  # highpass filter
        if locutoff / nyq < MINFREQ:
            raise ValueError("eegfilt() - highpass cutoff freq must be > {} Hz".format(MINFREQ * nyq))
        # print("eegfilt() - performing {}-point highpass filtering.".format(filtorder))
        if firtype == 'firls':
            f = [MINFREQ, locutoff * (1 - trans) / nyq, locutoff / nyq, 1]
            # print("eegfilt() - highpass transition band width is {:.1f} Hz.".format((f[2] - f[1]) * srate / 2))
            m = [0, 0, 1, 1]
        elif firtype == 'fir1':
            from scipy.signal import firwin
            filtwts = firwin(filtorder + 1, locutoff / (srate / 2), pass_zero=False)
    elif hicutoff > 0:  # lowpass filter
        if hicutoff / nyq < MINFREQ:
            raise ValueError("eegfilt() - lowpass cutoff freq must be > {} Hz".format(MINFREQ * nyq))
        # print("eegfilt() - performing {}-point lowpass filtering.".format(filtorder))
        if firtype == 'firls':
            f = [MINFREQ, hicutoff / nyq, hicutoff * (1 + trans) / nyq, 1]
            # print("eegfilt() - lowpass transition band width is {:.1f} Hz.".format((f[2] - f[1]) * srate / 2))
            m = [1, 1, 0, 0]
        elif firtype == 'fir1':
            from scipy.signal import firwin
            filtwts = firwin(filtorder + 1, hicutoff / (srate / 2))
    else:
        raise ValueError('You must provide a non-0 low or high cut-off frequency')


    if revfilt:
        if firtype == 'fir1':
            raise ValueError("Cannot reverse filter using 'fir1' option")
        else:
            m = [not x for x in m]

    if firtype == 'firls':
        from scipy.signal import firls
        filtwts = firls(filtorder, f, m)

    smoothdata = np.zeros(frames)

    for e in range(epochs):  # Filter each epoch, channel
        if causal:
            smoothdata[ (e) * epochframes:(e+1) * epochframes] = signal.lfilter(filtwts, 1, data[(e ) * epochframes:(e+1) * epochframes])
        else:
            smoothdata[ (e) * epochframes:(e+1) * epochframes] = signal.filtfilt(filtwts, 1, data[(e) * epochframes:(e+1) * epochframes])
    
    return smoothdata 

def bandpass_filter(xs, norder, f_range, fs=1e4): 
    sos = scipy.signal.butter(N=norder, Wn=f_range, btype="bandpass", fs=fs, output="sos")
    return scipy.signal.sosfiltfilt(sos, xs)

def bandpass_filter_and_hilbert_transform(sig, fs, f0, df, norder):
    scale = sig.mean() 
    f_range = (f0-df, f0+df)
    sig_band = bandpass_filter(sig, norder, f_range, fs)
    phase, envelope= hilbert_phase(sig_band)
    return sig_band+scale, envelope+scale,  phase

def bandpass_filter_and_hilbert_transform2(sig, fs, f0, df, norder):
    scale = sig.mean() 
    f_range = (f0-df, f0+df)
    sig_band = bandpass_filter(sig, norder, f_range, fs)
    phase, envelope= hilbert_phase(sig_band)
    return sig_band+scale, envelope,  phase

def bandpass_filter_and_hilbert_transform3(sig, fs, f0, df, norder):
    scale = sig.mean() 
    f_range = (f0-df, f0+df)

    sig_band = eegfilt(sig-scale, fs, f_range[0],f_range[-1])
    phase, envelope= hilbert_phase(sig_band)
    
    return sig_band+scale, envelope,  phase

def get_spectra(ts, fs, fmin, fmax,choise,znorm=False):
    if(znorm):
        signal = (ts-ts.mean())/ts.std()
    else:
        signal = ts-ts.mean()

    if(choise=="FFT"):
        power = np.abs(np.fft.fft(signal))**2/len(signal)**2
        fr = np.fft.fftfreq(len(signal), 1/fs)
        inds=(fr<fmax)&(fr>fmin)
        return fr[inds],power[inds]
    elif(choise=="MULTI"):
        ch_names = ['Signal']
        ch_types = ['eeg']
        info = create_info(ch_names=ch_names, sfreq=fs, ch_types=ch_types)
        # Create an MNE EpochsArray object
        epochs_data = np.expand_dims(signal, axis=0)  # Add a dimension to match the shape
        epochs_data = np.expand_dims(epochs_data, axis=0)# Add a dimension to match the shape
        epochs = EpochsArray(data=epochs_data, info=info,verbose=False)
        power_spectrum,frequencies = mne.time_frequency.psd_multitaper(epochs, fmin=fmin, fmax=fmax,verbose=False)
        power_spectrum= power_spectrum[0,0,:]
        return frequencies,power_spectrum
    elif(choise=="WELCH"):
        nfft = 8*2**np.ceil(np.log2(fs))
        fr,power = np.abs(scipy.signal.csd(signal, signal, fs=fs, window='hann', nperseg=nfft, noverlap=nfft//2, nfft=nfft))
        inds=(fr<fmax)&(fr>fmin)
        return fr[inds],power[inds]

def get_nfft(m):
    nfft = 2
    k = 0 
    while m > nfft:
        k+=1 
        nfft = 2**k
    return nfft

def surrogate_generation(timeseries, nsurrogates, method ="block boostrapping" , replace = False, nsplits = 99): 
    surrogate_list = [] 
    n = len(timeseries)
    
    if method == "2block boostrapping":
        for _ in range(nsurrogates):
            ir = np.random.randint(low=0,high=len(timeseries))
            surrogate = np.array( list(timeseries[ir:])+list(timeseries[:ir]) )
            surrogate_list.append( surrogate )
            
    if method == "sampling":
        for _ in range(nsurrogates):
            surrogate = np.random.choice(timeseries, size=n, replace=replace)
            surrogate_list.append( surrogate )
    
    if method == "fourier":
        for _ in range(nsurrogates):
            fourier = np.fft.fft(timeseries)
            random_phases = 1j*np.random.choice(fourier.imag[:n//2], size=n//2, replace=replace)
            random_phases = np.concatenate(( random_phases, random_phases[::-1] ))
            fourier_new = fourier*random_phases
            surrogate = np.real( np.fft.ifft(fourier_new))
            surrogate_list.append( surrogate )
            
    if method == "block boostrapping":
        timeseries_splitted = np.array_split(timeseries, nsplits) 
        for _ in range(nsurrogates):
            indexes_surrogate = np.random.choice( np.arange(nsplits), size = nsplits, replace=replace)
            surrogate = []
            for index in indexes_surrogate: 
                surrogate.append(  timeseries_splitted[index] )
            surrogate_list.append( np.concatenate(surrogate) )
            
    return surrogate_list 


#################################################################################################
def compute_coherence_measurements_diff_high_low(x1,x2, fs, nperseg, noverlap, ftheta_min=4, ftheta_max=20, fgamma_min=25, fgamma_max=100, dfgamma=2,Filter_gamma=10,
                                   norder=4, nfft=2048,surrogate_test=True, nsurrogates=100,choiseSurr = 0,ibeta=4,filter_func=1):
    if(filter_func==1):
        filter_function = bandpass_filter_and_hilbert_transform3
    elif(filter_func==2):
        filter_function = bandpass_filter_and_hilbert_transform2
    psi, nu = [],[] 
    psi_surrogates = [ [] for i in range(nsurrogates) ]
    nu_surrogates  = [ [] for i in range(nsurrogates) ]
    
    fgamma = np.arange(fgamma_min,fgamma_max,dfgamma)
    #scipy.signal.csd: cross power spectral density 
    fr, pxx  = scipy.signal.csd( x1, x1, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
    
    np.random.seed(12)
    for fg in fgamma:
        _, y, _ = filter_function(x2, norder=norder, fs=fs, f0=fg, df=Filter_gamma)
        fr, pyy  = scipy.signal.csd( y, y, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
        fr, pxy  = scipy.signal.csd( y, x1, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap) 
        complex_coherence = pxy / np.sqrt( np.abs(pxx)*np.abs(pyy) ) 
        
        # compute of psi #################################################
        imin = int( (ftheta_min-fr[0])/np.diff(fr)[0] )
        imax = int( (ftheta_max-fr[0])/np.diff(fr)[0] )
        cfd = []
        # module = np.abs(complex_coherence)
        # phase  = np.angle(complex_coherence)
        for i in range(imin, imax):  
            # factor1 = module[i-ibeta//2:i+ibeta//2]
            # factor2 = module[i+1-ibeta//2:i+1+ibeta//2]
            # factor3 = np.sin( phase[i+1-ibeta//2:i+1+ibeta//2]-phase[i-ibeta//2:i+ibeta//2])
            # cfd.append( np.sum(factor1*factor2*factor3) ) 
            factor1 = np.conj( complex_coherence[ i-ibeta//2: i+1+ibeta//2 ] )
            factor2 = complex_coherence[ i+1-ibeta//2 : i+2+ibeta//2 ]
            cfd.append(  np.imag( np.sum(factor1*factor2) ) )
            
        ##################################################################
        psi.append( np.array( cfd ) )
        nu.append( np.abs( complex_coherence[imin:imax] ) ) 

        if surrogate_test:
            if(choiseSurr==0):
                surrogates = surrogate_generation( x1, nsurrogates = nsurrogates)
            elif(choiseSurr==1):
                surrogates = surrogate_generation( x1, nsurrogates = nsurrogates,nsplits=int(len(x)/nperseg))
            elif(choiseSurr==2):
                surrogates = surrogate_generation( x1, nsurrogates = nsurrogates,method="2block boostrapping")
            for l, surrogate in enumerate(surrogates): 
                yf, y, _ = filter_function(x2, norder=norder, fs=fs, f0=fg, df=10)
                
                fr, pyy  = scipy.signal.csd( y, y, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
                fr, pxy  = scipy.signal.csd( y, surrogate, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap) 

                complex_coherence = pxy / np.sqrt( pxx*pyy ) 
                
                # compute psi ######################################################################
                cfd = []
                module = np.abs(complex_coherence)
                phase  = np.angle(complex_coherence)
                
                for i in range(imin, imax):  
                    factor1 = np.conj( complex_coherence[ i-ibeta//2: i+1+ibeta//2 ] )
                    factor2 = complex_coherence[ i+1-ibeta//2 : i+2+ibeta//2 ]
                    cfd.append(  np.imag( np.sum(factor1*factor2) ) )
                    
                    # factor1 = module[i-ibeta//2:i+ibeta//2]
                    # factor2 = module[i+1-ibeta//2:i+1+ibeta//2]
                    # factor3 = np.sin( phase[i+1-ibeta//2:i+1+ibeta//2]-phase[i-ibeta//2:i+ibeta//2])
                    # cfd.append( np.sum(factor1*factor2*factor3) ) 
        
                ####################################################################################
                psi_surrogates[l].append( np.array(cfd) )
                nu_surrogates[l].append( np.abs(complex_coherence[imin:imax]) )
                
    ftheta = fr[imin:imax]
    return ftheta, fgamma, psi, nu, psi_surrogates, nu_surrogates

#################################################################################################
def compute_coherence_measurements(x, fs, nperseg, noverlap, ftheta_min=4, ftheta_max=20, fgamma_min=25, fgamma_max=100, dfgamma=2,Filter_gamma=10,
                                   norder=4, nfft=2048,surrogate_test=True, nsurrogates=100,choiseSurr = 0,ibeta=4,filter_func=1):
    if(filter_func==1):
        filter_function = bandpass_filter_and_hilbert_transform3
    elif(filter_func==2):
        filter_function = bandpass_filter_and_hilbert_transform2
    psi, nu = [],[] 
    psi_surrogates = [ [] for i in range(nsurrogates) ]
    nu_surrogates  = [ [] for i in range(nsurrogates) ]
    
    fgamma = np.arange(fgamma_min,fgamma_max,dfgamma)
    #scipy.signal.csd: cross power spectral density 
    fr, pxx  = scipy.signal.csd( x, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
    
    np.random.seed(12)
    for fg in fgamma:
        _, y, _ = filter_function(x, norder=norder, fs=fs, f0=fg, df=Filter_gamma)
        fr, pyy  = scipy.signal.csd( y, y, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
        fr, pxy  = scipy.signal.csd( y, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap) 
        complex_coherence = pxy / np.sqrt( np.abs(pxx)*np.abs(pyy) ) 
        
        # compute of psi #################################################
        imin = int( (ftheta_min-fr[0])/np.diff(fr)[0] )
        imax = int( (ftheta_max-fr[0])/np.diff(fr)[0] )
        cfd = []
        # module = np.abs(complex_coherence)
        # phase  = np.angle(complex_coherence)
        for i in range(imin, imax):  
            # factor1 = module[i-ibeta//2:i+ibeta//2]
            # factor2 = module[i+1-ibeta//2:i+1+ibeta//2]
            # factor3 = np.sin( phase[i+1-ibeta//2:i+1+ibeta//2]-phase[i-ibeta//2:i+ibeta//2])
            # cfd.append( np.sum(factor1*factor2*factor3) ) 
            factor1 = np.conj( complex_coherence[ i-ibeta//2: i+1+ibeta//2 ] )
            factor2 = complex_coherence[ i+1-ibeta//2 : i+2+ibeta//2 ]
            cfd.append(  np.imag( np.sum(factor1*factor2) ) )
            
        ##################################################################
        psi.append( np.array( cfd ) )
        nu.append( np.abs( complex_coherence[imin:imax] ) ) 

        if surrogate_test:
            if(choiseSurr==0):
                surrogates = surrogate_generation( x, nsurrogates = nsurrogates)
            elif(choiseSurr==1):
                surrogates = surrogate_generation( x, nsurrogates = nsurrogates,nsplits=int(len(x)/nperseg))
            elif(choiseSurr==2):
                surrogates = surrogate_generation( x, nsurrogates = nsurrogates,method="2block boostrapping")
            for l, surrogate in enumerate(surrogates): 
                yf, y, _ = filter_function(surrogate, norder=norder, fs=fs, f0=fg, df=10)
                
                fr, pyy  = scipy.signal.csd( y, y, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
                fr, pxy  = scipy.signal.csd( y, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap) 

                complex_coherence = pxy / np.sqrt( pxx*pyy ) 
                
                # compute psi ######################################################################
                cfd = []
                module = np.abs(complex_coherence)
                phase  = np.angle(complex_coherence)
                
                for i in range(imin, imax):  
                    factor1 = np.conj( complex_coherence[ i-ibeta//2: i+1+ibeta//2 ] )
                    factor2 = complex_coherence[ i+1-ibeta//2 : i+2+ibeta//2 ]
                    cfd.append(  np.imag( np.sum(factor1*factor2) ) )
                    
                    # factor1 = module[i-ibeta//2:i+ibeta//2]
                    # factor2 = module[i+1-ibeta//2:i+1+ibeta//2]
                    # factor3 = np.sin( phase[i+1-ibeta//2:i+1+ibeta//2]-phase[i-ibeta//2:i+ibeta//2])
                    # cfd.append( np.sum(factor1*factor2*factor3) ) 
        
                ####################################################################################
                psi_surrogates[l].append( np.array(cfd) )
                nu_surrogates[l].append( np.abs(complex_coherence[imin:imax]) )
                
    ftheta = fr[imin:imax]
    return ftheta, fgamma, psi, nu, psi_surrogates, nu_surrogates