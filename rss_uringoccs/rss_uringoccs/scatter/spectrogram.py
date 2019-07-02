import numpy as np
from scipy.signal import spectrogram
from ..tools.write_output_files import construct_filepath

# Function to compute the continuous STFT
#   positional arguments:
#      time (*np.ndarray*): one-dimensional array of times at which
#                       signal is sampled -- MUST BE UNIFORM
#       signal (*np.ndarray*): the uniformly-sampled signal
#   keyword arguments:
#       numpts (*int*): number of points per STFT segment
#       nsegs (*int*): number of segments in STFT
#   Notes:
#        both numpts and nsegs MUST be smaller than the length of time
def cont_stft(time,signal,numpts=int(1e3),nsegs=int(5e2)):
    # compute sample spacing
    dt = (time[-1]-time[0])/float(len(time))
    # compute frequency values for each segment assuming frequency
    # range is -1/(2*sample rate) to 1/(2*sample rate)
    f = np.linspace(-0.5/dt,0.5/dt,numpts)
    # compute width of time bin for each segment
    segwid = dt*float(numpts)/2
    #print(segwid)
    # compute time segment values
    t = np.linspace(time[0]+segwid,time[-1]-segwid,nsegs)

    # compute continuous STFT
    Sxx = []
    for ti in t:
        # set time limits for bin
        tmin = ti-segwid
        tmax = ti+segwid
        # create mask array to select data within bin
        mask = [(time>=tmin)&(time<=tmax)]
        # compute continuous FT
        Sxi = [abs(np.sum(signal[mask]*np.exp(-2j*np.pi*fi*time[mask])))**2. for fi in f]
        # store
        Sxx += [Sxi]
        #print(t[0],ti,t[-1])
        pbar((ti-t[0])/(t[-1]-t[0]))

    return t,f,np.array(Sxx).T

# stack spectrogram
#   time (np.ndarray): Jx1 array of times at which the spectrogram was computed
#   Sxx (np.ndarray): IxJ array of spectrogram power values
#   N (int): number of FFT segments to include in each bin,
#               must be less than J, default is 32
def stack_spec(time,Sxx,N=int(32)):
    t = []
    S = []
    for i in range(0,len(time)-N,N):
        Si = np.zeros(len(Sxx[:,i]))
        for j in range(N):
            Si += Sxx[:,i+j]
        S += [Si]
        t += [time[i+int(N/2)]]
    t = np.array(t)
    S = np.array(S).T
    return t,S
# compute spectrogram from signal
#   time (np.ndarray): array of times at which signal was measured
#   signal (np.ndarray): array of signal values for which to compute
#                   the spectrogram, must match length of times
#   stack (bool): boolean specifying whether to stack the spectrogram
#                   to improve signal-to-noise, default is True
#   nstack (int): number of FFT segments to include in each bin,
#                   default is 32
def spectro(time,signal,stack=True,nstack=int(16),hires=False,numpts=None,nsegs=None):
    # get sampling frequency
    df = float(len(time))/(time[-1]-time[0])
    # set frequency sampling in spectrogram to one tenth true sampling
    nperseg = int(df/10.)
    twid = float(nperseg)/df
    # hi res option -- slow!
    if hires:
        # correct hi res STFT options if not properly set
        if numpts == None or numpts > len(time) :
            numpts = nperseg
        if nsegs == None or nsegs > len(time) :
            nsegs = int(len(time)/nperseg)
        # compute continuous STFT
        time,freqs,Sxx = cont_stft(time,signal,numpts=numpts,nsegs=nsegs)
    else:
        # compute two-sided power spectrogram using a Hamming window
        f,t,S = spectrogram(signal,df,return_onesided=False,scaling='spectrum',
                        nperseg=nperseg,mode='magnitude',detrend=False,window='hamming')
        # sort and scale spectrogram values
        time = t+time[0]        # scale times to starting time value
        a = np.argsort(f)       # sort frequencies from smallest to highest
        freqs = -f[a]           # sort frequencies and introduce a sign
                                # change to match Gresh 1989
        Sxx = abs(S[a,:])**2    # sort power along frequency axis to match
                                # the frequency sorting
    # stack if desired
    if stack:
        time,Sxx = stack_spec(time,Sxx,N=nstack)
    # return time, frequency, and power
    return time,freqs,Sxx
# run spectrogram code and output results to file
#   rsr_inst (*obj*): instance of RSR reader
#   geo_inst (*obj*): instance of Geometry
#   cal_inst (*obj*): instance of Calibration
#   rho_limits (*list*): 2x1 list of radii boundaries over which to
#                       compute the spectrogram
#   stack (*bool*): specifying whether to stack the resulting spectrogram
#   nstack (*int*): number of spectrogram slices to stack in each bin
#   hires (*bool*): specifying whether to compute spectrogram "manually"
#                       in a hi-res time sampling mode
def Scatter(spm_raw,t_oet,rho_km,spm,IQ,pnorm,rev_info,rho_limits=[4.5e4,6e4],
            stack=True,nstack=int(16),hires=False,numpts=None,nsegs=None):

    # obtain relevant info from each instance
    rho = np.interp(spm,t_oet,rho_km)
    # clip to relevant region
    clp = [(rho>=rho_limits[0])&(rho<=rho_limits[1])]
    # obtain stacked spectrogram
    t,f,S = spectro(spm[clp],IQ[clp],stack=stack,
                        nstack=nstack,hires=hires,numpts=numpts,nsegs=nsegs)
    # radii corresponding to spectrogram times
    r = np.interp(t,t_oet,rho_km)
    # freespace power fit at spectrogram times
    Snorm = np.interp(t,spm,pnorm)
    #generate plot file names
    outfile,outdir = construct_filepath(rev_info,'SCATTER')
    for file,dir in zip(outfile,outdir):
        # output results to file
        out = open(dir+file+'.TAB','w')
        f0 = 'f0 (Hz): '.rjust(18)+str(f[0])
        t0 = 't0 (SPM): '.rjust(18)+str(t[0])
        nf = 'n_freq: '.rjust(18)+str(len(f))
        nt = 'n_time: '.rjust(18)+str(len(t))
        df = 'Delta t (sec): '.rjust(18)+str((f[-1]-f[0])/float(len(f)))
        dt = 'Delta f (Hz): '.rjust(18)+str((t[-1]-t[0])/float(len(t)))
        header = '\n'.join(val for val in [f0,nf,df,t0,nt,dt])
        out.write(header+'\n')
        # write spectrogram
        for i in range(len(t)):
            outrow = ','.join(str(val) for val in [r[i],*S[:,i]/Snorm[i]])+'\n'
            out.write(outrow)
        out.close()

    return None
