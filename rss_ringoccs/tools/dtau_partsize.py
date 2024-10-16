'''
Purpose:
    Compute grids of predicted values of :math:`\\Delta\\tau/\\tau`
    for Ka, X, and S bands. Predictions are based on Mie scattering
    physics calculated from user-specified refractive index :math:`\\bar{m}`
    and an assumed power law particle size distribution
    :math:`N\\propto a^{-q}` for radius :math:`a` and power law index
    :math:`q`. These predictions are dependent on the PyMieScatt
    software package, which can be obtained by pip install. See their
    `online documentation <https://pymiescatt.readthedocs.io/en/latest/#>`_
    for details.

Notes:
    #.  PyMieScatt is a dependency for these predictions. No other part
        of ``rss_ringoccs`` is dependent on PyMieScatt.
    #.  Plotting will `not` produce an output file.
    #.  Suggested values for input variables for the particle
        size distribution and composition for Saturnian rings
        are as follows
            ============     ======================================
            Variable         Value
            ============     ======================================
            q                np.arange(2.8,3.6,0.2)
            a_min_vals       np.logspace(6,9,301)
            a_max_vals       np.array([3e9,1e10])
            m                1.28+0j (fluffy water ice), 1.78+0j (solid water ice)
            ============     ======================================
'''
###
### Import libraries
###
import numpy as np
import os
from scipy.integrate import simpson as simps
import PyMieScatt as ps
###
###   ~~ SUGGESTED INPUT MODEL OPTIONS ~~
###
# power law slope values
q_vals = [2.8,3.0,3.2,3.4]
# refractive index (1.78 = solid, 1.33 = composite, 1.28 = fluffy)
m = 1.78+0j
# min particle sizes in nanometers
a_min_vals = np.logspace(6,9,301)
# max particle sizes in nanometers
a_max_vals = [3e9,1e10]

def n(a,q):
    '''
    Purpose:
        Compute number of particles of given radius  assuming an
        unnormalize power law particle size distribution with
        power law index of q.

    Arguments:
            :a (*float* or *np.ndarray*): particle size
            :q (*float*):                 power law index
    '''
    return a**-q
def tau_int(a,Q,q=3.1):
    '''
    Purpose:
        Compute the integrand of :math:`\\tau` for given particle size
        and scattering efficiency.

    Arguments:
        :a (*np.ndarray*):  particle size
        :Q (*np.ndarray*):  scattering efficiency computed by PyMieScatt
    '''

    return a**2. * Q * n(a,q)
def dtau_int(a,dQ,q=3.1):
    '''
    Purpose:
        Compute the integrand of :math:`\\Delta\\tau` for given particle
        size and differential scattering efficiency.

    Arguments:
        :a (*np.ndarray*):  particle size
        :dQ (*np.ndarray*): differential scattering efficiency computed by PyMieScatt
    '''
    return a**2. * dQ * n(a,q)
def calc_dtau_pred(m,q_vals,a_min_vals,a_max_vals,w_K=9e6,w_X=3.6e7,w_S=13.e7,write=True):
    '''
    Purpose:
        Predict a set of :math:`\\Delta\\tau/\\tau` values for given
        particle compositions and size distribution

    Arguments:
        :m (*complex*):         refractive index :math:`\\bar{m}`
        :q_vals (*list*):       particle size distribution slope
        :a_min_vals (*list*):   range of minimum particle sizes in nanometers
        :a_max_vals (*list*):   range of maximum particle sizes in nanometers
        :w_K (*float*):         lowest wavelength reference band (Ka) in nanometers
        :w_X (*float*):         middle wavelength reference band (X) in nanometers
        :w_S (*float*):         highest wavelength reference band (S) in nanometers

    Keyword Arguments:
        :w_K (*float*):         smallest wavelength, default is radio Ka at 9e6 nm (0.9 cm)
        :w_X (*float*):         median wavelength, default is radio X at 3.6e7 nm (3.6 cm)
        :w_S (*float*):         largest wavelength, default is radio S at 1.3e8 nm (13 cm)
        :write (*bool*):        boolean specifying whether to write predictions to
                                a CSV file in the current working directory from
                                which this function is called. Files for different
                                refractive indices are differentiated using the
                                index of refraction :math:`\\bar{m}` with real
                                component R.RR and imaginary component I.II such
                                that the file name is
                                ./dtau_miescatt_partsize_grid_mRRRiIII.csv
                                For example, solid ice has a refractive index of
                                :math:`\\bar{m}=1.78+i0.00`, which gives a file
                                name of
                                ./dtau_miescatt_partsize_grid_m178i000.csv

    Returns:
        :dtau_KX (*np.ndarray*): N x M array of the differential optical depth
                            for Ka and X band, where N = length of q_vals times
                            the length of a_max_vals, and M = length of a_min_vals
        :dtau_XS (*np.ndarray*): N x M array of the differential optical depth
                            for Ka and X band, where N = length of q_vals times
                            the length of a_max_vals, and M = length of a_min_vals
    '''

    # storage lists for plotting
    tau_ratio_KX = [] # (tau_K - tau_X)/tau_x
    tau_ratio_XS = [] # (tau_X - tau_S)/tau_X

    ## Get differential Mie extinction efficiencies for full range of radii
    # set PyMieScatt keywargs
    drange = (2*a_min_vals[0],2*a_max_vals[-1])
    kwargs = {'nMedium':1.0,'diameterRange':drange,'nd':1000,'logD':True}
    # use PyMieScatt to compute Mie extinction efficiencies at three wavelengths
    d,Q_K = ps.MieQ_withDiameterRange(m,w_K,**kwargs)[0:2]
    d,Q_X = ps.MieQ_withDiameterRange(m,w_X,**kwargs)[0:2]
    d,Q_S = ps.MieQ_withDiameterRange(m,w_S,**kwargs)[0:2]
    # differential extinction efficiencies from extinction efficiencies
    dQ_KX = Q_K - Q_X
    dQ_XS = Q_X - Q_S
    # particle radii from diameters
    a_vals = d/2

    # output file
    if write:
        #current working directory
        cwd = os.getcwd()
        # strings for index of refraction
        mreal = str(int(np.real(m)*100))
        if len(mreal) == 2:
            mreal = mreal + '0'
        elif len(mreal) == 1:
            mreal = mreal + '00'
        mimag = str(int(np.imag(m)*100))
        if len(mimag) == 2:
            mimag = mimag + '0'
        elif len(mreal) == 1:
            mimag = mimag + '00'
        # create file
        out = open(cwd+'/dtau_miescatt_partsize_grid_m'+mreal+'i'+mimag+'.csv','w')
        heads = ['log10(a_min (mm) )','log10(a_max (mm) )','q','Delta tau KX','Delta tau XS']
        out.write(",".join(h.ljust(12) for h in heads)+"\n")

    # compute integrals for different radii ranges
    for a_max in a_max_vals:

        # compute integrals for different power law slopes
        for qi in q_vals:

            # new empty lists for storage
            Tau_X = []
            DeltaTau_KX = []
            DeltaTau_XS = []

            # compute dtau/tau for a range of minimum particle radii
            for a_min in a_min_vals:

                # clipping mask array based on particle radius
                rclip = [(a_vals>=a_min)&(a_vals<=a_max)]
                # compute tau integral with Simpson's rule and store
                tx = simps(tau_int(a_vals[rclip],Q_X[rclip],q=qi),x=a_vals[rclip])
                Tau_X += [tx]
                # compute Delta tau integrals with Simpson's rule and store in lists
                tkx = simps(dtau_int(a_vals[rclip],dQ_KX[rclip],q=qi),x=a_vals[rclip])
                txs = simps(dtau_int(a_vals[rclip],dQ_XS[rclip],q=qi),x=a_vals[rclip])
                DeltaTau_KX += [tkx]
                DeltaTau_XS += [txs]
                if write:
                    row = [np.log10(a_min)-6,np.log10(a_max)-6,qi,tkx/tx,txs/tx]
                    out.write(",".join(str(round(val,8)) for val in row)+"\n")

            # convert lists to numpy arrays
            Tau_X = np.array(Tau_X)
            DeltaTau_KX = np.array(DeltaTau_KX)
            DeltaTau_XS = np.array(DeltaTau_XS)

            # store arrays for later
            tau_ratio_KX += [DeltaTau_KX/Tau_X]
            tau_ratio_XS += [DeltaTau_XS/Tau_X]

    if write:
        out.close()

    return np.array(tau_ratio_KX),np.array(tau_ratio_XS)


def plot_dtau_grid(m,q_vals,a_min_vals,a_max_vals):
    '''
    Purpose:
        Plot a set of :math:`\\Delta\\tau/\\tau` predictions for given
        particle compositions and size distribution and save
        plot to current working directory.

    Arguments:
        :m (*complex*):         refractive index
        :q_vals (*list*):       particle size distribution slope
        :a_min_vals (*list*):   range of minimum particle sizes in nanometers
        :a_max_vals (*list*):   range of maximum particle sizes in nanometers
        :w_K (*float*):         lowest wavelength reference band (Ka) in nanometers
        :w_X (*float*):         middle wavelength reference band (X) in nanometers
        :w_S (*float*):         highest wavelength reference band (S) in nanometers

    Returns:
        None
    '''
    ###
    ###     ~~ COMPUTE MODEL GRIDS ~~
    ###
    dtau_KX,dtau_XS = calc_dtau_pred(m,q_vals,a_min_vals,a_max_vals,write=False)
    ###
    ###    ~~ PLOT RESULTS! ~~
    ###
    import matplotlib.pyplot as plt
    # plot zero-points for reference
    plt.plot([-100,100],[0,0],dashes=[12,4],color='0.5')
    plt.plot([0,0],[-100,100],dashes=[12,4],color='0.5')
    # set up color and labels
    colors = []
    labels = []
    for i,a_max in enumerate(a_max_vals):
        for j,q in enumerate(q_vals):
            labels += [r'$a_{max}=$'+str(a_max/1e9)+' m']
            colors += [str(float(i)/float(len(a_max_vals)))]
    # iterate over all model power laws and max particle sizes
    for i in range(len(q_vals)*len(a_max_vals)):
        # grid line
        if labels[i] != labels[i-1]:
            pkwargs = {'color':colors[i],'label':labels[i],'lw':1}
        else:
            pkwargs = {'color':colors[i],'lw':1}
        plt.plot(100.*dtau_KX[i],100.*dtau_XS[i],**pkwargs)
        # a_min markers
        pkwargs = {'color':'none','ls':'none','marker':'o','markeredgecolor':colors[i]}
        for j in range(0,len(dtau_KX[i]),50):
            plt.plot(100.*dtau_KX[i][j],100.*dtau_XS[i][j],**pkwargs)
        # q label
        if i < len(q_vals):
            tx = 100.*dtau_KX[i][0]+5
            ty = 100.*dtau_XS[i][0]-2
            plt.text(tx,ty,'q='+str(q_vals[i]))
    # labels and output
    plt.xlim(-30,100)
    plt.ylim(-30,50)
    plt.xlabel(r'$\Delta\tau_{Ka-X} / \tau_X$ as %')
    plt.ylabel(r'$\Delta\tau_{X-S} / \tau_X$ as %')
    plt.legend(loc=4,frameon=False,numpoints=1)
    plt.title(r'Scattering Model with $\overline{m}=$'+str(m))
    cwd = os.getcwd()
    plt.savefig(cwd+'/dtau_dtau_m'+str(np.real(m))+'i'+str(np.imag(m))+'.png',dpi=128)
    plt.close()
