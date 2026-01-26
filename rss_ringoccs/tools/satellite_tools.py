import rss_ringoccs_global_paths # required to allow program to run from other directories

from scipy import interpolate

global global_path_to_rss_ringoccs
global_path_to_rss_ringoccs = rss_ringoccs_global_paths.global_path_to_rss_ringoccs

global global_path_to_kernels
global_path_to_kernels = global_path_to_rss_ringoccs + 'kernels/'

# # from Stefan Renner, 2009 March09 Stefan.Renner@imcce.fr
# # Please credit Stefan when using this code.

# # Translation of posi_geo.f

# # ****************
import numpy as np
import spiceypy as spice

def moyen_mouvement(GM_p,R_p,J2,J4,J6,a,e,i):
    n= ( np.sqrt( GM_p/a**3 ) )*(	 
            1.0 +    (3.0/4.0)*J2     *((R_p/a)**2)	 
               -  (15.0/16.0)*J4     *((R_p/a)**4)	 
               +  (35.0/32.0)*J6     *((R_p/a)**6)	 
               -   (9.0/32.0)*(J2**2)*((R_p/a)**4)	 
               +  (45.0/64.0)*J2*J4 *((R_p/a)**6)	 
               + (27.0/128.0)*(J2**3)*((R_p/a)**6)	 
                   +        3.0 *J2*(e**2)*((R_p/a)**2) 	 
                   -       12.0 *J2 *(i**2)*((R_p/a)**2) )
    
    return n
	
# *******************

def precess(GM_p,R_p,J2,J4,J6,a,e,i):
	pperi=  ( np.sqrt(GM_p/a**3) )*(	 
      		    (3.0/2.0)* J2     *((R_p/a)**2)	 
      		-  (15.0/4.0)* J4     *((R_p/a)**4)	 
      		+ (105.0/16.0)* J6     *((R_p/a)**6)	 
      		- (45.0/32.0)* J2*J4 *((R_p/a)**6)	 
      		+ (27.0/64.0)*(J2**3) *((R_p/a)**6)	 
                +        3.0* J2 *(e**2)*((R_p/a)**2) 	 
                -        3.0* J2 *(i**2)*((R_p/a)**2) )	
     

	pnoeu=  ( np.sqrt(GM_p/a**3) )*(	 
      		-   (3.0/2.0)* J2     *((R_p/a)**2)	 
      		+  (15.0/4.0)* J4     *((R_p/a)**4)	 
      		- (105.0/16.0)* J6     *((R_p/a)**6)	 
      		+   (9.0/4.0)*(J2**2) *((R_p/a)**4)	 
      		-(315.0/32.0)* J2*J4 *((R_p/a)**6)	 
      		-(351.0/64.0)*(J2**3) *((R_p/a)**6)	 
                -        3.0* J2 *(e**2)*((R_p/a)**2)	 
                +   (3.0/4.0)* J2 *(i**2)*((R_p/a)**2))	

	return pperi, pnoeu

# ****************

def freq(GM_p,R_p,J2,J4,J6,a,e,i):

	eta2=  ( GM_p/a**3 )*( 1.0 -	 
      		    2.0* J2     *((R_p/a)**2) 	 
      		+  (75.0/8.0)* J4     *((R_p/a)**4) 	 
      		- (175.0/8.0)* J6     *((R_p/a)**6)  )
     
	chi2=  ( GM_p/a**3 )*( 1.0 +	 
      		    (15.0/2.0)* J2     *((R_p/a)**2) 	 
      		-  (175.0/8.0)* J4     *((R_p/a)**4) 	 
      		+ (735.0/16.0)* J6     *((R_p/a)**6)  )

	return eta2,chi2
	
# ****************

def freq_n0(GM_p,R_p,J2,J4,J6,a):
	
	n02= ( GM_p/a**3 )*( 1.0 +  
      		    1.50* J2     *((R_p/a)**2)  
      		-  (15.0/8.0)* J4     *((R_p/a)**4)  
      		+ (35.0/16.0)* J6     *((R_p/a)**6)  )

	return n02

# ****************
	
def arctan(x,y):
	# pi= np.pi
	
	# if (x !=0):
	# 	arctan= np(y/x)
	# 	if (x lt0.0) then      arctan= arctan +      pi
	# 	if (arctan lt0.0) then arctan= arctan + 2.0*pi
	# endif

	# if (x eq0.0) then begin
	# 	if (y gt0.0) then arctan=  pi/2.0
	# 	if (y lt0.0) then arctan= -pi/2.0
	# 	if (y eq0.0) then arctan=0.0 # rfrench addition
	# endif

	# return,arctan
    return np.arctan2(y,x) # note reversed arguments
	

def amodulo(x,y):
    amodulo= x % y
    if (x < 0.0):
        amodulo= amodulo + y
    return amodulo

# ****************

def geo_posi(GM,R_eq,J2,J4,J6,elemINPUT):# ,svOUTPUT # rfrench 
##	input: GM of central body (km**3/sec**2)
#	       R_eq equatorial radius of central body (km)
#	       J2,J4,J6 potential coeff
#              svOUTPUT[0:5] position (km) and velocity (km/sec)
#	NB - to minimize the possibility of index problems, 
#	       this is padded to 
#	       sv[1:6] position (km) and velocity (km/sec) 
##	input:  elemINPUT
#		a semi major axis (km) 
#               e eccentricity 
#               i INCLINATION (deg) 
#		gomega ascending node (deg) 
#		tomega periapsis longitude (deg)
#               lambda mean longitude (deg)		
#	NB - to minimize the possibility of index problems, 
#	       this is padded to elem[1:6] 

	elem	= [-999.0,elemINPUT]
	svOUTPUT	= np.zeros(7)
	
	pi=np.dpi
	ppi=2.0*pi
	rad= pi/180.0
	a = elem[1]
	e = elem[2]
	i = np.radians(elem[3])
	gomega = np.radians(elem[4])
	tomega = np.radians(elem[5])
	lambda_ = np.radians(elem[6])
	
	n = moyen_mouvement(GM,R_eq,J2,J4,J6,a,e,i)
	pperi,pnoeu = precess(GM,R_eq,J2,J4,J6,a,e,i)
	eta2,chi2= freq(GM,R_eq,J2,J4,J6,a,e,i)
	kappa= n - pperi
	nu =n - pnoeu
	
	alpha1= (1.0/3.0)*(2.0*nu + kappa)
	alpha2=2.0*nu - kappa
	alph2 = alpha1*alpha2

# NB alph2 is not alpha2!

# now use eqs. 2-7 from Renner and Sicardy
# to define intermediate variables needed to get 
# state vector

	kappa2 = kappa**2
# eq 2:
	r = a * (1.0 - e*np.cos(lambda_-tomega)   
		+ e**2 * ((3.0*eta2)/(2.0*kappa2) - 1.0 -\
			(eta2/(2.0*kappa2))*np.cos(2.0*(lambda_-tomega)))\
		+ i**2 * ((3.0*chi2)/(4.0*kappa2) - 1.0 +\
			(chi2/(4.0*alph2))*np.cos(2.0*(lambda_-gomega))) )
# eq 3:
	L = lambda_ + 2.0 * e * (n/kappa)*np.sin(lambda_-tomega)\
		+ e**2 *(3.0/4.0 + eta2/(2.0*kappa2))*\
			(n/kappa)*np.sin(2.0*(lambda_-tomega))\
		- i**2 * (chi2/(4.0*alph2)) * (n/nu) * np.sin(2.0*(lambda_-gomega))

# eq 4:
	z = a * i * ( np.sin(lambda_-gomega)\
		+ e*(chi2/(2.0*kappa*alpha1)) * np.sin(2.0*lambda_-tomega-gomega)\
		- e*((3.0 * chi2)/(2.0*kappa*alpha2)) * np.sin(tomega-gomega))

# eq 5:
	vr = a * kappa *\
		( e*np.sin(lambda_-tomega)\
		+ e**2 * (eta2/kappa2)*np.sin(2.0*(lambda_-tomega))\
		- i**2 * ((chi2/(2.0*alph2))*(nu/kappa)*np.sin(2.0*(lambda_-gomega))))

# eq 6:
	vL = n * ( 1.0 + 2.0 * e * np.cos(lambda_-tomega)\
		+ e**2 * (7.0/2.0 - 3.0*eta2/kappa2 - kappa2/(2.0*n**2)\
			+ (3.0/2.0 + eta2/kappa2) * np.cos(2.0*(lambda_-tomega)))\
		+ i**2 * (2.0 - kappa2/(2.0*n**2) - (3.0/2.0)*(chi2/kappa2)\
			- (chi2/(2.0*alph2)) *np.cos(2.0*(lambda_-gomega))))

# eq 7:
	vz = a * i * nu * ( np.cos(lambda_-gomega)\
		+ e * ((chi2*(kappa+nu))/(2.0*kappa*alpha1*nu)) *\
			cos(2.0*lambda_ - tomega - gomega)\
		+ e * ((3.0/2.0) * (chi2*(kappa-nu)) /(kappa*alpha2*nu)) *\
			cos(tomega-gomega) )

# convert to state vector

	x = r * np.cos(L)
	y = r * np.sin(L)
#	z = z, defined above
	vx = vr*np.cos(L) - r*vL*np.sin(L)
	vy = vr*np.sin(L) + r*vL*np.cos(L)
#	vz = vz, defined above
	sv= [-999.0,x,y,z,vx,vy,vz]
	svOUTPUT = sv[1:6]
	
	return svOUTPUT
# ************************ end of code by rfrench

def posi_geo(GM,R_eq,J2,J4,J6,svINPUT,verbose=True):#,elemOUTPUT # Stefan Renner code
# ## #	input: GM of central body (km**3/sec**2)
# #	       R_eq equatorial radius of central body (km)
# #	       J2,J4,J6 potential coeff
# #              svINPUT[0:5] position (km) and velocity (km/sec)
# #	NB - to minimize the possibility of index problems, 
# #	       this is padded to 
# #	       sv[1:6] position (km) and velocity (km/sec) 
# ## #	output: a semi major axis (km) 
# #               e eccentricity 
# #               i INCLINATION (deg) 
# #		gomega ascending node (deg) 
# #		tomega periapsis longitude (deg)
# #               lambda mean longitude (deg)		
# #	NB - to minimize the possibility of index problems, 
# #	       this is padded to elem[1:6] and returned as elemOUTPUT
# #		[0:5]
# ## #	NB: iterative computation because n, kappa,... a priori unknown
# #       
# #	implicit none
	
# #	double precision GM,R_eq,J2,J4,J6
# #	double precision sv(6),elem(6)
# #	double precision n,pperi 
# #	double precision pi,ppi,rad,epsilon
# #	double precision x,y,z,vx,vy,vz
# #	double precision r,L,vr,vL,arctan	
# #	double precision a,e,i,pnoeu,eta2,chi2,kappa,nu
# #	double precision alpha1,alpha2,alph2
# #	double precision a0,rcor,Lcor,vrcor,vLcor,zcor,vzcor
# #	double precision lambda,aaa,tomega,bbb,gomega
# #	double precision amodulo	
# #	double precision hz,r0,r0c,n0,n02

    if verbose:
        print('posi_geo TP1: GM,R_eq,J2,J4,J6,svINPUT',GM,R_eq,J2,J4,J6,svINPUT)
#    sv = [-999.0,svINPUT] # to give fortran-style indexing
#    sv = np.zeros(7)
    #	integer iter
    
    #	pi= dacos(-1.0)
    # pi=!dpi
    # ppi=2.0*pi
    # rad= pi/180.0
    # epsilon= 1.d-7		# precision on semi major axis (km)
    epsilon= 1.e-8		# precision on semi major axis (km)
    
    # x= sv[1]
    # y= sv[2]
    # z= sv[3]
    # vx= sv[4]
    # vy= sv[5]
    # vz= sv[6]
    x= svINPUT[0]
    y= svINPUT[1]
    z= svINPUT[2]
    vx= svINPUT[3]
    vy= svINPUT[4]
    vz= svINPUT[5]
    if verbose:
        print('posi_geo TP2: x,y,z,vx,vy,vz',x,y,z,vx,vy,vz)
    
    hz=x*vy-y*vx
     
    r= np.sqrt( x**2 + y**2 )
    L= np.arctan2(y,x) # arctan(X,Y)
    vr=   vx*np.cos(L) + vy*np.sin(L) 
    vL= (-vx*np.sin(L) + vy*np.cos(L))/r
                  
    #  initial computation of a, e, I, n, kappa, nu  
    a=r
    e=0.0
    i=0.0
    
    n = moyen_mouvement(GM,R_eq,J2,J4,J6,a,e,i)
    pperi,pnoeu = precess(GM,R_eq,J2,J4,J6,a,e,i)
    eta2,chi2 = freq(GM,R_eq,J2,J4,J6,a,e,i)

    if verbose:
        print('posi_geo TP3: n,pperi,pnoeu,eta2,chi2:',n,pperi,pnoeu,eta2,chi2)

    kappa= n - pperi
    nu =n - pnoeu
    
    alpha1= (1.0/3.0)*(2.0*nu + kappa)
    alpha2=2.0*nu - kappa
    alph2 = alpha1*alpha2
    
    a0=0.0
    iter_=0
    rcor=0.0
    Lcor=0.0
    vrcor=0.0
    vLcor=0.0
    zcor=0.0
    vzcor=0.0

	
# iterative computation of geometric elements
    while True:
        
        if verbose:
            print('posi_geo TP4: a,a0,a-a0,epsilon',a,a0,a-a0,epsilon)
        a= ( r - rcor )/( 1.0 - (vL - n - vLcor)/(2.0*n) )
        
        e= np.sqrt( ((vL - n - vLcor)/(2.0*n))**2 +  
                  ((vr - vrcor)/(a*kappa))**2 )
        
        lambda_= L - Lcor - 2.0*(n/kappa)*(vr -  vrcor)/(a*kappa)
        
        aaa= arctan( 1.0 - (r - rcor)/a , (vr - vrcor)/(a*kappa) )
        tomega= lambda_ - aaa
        
        i = np.sqrt( ( (z-zcor)/a )**2 + ( (vz-vzcor)/(a*nu) )**2 )
        
        bbb= arctan( (vz - vzcor)/nu , z-zcor )
        gomega=lambda_ - bbb
        
        rcor= a*(e**2)*( (3.0/2.0)*(eta2/(kappa**2)) - 1.0 -\
              0.50*(eta2/(kappa**2))*np.cos(2.0*aaa) )
        
        rcor= rcor + a*(i**2)*( (3.0*chi2)/(4.0*(kappa**2)) - 1.0\
                + (chi2/(4.0*alph2))*np.cos(2.0*bbb) )
         
        Lcor=  (3.0/4.0 + eta2/(2.0*kappa**2))*\
                    (n/kappa)*(e**2)*np.sin(2.0*aaa) 
        
        Lcor = Lcor - (1.0/4.0)*\
                (chi2/alph2)*(n/nu)*(i**2)*np.sin(2.0*bbb)
        
        vrcor= a*kappa*( (eta2/(kappa**2))*(e**2)*np.sin(2.0*aaa) )
        
## [vr - a*(I**2)...] correct (check formula (65) of Borderies 1994, + sign wrong) 	
        vrcor= vrcor - a*(i**2)*((chi2*nu)/(2.0*alph2))*np.sin(2.0*bbb)

        vLcor= n*(e**2)*( 7.0/2.0 - 3.0*(eta2/(kappa**2)) -\
      		     (kappa**2)/(2.0*n**2) +\
      		(3.0/2.0 + eta2/(kappa**2))*np.cos(2.0*aaa) )
     
        vLcor = vLcor + n*(i**2)*( 2.0 - (kappa**2)/(2.0*n**2)\
      		           - 1.50*(chi2/(kappa**2))\
      		           -0.50*(chi2/alph2)*np.cos(2.0*bbb) )
     
        zcor = a*i*e*(0.50*(chi2/(kappa*alpha1))*np.sin(aaa+bbb)\
            - 1.50*(chi2/(kappa*alpha2))*np.sin(bbb-aaa) )
     
        vzcor = a*i*e*(0.50*( (chi2*(kappa+nu))/(kappa*alpha1) )\
      			*np.cos(aaa+bbb)  
         + 1.50*( (chi2*(kappa-nu))/(kappa*alpha2) )\
      			*np.cos(bbb-aaa) ) 

# refine n, kappa, nu values
        n = moyen_mouvement(GM,R_eq,J2,J4,J6,a,e,i)
        pperi,pnoeu = precess(GM,R_eq,J2,J4,J6,a,e,i)
        eta2,chi2 = freq(GM,R_eq,J2,J4,J6,a,e,i)
             	
        kappa= n - pperi
        nu= n - pnoeu
        
        alpha1= (1.0/3.0)*(2.0*nu + kappa)
        alpha2=2.0*nu - kappa
        alph2 = alpha1*alpha2
        # THIS DOESN'T work
        # if (np.abs(a-a0) < epsilon):
        #     break
        a0=np.copy(a)
        iter_=iter_+1

        if iter_ > 100:  # can't get the above inequality to work - bails out right away
            break
        
    if verbose:
        print('posi_geo: converged TP5 n,pperi,pnoeu,eta2,chi2',n,pperi,pnoeu,eta2,chi2)

# converged	

    gomega= (np.degrees(gomega) + 360.0)%360.
    tomega= (np.degrees(tomega) + 360.0)%360.
    lambda_= (np.degrees(lambda_)+ 360.0)%360.
    if verbose:
        print('posi_geo: TP6 iter_,gomega,tomega,lambda_',iter_,gomega,tomega,lambda_)

#c		elem(1) = a
    elem = np.zeros(7)
    elem[2] = e
    elem[3] = np.degrees(i)
    elem[4] = amodulo(gomega,360.0)
    elem[5] = amodulo(tomega,360.0)
    elem[6] = amodulo(lambda_,360.0)

	
# semi major axis deduced from vertical angular momentum  
#STEP2:
    r0 = r
    n02 = freq_n0(GM,R_eq,J2,J4,J6,r0)
    n0 = np.sqrt(n02)
    
    r0c=0.0
    iter_=0
	
#STEP3:

    while True:#(np.abs(r0c-r0) >= epsilon):
        r0=np.sqrt(hz/n0)
        # if np.abs(r0c-r0) < epsilon:   
        #     break
        
        #	if (np.abs(r0c-r0) < epsilon) then goto,STEP4
        
        r0c = r0
        iter_=iter_+1	  
        n02 = freq_n0(GM,R_eq,J2,J4,J6,r0)
        n0 = np.sqrt(n02)
        if verbose:
            print('posi_geo: TP7 iter_,r0',iter_,r0)
        if iter_ > 100:
#        if np.abs(r0c-r0) < epsilon:   
            break

#	goto,STEP3

#STEP4:
    elem[1]= r0*(1.0 + e*e + i*i)	
    #  	elem(1)= a # TEST **********
    # test differed only by .014 km 
    elemOUTPUT = elem[1:]
    
    return elemOUTPUT
    
def calc_titan_longitude(UTCstr,MEAN_LON=True,USE_KERNELS_IF_PRESENT=True):
    ETsec = spice.str2et(UTCstr)
    return calc_titan_longitude_ET(ETsec,MEAN_LON=MEAN_LON,USE_KERNELS_IF_PRESENT=USE_KERNELS_IF_PRESENT)
        
def calc_titan_longitude_ET(ET,MEAN_LON=True,USE_KERNELS_IF_PRESENT=True):
    kernels_dir = '/Volumes/dione_raid2/Research/kernels/'
    kernels         =  ['cpck16Apr2014.tpc','Saturn_ring_plane.tf','naif0012.tls','sat425.bsp']
    count = 0 # default is to load kernels
    if USE_KERNELS_IF_PRESENT:
        count =spice.ktotal('ALL')
    if count == 0:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernels_dir + kernel)
    	
    # find GM, J2, J4, J6
    GM = spice.bodvar(699,'GM',1)
    GM = GM[0]
    jcoef = spice.bodvar(699,'JCOEF',6)
    J2 = jcoef[1]
    J4 = jcoef[3]
    J6 = jcoef[5]
    RADII = spice.bodvar(699,'RADII_RS',3)#,RADII # gives 60330.
    R_eq = RADII[0]
    
#    ET = spice.str2et(UTCstr)
    
    londegtitan = 0.0
    meanlontitan = 0.0

    try:
        starg,ltime = spice.spkezr('titan',ET,'SATURN_RING_PLANE','NONE','Saturn')
    except:
        spice.furnsh(kernels_dir+'Saturn_ring_plane.tf')
        try:
            starg,ltime = spice.spkezr('titan',ET,'SATURN_RING_PLANE','NONE','Saturn')
        except:
            spice.furnsh(kernels_dir+'sat425.bsp')
            try:
                starg,ltime = spice.spkezr('titan',ET,'SATURN_RING_PLANE','NONE','Saturn')
            except:
                spice.furnsh(kernels_dir+'naif0012.tls')
                try:
                    starg,ltime = spice.spkezr('titan',ET,'SATURN_RING_PLANE','NONE','Saturn')
                except:
                    spice.furnsh(kernels_dir+'cpck16Apr2014.tpc') # last resort since another tpc may be expected
                    starg,ltime = spice.spkezr('titan',ET,'SATURN_RING_PLANE','NONE','Saturn')
        
    elements = posi_geo(GM,R_eq,J2,J4,J6,starg,verbose=False)

    a=elements[0]
    e=elements[1]
    i=np.radians(elements[2])
    londegtitan = np.degrees(np.arctan2(starg[1],starg[0]))
    londegtitan = (londegtitan + 360.0) % 360.0
    meanlontitan = elements[5] % 360.0
    
    n = moyen_mouvement(GM,R_eq,J2,J4,J6,a,e,i)

    varpi_deg	 	= elements[4] % 360
    pperi,pnoeu  = precess(GM,R_eq,J2,J4,J6,a,e,i)
    varpi_dot_degday	= pperi * spice.dpr() *  spice.spd() # deg/day
    mean_motion_degday	= n * spice.dpr() *  spice.spd() # deg/day 

    if MEAN_LON:
        return meanlontitan
    else:
        return londegtitan

def calc_dr_Titan(ETsec,londeg,rkm,MEAN_LON=True,min_dist_from_res=np.sqrt(382.)):
    '''
    Translated from function calc_titan_droff.pro, but with saturation value
    This is for a scalar value of ETsec, londeg, rkm only
    '''
    londeg_Titan = calc_titan_longitude_ET(ETsec,MEAN_LON=MEAN_LON,USE_KERNELS_IF_PRESENT=True)
    dlondeg = londeg - londeg_Titan
    
    Smode = 382. # km^2 Nicholson et al Icarus 241 (214) 373 Table 10
    r_res = 77861.5
        # PDN Email Jan 1, 2018: From Paper II (2014), Table 10 gives the predicted resonance location
        # as 77857.4 km (using Jacobson's 2008 gravity solution). The text (Sec 12)
        # gives a revised location, based on the fitted eccentricity of the Titan
        # ringlet, of 77861.5\pm0.2 km.  So this is probably what we should be
        # using. I used J2 - J12 for my prediction (see Table 8).
    dr_res = rkm - r_res
    dr_res_mag = np.max([min_dist_from_res,np.abs(dr_res)]) # to avoid singularity
    # see Nicholson et al for limiting validity of test particle model to sqrt(S) ~ 17  km for Titan resonance

    if dr_res > 0:
        dphideg = 180.0
    else:
        dphideg = 0.0
    dr_Titan     = (Smode/dr_res_mag) * np.cos(np.radians(dlondeg - dphideg))

    return dr_Titan

def calc_dr_Titan_vec(ETsecs,londegs,rkms,MEAN_LON=True,min_dist_from_res=np.sqrt(382.),NMAX=500):
    Npts = len(ETsecs)
    if Npts > NMAX:
        ETsecs_ = np.linspace(ETsecs[0],ETsecs[-1],NMAX)
        
        flondeg_of_et = interpolate.interp1d(ETsecs,londegs,                                          
                                             bounds_error=False,fill_value='extrapolate')
        londegs_ = flondeg_of_et(ETsecs_)
        
        frkm_of_et = interpolate.interp1d(ETsecs,rkms,                                          
                                             bounds_error=False,fill_value='extrapolate')
        rkms_ = frkm_of_et(ETsecs_)
    else:
        ETsecs_ = ETsecs
        londegs_ = londegs
        rkms_ = rkms
        
    dr_Titan_vec_ = []
    for ETsec,londeg,rkm in zip(ETsecs_,londegs_,rkms_):
        dr_Titan_vec_.append(calc_dr_Titan(ETsec,londeg,rkm,MEAN_LON=MEAN_LON,min_dist_from_res=min_dist_from_res))

    if Npts > NMAX:
        fdr_of_et = interpolate.interp1d(ETsecs_,dr_Titan_vec_,
                                         bounds_error=False,fill_value='extrapolate')
        dr_Titan_vec = fdr_of_et(ETsecs)
    else:
        dr_Titan_vec = dr_Titan_vec_
    return dr_Titan_vec

def calc_dr_Titan_from_taufile(taufilepath,MEAN_LON=True,min_dist_from_res=np.sqrt(382.),NMAX=1000):
    '''
    Compute radial offset due to Titan across the radial range of a taufile. Note that no rss_ringoccs routines are required here.
    '''
    f = taufilepath
    REFERENCE_TIME = f[f.index('RSS_2')+4:f.index('RSS_2')+12]
    REFERENCE_TIME = REFERENCE_TIME.replace('_','-')+'T00:00'
    try:
        ETref = spice.str2et(REFERENCE_TIME)
    except: # likely error is absence of leap-seconds kernel
        tls_file = global_path_to_kernels+'naif/CASSINI/kernels/lsk/naif0012.tls'
        spice.furnsh(tls_file)
        ETref = spice.str2et(REFERENCE_TIME)

    rkms,londegs,SPM_RETs = np.loadtxt(f,delimiter=',',unpack=True,usecols = [0,3,10])
    ETsecs = ETref + SPM_RETs # ring event times

    try:
        dr_Titans = calc_dr_Titan_vec(ETsecs,londegs,rkms,MEAN_LON=MEAN_LON,min_dist_from_res=min_dist_from_res,NMAX=NMAX)
    except: # likely error is absence of planetary constants file
        tpc_file = global_path_to_kernels+'naif/CASSINI/kernels/pck/cpck26Feb2009.tpc'
        spice.furnsh(tpc_file)
        dr_Titans = calc_dr_Titan_vec(ETsecs,londegs,rkms,MEAN_LON=MEAN_LON,min_dist_from_res=min_dist_from_res,NMAX=NMAX)

    return dr_Titans

