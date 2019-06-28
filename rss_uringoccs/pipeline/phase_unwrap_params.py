##################################################################
###
###	CONTENTS
###
###	  Complex radio signal phase unwrapping parameters
###	  for 3rd-order steering done using correct_phase
###	  in the uring_dlp.py script. These parameters are
###   stored in this script file in a single dict which
###   is imported by uring_dlp.py and indexed by the
###   occultation profile direction and ring name.
###
###	  Note: this is after 1st-order steering has
###	  been performed in calibration.
###
##################################################################
#
# PROFILE DIRECTIONS
#
#   I = ingress
#   E = egress
#
dirs = ['I','E']
#
# RING NAMES
#
#   6 = ring 6
#   5 = ring 5
#   4 = ring 4
#   A = ring alpha
#   B = ring beta
#   N = ring eta
#   G = ring gamma
#   D = ring delta
#   L = ring lambda
#   E = ring epsilon
#
rings = ['6','5','4','A','B','N','G','D','L','E']
#
# PARAMETERS
#
#   T_OET_START  -- the start time of the time interval
# 	                for correcting the unwrapping results
#   T_OET_END    -- the end time of the time interval
# 	                for correcting the unwrapping results
#   PHASE_RANGE  -- the maximum or minimum limit of phase
# 	                values to correct, depends on sign of 2PI_INCR
#   2PI_INCR     -- the integer increment of 2pi by which
# 	                to correct the phase
#
pars = ['T_OET_START','T_OET_END','PHASE_RANGE','2PI_INCR']
##################################################################
###  CREATE DICTIONARY FOR BOTH PROFILE DIRECTIONS AND ALL RINGS
##################################################################
params = {}
for dir in dirs:
    params[dir] = {}
    for ring in rings:
        params[dir][ring] = {}
        for par in pars :
            params[dir][ring][par] = []
##################################################################
###   SET PARAMETERS FOR EACH RING PROFILE
##################################################################
#
# INGRESS
#
#   6 -- no corrections necessary, unwrapping was sufficient
#   5 -- no corrections necessary, unwrapping was sufficient
#   4 -- no corrections necessary, unwrapping was sufficient
#   ALPHA
params['I']['A']['T_OET_START'] += [81610.3,81644.5,81648.7]
params['I']['A']['T_OET_END']   += [81644.5,81648.7,81750]
params['I']['A']['PHASE_RANGE'] += [300,292.5,285]
params['I']['A']['2PI_INCR']    += [-1,-1,1]
#   BETA
params['I']['B']['T_OET_START'] += [81512.3,81512.25,81610]
params['I']['B']['T_OET_END']   += [81610.4,81512.40,81640]
params['I']['B']['PHASE_RANGE'] += [4,509,1140]
params['I']['B']['2PI_INCR']    += [-1,1,-2]
#  ETA
params['I']['N']['T_OET_START'] += [81281.1]+11*[81294.9]+2*[81295.5]+[81315,81385,81440]
params['I']['N']['T_OET_END']   += [81295.5]+11*[81295.5]+2*[81315]+[81385,81450,81450]
params['I']['N']['PHASE_RANGE'] += [0]+11*[6.28]+2*[6.28]+[6.28,-1000,-650]
params['I']['N']['2PI_INCR']    += [1]+3*[1]+8*[-1]+[-5,-1]+[-5,-5,-1]
#   GAMMA
params['I']['G']['T_OET_START'] += [81211.7,81211.7,81212.0,81212.0,81281.0,81295,81295.04,81295.04,81295.04,81295.1,81295.1,81295.0,81295.0,81335]
params['I']['G']['T_OET_END']   += [81295.4,81295.4,81295.4,81295.4,81295.4,81295.0375,81295.075,81295.075,81295.075,81295.2,81295.2,81335,81335,81400]
params['I']['G']['PHASE_RANGE'] += [2,2,2,2,2,2,2,2,2,10,10,10,10,-1000]
params['I']['G']['2PI_INCR']    += [1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1]
#   DELTA
params['I']['D']['T_OET_START'] += 5*[81211.6]+[81220.0]+5*[81294.6]
params['I']['D']['T_OET_END']   += 5*[81295.4]+[81295.5]+5*[81295.6]
params['I']['D']['PHASE_RANGE'] += 5*[-5]+[-2]+4*[-5]+[3]
params['I']['D']['2PI_INCR']    += 5*[1]+[1]+4*[1]+[-1]
#   EPSILON
params['I']['E']['T_OET_START'] += 8*[80901]+26*[80905]
params['I']['E']['T_OET_END']   += 8*[81010]+26*[80910]
params['I']['E']['PHASE_RANGE'] += 8*[3.5]+16*[3.5]+10*[-2.5]
params['I']['E']['2PI_INCR']    += 8*[-1]+16*[-1]+10*[1]
#
# EGRESS
#
#   6
params['E']['6']['T_OET_START'] += [4610]
params['E']['6']['T_OET_END']   += [4747]
params['E']['6']['PHASE_RANGE'] += [-2]
params['E']['6']['2PI_INCR']    += [1]
#   5
params['E']['5']['T_OET_START'] += [4745.9]
params['E']['5']['T_OET_END']   += [6000]
params['E']['5']['PHASE_RANGE'] += [3]
params['E']['5']['2PI_INCR']    += [-1]
#   4 -- no corrections necessary, unwrapping was sufficient
#   ALPHA
params['E']['A']['T_OET_START'] += [5055.75,5055.85,5055.93,5055.94,5056.14,5145.,5120.]
params['E']['A']['T_OET_END']   += [6000]*7
params['E']['A']['PHASE_RANGE'] += [-12,-12,-12,-12,-10,-8,-8]
params['E']['A']['2PI_INCR']    += [1]*6+[-1]
#   BETA
params['E']['B']['T_OET_START'] += [5145.1,5217.7,5220.8]
params['E']['B']['T_OET_END']   += [5153.6,5244.3,5231]
params['E']['B']['PHASE_RANGE'] += [2,2,2]
params['E']['B']['2PI_INCR']    += [1]*3
#   ETA
params['E']['N']['T_OET_START'] += [5412.5]
params['E']['N']['T_OET_END']   += [6000]
params['E']['N']['PHASE_RANGE'] += [2]
params['E']['N']['2PI_INCR']    += [-1]
#   GAMMA
params['E']['G']['T_OET_START'] += [5375,5400]
params['E']['G']['T_OET_END']   += [6000]*2
params['E']['G']['PHASE_RANGE'] += [-2,-1]
params['E']['G']['2PI_INCR']    += [-1]*2
#   DELTA
params['E']['D']['T_OET_START'] += [5412.5,5492.4,5493.3]
params['E']['D']['T_OET_END']   += [5481.1,5493.3,5493.4]
params['E']['D']['PHASE_RANGE'] += [7.5,2,7.5]
params['E']['D']['2PI_INCR']    += [-1,1,-1]
#   EPSILON
params['E']['E']['T_OET_START'] += [5820]*20+[5800]*10#[5825.2,5855.5,5855.56,5855.58,5855.6,5855.6,5855.6,5855.6,5855.6,5856.9,5856.9]
params['E']['E']['T_OET_END']   += [6000]*20+[5875]*10#[5855.3,5855.56,5855.58,5855.6,5855.8,5856.4,5856.4,5856.4,5856.4,5858.9,5858.9]
params['E']['E']['PHASE_RANGE'] += [3]*20+[10]*10#[6,0,7,7,7,5,1,5,1,5,5]
params['E']['E']['2PI_INCR']    += [1]*20+[-1]*10#[-1,1,-1,-1,-1,-1,1,-1,1,-1,-1]
