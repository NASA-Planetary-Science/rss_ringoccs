"""

get_freespace.py

Purpose: Calculate SPM times at which freespace regions are expected.

Revisions:
    2018 Nov 19 - jfong - original
"""

from .find_gaps import find_gaps
from scipy.interpolate import interp1d
import pdb
import numpy as np

def get_freespace(t_ret_spm_vals, year, doy, rho_km_vals, 
        phi_rl_deg_vals, t_oet_spm_vals, atmos_occ_spm_vals,
        niter=int(100), tolerance=0.001,
        t0=2454467.000000, kernels=None, split_ind=None):

    if split_ind:
        t_ret_ing, t_oet_ing, phi_rl_ing, rho_ing = split_chord_arr(
                t_ret_spm_vals, t_oet_spm_vals, atmos_occ_spm_vals,
                phi_rl_deg_vals, rho_km_vals,
                split_ind, '"INGRESS"')
        
        gaps_km_ing = get_freespace_km(t_ret_ing, year, doy, rho_ing, 
                phi_rl_ing)

        
        

        t_ret_egr, t_oet_egr, phi_rl_egr, rho_egr = split_chord_arr(
                t_ret_spm_vals, t_oet_spm_vals, atmos_occ_spm_vals,
                phi_rl_deg_vals, rho_km_vals,
                split_ind, '"EGRESS"')

        rho_to_spm_ing = interp1d(rho_ing, t_oet_ing, fill_value='extrapolate')
        rho_to_spm_egr = interp1d(rho_egr, t_oet_egr, fill_value='extrapolate')
        

        gaps_km_egr = get_freespace_km(t_ret_egr, year, doy, rho_egr, 
                phi_rl_egr)
        gaps_spm_ing = rho_to_spm_ing(gaps_km_ing)
        gaps_spm_egr = rho_to_spm_egr(gaps_km_egr)

        gaps_km = {'INGRESS': gaps_km_ing,
                'EGRESS': gaps_km_egr}
        gaps_spm = gaps_spm_ing.tolist() + gaps_spm_egr.tolist()

    else:
        t_ret_out, t_oet_out, phi_rl_out, rho_out = remove_blocked(
                t_oet_spm_vals, atmos_occ_spm_vals, t_ret_spm_vals,
                phi_rl_deg_vals, rho_km_vals)

        gaps_km = get_freespace_km(t_ret_out, year, doy, rho_out,
                phi_rl_out)

        rho_to_spm = interp1d(rho_out, t_oet_out, fill_value='extrapolate')



        gaps_spm = rho_to_spm(gaps_km).tolist()

        # reverse list for ingress occ so that gaps_spm in increasing order
        if (rho_out[1]-rho_out[0]) < 0:
            # reverse each list in gaps_spm
            gaps_spm_1 = [[x[1],x[0]] for x in gaps_spm]
            gaps_spm = gaps_spm_1[::-1]

    return gaps_km, gaps_spm



        

def split_chord_arr(t_ret_spm_vals, t_oet_spm_vals,
        atmos_occ_spm_vals, phi_rl_deg_vals, rho_km_vals, ind, profdir):

    if profdir == '"EGRESS"':
        t_ret_spm_split = t_ret_spm_vals[ind:]
        t_oet_spm_split = t_oet_spm_vals[ind:]
        phi_rl_deg_split = phi_rl_deg_vals[ind:]
        rho_km_split = rho_km_vals[ind:]

    if profdir == '"INGRESS"':
        t_ret_spm_split = t_ret_spm_vals[:ind]
        t_oet_spm_split = t_oet_spm_vals[:ind]
        phi_rl_deg_split = phi_rl_deg_vals[:ind]
        rho_km_split = rho_km_vals[:ind]

    t_ret_out, t_oet_out, phi_rl_out, rho_out = remove_blocked(
            t_oet_spm_split, atmos_occ_spm_vals, t_ret_spm_split,
            phi_rl_deg_split, rho_km_split)

    return t_ret_out, t_oet_out, phi_rl_out, rho_out


def remove_blocked(t_oet_spm_vals, atmos_occ_spm_vals, t_ret_spm_vals,
        phi_rl_deg_vals, rho_km_vals):

    t_oet_spm_vals = np.asarray(t_oet_spm_vals)
    atmos_occ_spm_vals = np.asarray(atmos_occ_spm_vals)
    t_ret_spm_vals = np.asarray(t_ret_spm_vals)
    phi_rl_deg_vals = np.asarray(phi_rl_deg_vals)
    rho_km_vals = np.asarray(rho_km_vals)

    npts = len(t_oet_spm_vals)
    mask_not_blocked = np.array([False for i in range(npts)])

    for i in range(npts):
        if t_oet_spm_vals[i] not in atmos_occ_spm_vals:
            mask_not_blocked[i] = True


    t_ret_out = t_ret_spm_vals[mask_not_blocked]
    t_oet_out = t_oet_spm_vals[mask_not_blocked]
    phi_rl_out = phi_rl_deg_vals[mask_not_blocked]
    rho_out = rho_km_vals[mask_not_blocked]

    return t_ret_out, t_oet_out, phi_rl_out, rho_out



    

def get_freespace_km(ret_spm, year, doy, rho_km, phi_rl_deg):
    add_fsp1 = [[133500.0, 133650.0]]
    add_fsp2 = [[1.37e5, 1.43e5]]
    add_fsp3 = [[1.43e5, 1.7e5]]
    # check rho does not cover fsp1 and fsp2
    if add_fsp1[0][0] < min(rho_km) or add_fsp1[0][1] > max(rho_km):
        add_fsp1 = []
    if add_fsp2[0][0] < min(rho_km) or add_fsp2[0][1] > max(rho_km):
        add_fsp2 = []
    if add_fsp3[0][0] < min(rho_km) or add_fsp3[0][1] > max(rho_km):
        add_fsp3 = []


    freespace_km = (find_gaps(ret_spm, year, doy, rho_km, phi_rl_deg)
            + add_fsp1 + add_fsp2 + add_fsp3) #[[133500.0, 133650.0]] + [[1.37e5, 1.7e5]])

    ## add rho_atmo to 74000 as a freespace gap, if it exists
    #if max(rho_atmo) in rho_km:
    #    freespace_km = [[max(rho_atmo), 7.4e4]] + freespace_km


    #if profdir == '"INGRESS"':
    #    freespace_km = np.flip(freespace_km, 0)

    return freespace_km
#----------------------------------------------------------------- scrapped
#def get_chord_ind(self, spm=None):
#
#    # if spm given, interpolate rho_dot_kms_values to given spm values
#    #if spm:
#    #    spm_to_rho_dot = splrep(self.t_oet_spm_vals, self.rho_dot_kms_vals)
#    #    rho_dot_kms_vals = splev(spm, spm_to_rho_dot)
#    #else:
#    #    rho_dot_kms_vals = self.rho_dot_kms_vals
#
#
#    # add 2 for the index of the first element after the sign change
#    ind = np.argwhere(np.diff(np.sign(self.rho_dot_kms_vals)))
#
#    #print('PLOX STOP HERE')
#    #pdb.set_trace()
#    if len(ind) > 1:
#        print('WARNING! ring radius changes direction twice!')
#        pdb.set_trace()
#    elif len(ind) == 0:
#        ind = None
#    else:
#        ind = int(ind[0]) + 2
#
#     self.ingress_spm = t_oet_spm_vals[0:ind+2]
#     self.egress_spm = t_oet_spm_vals[ind+2:]
#    self.split_ind = ind
#
#    return ind
#
#def get_freespace_spm(self):
#    ind = self.get_chord_ind()
#    freespace_spm = []
#
#    year = self.rev_info['year']
#    doy = self.rev_info['doy']
#    spm_to_rho = splrep(self.t_oet_spm_vals, self.rho_km_vals)
#    rho_atmo = splev(self.atmos_occ_spm_vals, spm_to_rho)
#
#    if ind is not None:
#        (t_ret_spm_ing, t_oet_spm_ing, phi_rl_deg_ing, rho_km_ing,
#                rho_to_spm_ing) = self.split_chord_arr(ind, '"INGRESS"')
#
#        rho_to_spm_ing= splrep(np.sort(rho_km_ing),
#                np.flip(t_oet_spm_ing,0))
#
#        freespace_km_ing = self.get_freespace_km(
#                t_ret_spm_ing, year, doy, rho_km_ing, phi_rl_deg_ing,
#                '"INGRESS"')
#
#
#        #if min(self.atmos_occ_spm_vals) in t_oet_spm_ing:
#        #    add_fsp_spm_ing = [[float(splev(7.4e4, rho_to_spm_ing)),
#        #        min(self.atmos_occ_spm_vals)]]
#        #else:
#        #    add_fsp_spm_ing = None
#
#        for ni in freespace_km_ing:
#
#
#            if min(rho_km_ing)> ni[0] or min(rho_km_ing)>ni[1]:
#                print('gr')
#                break
#
#            fsp = [float(splev(ni[1], rho_to_spm_ing)),
#                float(splev(ni[0], rho_to_spm_ing))]
#            freespace_spm.append(fsp)
#
#        #if add_fsp_spm_ing:
#        #    freespace_spm = freespace_spm + add_fsp_spm_ing
#
#
#
#
#        (t_ret_spm_egr, t_oet_spm_egr, phi_rl_deg_egr, rho_km_egr,
#                rho_to_spm_egr) = self.split_chord_arr(ind, '"EGRESS"')
#
#        #pdb.set_trace()
#
#        rho_to_spm_egr= splrep(rho_km_egr, t_oet_spm_egr)
#
#        freespace_km_egr = self.get_freespace_km(
#                t_ret_spm_egr, year, doy, rho_km_egr, phi_rl_deg_egr,
#                '"EGRESS"')
#        #if min(self.atmos_occ_spm_vals) in t_oet_spm_egr and (rho_km_egr <7.4e4).all():
#        #    add_fsp_spm_egr= [[min(self.atmos_occ_spm_vals),
#        #        float(splev(7.4e4, rho_to_spm_egr))]]
#        #else:
#        #    add_fsp_spm_egr = None
#
#        for ne in freespace_km_egr:
#            print(ne)
#            if min(rho_km_egr)>ne[0]:
#                print('gr')
#                continue
#            fsp = [float(splev(ne[0], rho_to_spm_egr)),
#                float(splev(ne[1], rho_to_spm_egr))]
#            freespace_spm.append(fsp)
#        #if add_fsp_spm_egr:
#        #    freespace_spm = add_fsp_spm_egr + freespace_spm
#
#        # remove all spm gap ranges that include atmosphere
#        # note that each spm limit is from smallest spm to largest
#        #fsp_copy = freespace_spm.copy()
#        #pdb.set_trace()
#        #for ex in fsp_copy:
#        #    print(ex)
#        #    if (self.atmos_occ_spm_vals>ex[0]).any() and (self.atmos_occ_spm_vals<ex[1]).any():
#        #        print('found one!')
#        #        freespace_spm.remove(ex)
#
#        #pdb.set_trace()
#
#
#
#    elif self.rev_info['prof_dir'] == '"INGRESS"':
#        freespace_km = self.get_freespace_km(
#                self.t_ret_spm_vals, year, doy, self.rho_km_vals,
#                self.phi_rl_deg_vals, self.rev_info['prof_dir'])
#        rho_to_spm = splrep(np.sort(self.rho_km_vals),
#                np.flip(self.t_oet_spm_vals,0))
#        if min(self.atmos_occ_spm_vals) in self.t_oet_spm_vals:
#            add_fsp_spm = [[float(splev(7.4e4, rho_to_spm)),
#                min(self.atmos_occ_spm_vals)]]
#        else:
#            add_fsp_spm = None
#
#        #print(freespace_km)
#        #pdb.set_trace()
#
#        for n in freespace_km:
#            if min(self.rho_km_vals)> n[0]:
#                print('gr')
#                pdb.set_trace()
#                break
#
#            fsp = [float(splev(n[1], rho_to_spm)),
#                float(splev(n[0], rho_to_spm))]
#            print(fsp)
#            freespace_spm.append(fsp)
#        if add_fsp_spm:
#            freespace_spm = freespace_spm + add_fsp_spm
#
#    elif self.rev_info['prof_dir'] == '"EGRESS"':
#        freespace_km = self.get_freespace_km(
#                self.t_ret_spm_vals, year, doy, self.rho_km_vals,
#                self.phi_rl_deg_vals, self.rev_info['prof_dir'])
#                #rho_atmo)
#        self.freespace_km = freespace_km
#        rho_to_spm = splrep(self.rho_km_vals, self.t_oet_spm_vals)
#        if max(self.atmos_occ_spm_vals) in self.t_oet_spm_vals:
#                            add_fsp_spm = [[max(self.atmos_occ_spm_vals),
#                float(splev(7.4e4, rho_to_spm))]]
#        else:
#            add_fsp_spm = None
#
#        #print(freespace_km)
#        #pdb.set_trace()
#
#        for n in freespace_km:
#            if min(self.rho_km_vals)> n[0]:
#                print('gr')
#                pdb.set_trace()
#                continue
#
#            fsp = [float(splev(n[0], rho_to_spm)),
#                float(splev(n[1], rho_to_spm))]
#            print(fsp)
#            freespace_spm.append(fsp)
#
#        if add_fsp_spm:
#            freespace_spm = add_fsp_spm + freespace_spm
#
#
#    #print(freespace_spm)
#    #pdb.set_trace()
#        return freespace_spm
