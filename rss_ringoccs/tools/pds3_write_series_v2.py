
'''
pds3_write_series.py

Purpose: Functions for getting PDS3 label file information and writing PDS3
         label files.

'''
import numpy as np
import time

def pds3_write_series_lbl(str_lbl, out_lbl_file):

    all_kwd = str_lbl['keywords_values']

    n_kv_all = len(all_kwd)

    cr = '\r\n'
    blank_line = cr
    eq = '= '


    npad0 = str_lbl['alignment_column']
    npad1 = str_lbl['series_alignment_column']

    f = open(out_lbl_file, 'w')

    # Write keywords and values
    for n in range(n_kv_all):
        kwd_list = all_kwd[n]
        if kwd_list == 'BLANK':
            f.write(blank_line)
        else:
            key_order = kwd_list['key_order']
            for key in key_order:
                this_kwd = key
                this_val = kwd_list[this_kwd]
                f.write(this_kwd.ljust(npad0) + eq + this_val + cr)

    # Write NAIF toolkit version
    naif_all_kwd = str_lbl['keywords_NAIF_TOOLKIT_VERSION']
    naif_key_order = naif_all_kwd['key_order']
    naif_ver_kwd = naif_key_order[0].ljust(npad0)
    naif_ver_val = naif_all_kwd[naif_key_order[0]]

    if naif_ver_val != '':
        f.write(naif_ver_kwd + eq + naif_ver_val + cr)
    
    # Write kernels SPICE file names
    #   Add braces to first and last kernel
    kernel_kwd = naif_key_order[1]
    kernel_val_list = naif_all_kwd[kernel_kwd]

    if kernel_val_list != '':
        if isinstance(kernel_val_list, list) is False:
            kernel_val_list = [kernel_val_list]
        kernel_val_list[0] = '{' + kernel_val_list[0]
        kernel_val_list[-1] = kernel_val_list[-1] + '}'
        kern_pad = ' ' * (npad0+3)

        nkern = len(kernel_val_list)
        for i in range(nkern):
            if i == 0:
                if nkern == 1:
                    app = ''
                else:
                    app = ','
                f.write(kernel_kwd.ljust(npad0) + eq + kernel_val_list[i]
                        + app + cr)
            elif i == (nkern-1):
                f.write(kern_pad + kernel_val_list[i] + cr)
            else:
                f.write(kern_pad + kernel_val_list[i] + ',' + cr)

    f.write(blank_line)

    # Write file description
    desc_key = 'DESCRIPTION'.ljust(npad0)
    desc_val = str_lbl['description']
    sd = str_lbl['string_delimiter']
    
    desc_lines = desc_val.split(sd)
    nd_lines = len(desc_lines)

    for dline in range(nd_lines):
        if dline == 0:
            f.write(desc_key + eq + desc_lines[dline] + cr)
        else:
            f.write(desc_lines[dline] + cr)

    f.write(blank_line)

    # Write processing history description
    hist_key = 'PROCESSING_HISTORY_TEXT'.ljust(npad0)
    if len(str_lbl['history'])!=0:
        hist_desc_val = str_lbl['history']['description']
    
        hist_desc_lines = hist_desc_val.split(sd)
        hist_desc_nlines = len(hist_desc_lines)
        for hdline in range(hist_desc_nlines):
            if hdline == 0:
                f.write(hist_key + eq + '"' +  hist_desc_lines[hdline] + cr)
            else:
                f.write(hist_desc_lines[hdline] + cr)
    
        # Loop through all history
        f.write(blank_line)
        h_indent = '    '
        hpad1 = 18    
        hpad2 = 24
        hist_dict = str_lbl['history']
        hist_user_keys = str_lbl['history']['key_order0']
    
        # Write user history
        f.write('User history:' + cr)
        for ukey in hist_user_keys:
            f.write(h_indent + ukey.ljust(hpad2) + eq + hist_dict[ukey] + cr)
    
        hist_dict = str_lbl['history']
        hist_dict_keys = str_lbl['history']['key_order1']
        hist_list_keys = [str_lbl['history']['hist name']]
        hist_list_kwd = ['']
    
    
        for hkey in hist_dict_keys:
            this_kwd = hkey
            this_val = hist_dict[hkey]
            if (this_kwd == 'Positional Args'): 
                for varkey in this_val:
                    this_input = hist_dict[this_kwd][varkey]
                    if isinstance(this_input, dict):
                        hist_list_keys.append(varkey)
                        hist_list_kwd.append(this_kwd)
    
      # check if hist_list_keys covers all input instances from start
      #     of processing to this point in history
       
        n_extrahist = len(hist_list_keys)
        for extra in range(1,n_extrahist):
            exkey = hist_dict[hist_list_kwd[extra]][hist_list_keys[extra]][
                    'Keyword Args']
            exvar = hist_dict[hist_list_kwd[extra]][hist_list_keys[extra]][
                    'Positional Args']
            for exkeyval in exkey:
                if isinstance(exkey[exkeyval], dict):
                    if exkeyval not in hist_list_keys:
                        hist_list_keys.append(exkeyval)
                        hist_dict['Keyword Args'][exkeyval] = exkey[exkeyval]
            for exvarval in exvar:
                if isinstance(exvar[exvarval], dict):
                    if exvarval not in hist_list_keys:
                        hist_list_keys.append(exvarval)
                        hist_dict['Positional Args'][exvarval] = exvar[exvarval]
    
        n_extrahist_final = len(hist_list_keys)
    
    
    #        else:
        # Loop over all other instance histories
        write_history_text(f, hist_list_keys, hist_dict, hist_dict_keys)


    f.write('"' + cr)
    f.write(blank_line)



    # Write series info
    #   print keyword/value pairs, parsing value to see if it is a packed string
    #   that contains a string_delimiter marking the end of a line
    s_set = str_lbl['keywords_series']
    series_all_kwd = str_lbl['keywords_series']
    series_key_order = series_all_kwd['key_order']
    indent = '  '

    for skey in series_key_order:
        this_key = skey
        this_val = series_all_kwd[skey]
        if this_key == 'OBJECT':
            f.write(this_key.ljust(npad1+2) + eq + this_val + cr)
        elif this_val.find(sd) == -1:
            f.write(indent + this_key.ljust(npad1) + eq + this_val + cr)
        else:
            lines = this_val.split(sd)
            nlines = len(lines)
            for d in range(nlines):
                if d==0:
                    f.write(indent + this_key.ljust(npad1) + eq + lines[d]
                            + cr)
                else:
                    f.write(indent + lines[d] + cr)




    # Loop over all objects
    obj_keys = str_lbl['object_keys']
    obj_vals = str_lbl['object_values']

    
    nobjs = len(obj_vals[0])
    nkeys = len(obj_keys)
    indent = '    ' #* nkeys
    npad_obj = npad1-2

    for obj in range(nobjs):
        f.write(blank_line)
        objcol = '  ' + 'OBJECT'.ljust(npad1) + '= COLUMN' + cr 
        f.write(objcol)
#        f.write('  OBJECT                            = COLUMN' + cr)
        for k in range(nkeys):
            this_key = obj_keys[k].ljust(npad1-2)
            this_val = str(obj_vals[k][obj])
            if this_val != '':
                if this_val.find(sd) == -1:
                    f.write(indent + this_key + eq + this_val + cr)
                else:
                    lines = this_val.split(sd)
                    nlines = len(lines)
                    for d in range(nlines):
                        if d == 0:
                            f.write(indent + this_key + eq + lines[d] + cr)
                        else:
                            f.write(indent + lines[d] + cr)
        objend = '  ' + 'END_OBJECT'.ljust(npad1) + '= COLUMN' + cr
        f.write(objend)

    f.write(blank_line)
    f.write('END_OBJECT'.ljust(npad1+2) + eq + 'SERIES' + cr)
    f.write(blank_line)
    f.write('END' + cr)


    f.close()

    
    return None

def write_history_text(f, hist_list_keys, hist_dict, hist_dict_keys):
    cr = '\r\n'
    blank_line = cr
    eq = '= '
    h_indent = '    '
    hpad = 18
    n_histories = len(hist_list_keys)

    hists_written = []

    # Loop over all histories
    for ehist in range(n_histories):
        f.write(blank_line)
        this_inst = hist_list_keys[ehist]
        if ehist == 0:
            f.write(str(this_inst) + ':' + cr)
            this_hist = hist_dict
        else:
            f.write(str(this_inst) + ' history:' + cr)
            this_hist = hist_dict['Positional Args'][this_inst]
        # Loop over all keys in hist_dict_keys
        for hkey in hist_dict_keys:
            this_kwd = hkey
            this_val = this_hist[this_kwd]
            if this_val == '':
                continue
            n_thisval = len(this_val)
            if (this_kwd == 'Positional Args') or (this_kwd == 'Keyword Args') or (this_kwd == 'Additional Info'):
                n_input = 0
                # Loop over input var/kwd dictionary keys
                for varkey in this_val:
                    this_input = this_val[varkey]
                    if isinstance(this_input, dict):
                        #if (ehist == n_histories-1):
                        if str(varkey) in hists_written:
                            pstr = 'see history above'
                        else:
                            pstr = 'see history below'
                        if n_input == 0:
                            f.write(h_indent + this_kwd.ljust(hpad) + eq
                                    + '{' + str(varkey)+ ': ' + pstr + ',' + cr)
                        elif n_input == len(this_val)-1:
                            f.write(h_indent + ''.ljust(hpad+3) + str(varkey)
                                    + ': ' + pstr + '}' + cr)
                        else:
                            f.write(h_indent + ''.ljust(hpad+3) + str(varkey)
                                    + ': ' + pstr + ',' + cr)
                        n_input = n_input + 1

                    else: # if not a dictionary
                        space_avail = 40

                        nchars = len(str(varkey)) + 1 + len(str(this_input))
                        if n_input == 0:
                            if n_input == len(this_val)-1:
                                app = '}'
                            else:
                                app = ','
                            if nchars > space_avail:
                                history_loop_over_long_inputs(f, this_input,
                                        this_kwd, varkey, 'first')
                            else:
                                f.write(h_indent + this_kwd.ljust(hpad)
                                        + eq + '{' + str(varkey) + ': '
                                        + str(this_input) + app + cr)
                        elif n_input == len(this_val)-1:
                            if nchars > space_avail:
                                history_loop_over_long_inputs(f, this_input,
                                        this_kwd, varkey, 'last')
                            else:
                                f.write(h_indent + ''.ljust(hpad+3)
                                        + str(varkey) + ': '
                                        + str(this_input) + '}' + cr)
                        else:
                            if nchars > space_avail:
                                history_loop_over_long_inputs(f, this_input,
                                        this_kwd, varkey, '')
                            else:
                                f.write(h_indent + ''.ljust(hpad+3)
                                        + str(varkey) + ': '
                                        + str(this_input) + ',' + cr)
                        n_input = n_input+1
            else:
                f.write(h_indent + this_kwd.ljust(hpad) + eq + this_val + cr)
        hists_written.append(str(this_inst))
    return None

def history_loop_over_long_inputs(f, this_input, this_kwd, varkey, num):
    h_indent = '    '
    hpad = 18
    eq = '= '
    cr = '\r\n'
    blank_line = cr

    # split input into several lines
    kline = str(this_input).split(',')
    n_lines = len(kline)
    varlen = len(str(varkey)) + 2

    if varkey == 'freespace_spm':
        kline[0] = kline[0].replace('\n', ',\r\n'+' '*(hpad+varlen+7))
        #kline = (str(this_input)).replace('\n', '\n'+' '*hpad)

    nval = 0
    for nl in range(n_lines):
        if nval == 0 and num == 'first':
            if varkey == 'freespace_spm':
                app = ''
            else:
                app = ','
            f.write(h_indent + this_kwd.ljust(hpad) + eq + '{' + str(varkey)
                    + ': ' + kline[nl] + app + cr)
        elif nval == 0 and num != 'first':
            f.write(h_indent + ''.ljust(hpad+3) + str(varkey) + ': '
                    + kline[nl] + ',' + cr)
        elif nval == n_lines-1 and num == 'last':
            f.write(h_indent + ''.ljust(hpad+3+varlen) + kline[nl] + '}' + cr)
        else:
            f.write(h_indent + ''.ljust(hpad+3+varlen) + kline[nl] + ',' + cr)
        nval = nval+1

    return None

def get_ring_obs_id(year, doy, band, dsn):
    band = band[1]
    dsn = dsn.split('-')[-1]
    out_str = '"S/OCC/CO/RSS/' + year + '-' + doy + '/' + band + dsn + '"'
    return out_str

def get_ring_profile_direction(rho_vals):
    r = rho_vals
    dr_start = r[1] - r[0]
    dr_end = r[-1] - r[-2]
    
    if dr_start < 0 and dr_end < 0:
        occ_dir = '"INGRESS"'
        occ_char = 'I'
    if dr_start > 0 and dr_end > 0:
        occ_dir = '"EGRESS"'
        occ_char = 'E'
    if (dr_start < 0 and dr_end > 0) or (dr_start > 0 and dr_end < 0):
        occ_dir = '"BOTH"'
        occ_char = 'C'
    return occ_dir, occ_char

        
def get_record_bytes(fmts):
    bytes_list = []
    ncols = len(fmts)
    for n in range(ncols):
        fmt = fmts[n]
        ind1 = 2
        ind2 = fmt.index('.')
        bytes_list.append(int(fmt[ind1:ind2]))

    bytes_list_wc = [x+1 for x in bytes_list]
    bytes_list_wc.insert(0,1)
    start_bytes_list = (list(np.cumsum(bytes_list_wc)))[:-1]

    # ncols-1 for number of bytes used for commas
    # additional 2 for special characters (new line + carriage return)
    record_bytes = str(np.sum(bytes_list) + (ncols-1) + 2)

    return record_bytes, bytes_list, start_bytes_list


        


def get_ISOD_str(spm_start, year, doy):
    # sample START_TIME                      = 2005-123T07:40:00.000

    hrs_float = spm_start/60./60.
    mins_float = (hrs_float % 1) * 60.
    secs_float = (mins_float % 1) * 60.

    hr_start_str = str(int(hrs_float)).zfill(2)
    min_start_str = str(int(mins_float)).zfill(2)
    sec_start_str = str(int(secs_float)).zfill(2)

    out_str = (year + '-' + doy + 'T' + hr_start_str + ':' + min_start_str + ':' + sec_start_str)
    return out_str

def get_sampling_interval(sampling_param):
    dr_start = sampling_param[1] - sampling_param[0]
    dr_end = sampling_param[-1] - sampling_param[-2]
    err = 0.001
    dr_end_bounds = [dr_end - err, dr_end + err]
    if dr_start > dr_end_bounds[0] and dr_start < dr_end_bounds[1]:
        dr = dr_start
    else:
        print('WARNING (pds3_write_series_v2): Sampling interval is not constant!\n\tdr_start, dr_end = ', dr_start, dr_end, ' Using dr_start.')

#    return str(dr)
    return f'{dr:0.4f}'

def get_spm_ref_time(year, doy):
    return year + '-' + doy + 'T00:00:00'


"""
Revisions:
    2018 Jul 23 - jfong - copied from jwf_pds3_write_series_v3.py
                        - copy over write_history_text()
    2018 Jul 27 - jfong - get rid of comma for kernels if only one listed
                        - debug history loop (extra instance)
                        - add error range for constant sampling interval check
    2018 Sep 20 - jfong - remove kernels str split when varkey=='kernels'
    2018 Nov 13 - jfong - redo processing history to have user history only
                          once
"""
