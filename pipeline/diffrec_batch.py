import e2e_batch_args as args
import sys, os
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import traceback
import time
import numpy as np

files = [args.mpath+line.strip('\n') for line in open(
                args.rsr_file_list,'r').readlines()]
path = '../output/'
occdirs = [path+dir1+'/'+dir2+'/'+dir3+'/'
            for dir1 in os.listdir(path)
            if '.' not in dir1
            for dir2 in os.listdir(path+dir1)
            if '.' not in dir1+'/'+dir2
            for dir3 in os.listdir(path+dir1+'/'+dir2)
            if '.' not in dir1+'/'+dir2+'/'+dir3]

err_file = sys.argv[0].split('.')[0] + time.strftime("_%Y%m%d-%H%M%S") + '.err'
fail_file = open('../output/' + err_file, 'w')

def filename(file,type):
    # directory
    dir = file.split('/RSS')[0]+'/'
    file = file[len(dir):]
    # split up file name
    splts = np.array([splt for splt in file.split('_')],dtype=str)
    # find corresponding file
    file = "_".join(splt for splt in splts[:5])+'_'+type+'_'+"_".join(splt for splt in splts[7:])
    if not os.path.exists(dir+file):
        if splts[-1][-5] == '1':
            splts[-2] = str(int(splts[-2])-1)
        else:
            splts[-1] = '000'+str(int(splts[-1][-5])-1)+'.TAB'
        file = "_".join(splt for splt in splts[:5])+'_'+type+'_'+"_".join(splt for splt in splts[7:])

    return dir+file

read_time = 0.0
diffcorr_time = 0.0

N = int(0)
for n,dir in enumerate(occdirs):

    for m,file in enumerate(os.listdir(dir)):

        if 'DLP' in file and 'TAB' in file:

            try:
                N += 1
                print('Directory ',n,'File ',N)

                # copy file name for dlp file
                dlp_file = dir+"".join(f for f in file)
                # find corresponding geo file
                geo_file = filename(dlp_file,'GEO')
                # find corresponding cal file
                cal_file = filename(dlp_file,'CAL')

                # build DLP instance
                t1 = time.time()
                dlp_inst = rss.tools.ExtractCSVData(geo_file, cal_file, dlp_file,
                        verbose=args.verbose)
                t2 = time.time()
                read_time += t2-t1
                # diffraction reconstruction
                t1 = time.time()
                tau_inst = rss.diffrec.DiffractionCorrection(
                    dlp_inst, args.res_km,
                    rng=args.inversion_range, res_factor=args.res_factor,
                    psitype=args.psitype, wtype=args.wtype, fwd=args.fwd,
                    norm=args.norm, bfac=args.bfac, write_file=True,
                    verbose=args.verbose)
                t2 = time.time()

                diffcorr_time += t2-t1

            except KeyboardInterrupt:
                sys.exit()

            except:
                tb = traceback.format_exc()
                fail_file.write('-'*72+'\n')
                print(tb)
                fail_file.write('  '+dir+'\n'+'-'*72+'\n\n'+ 'n='+
                                 str(n)+'\n'+ tb+'\n')

print('-'*48)
print('Processing time in seconds for each step')
print('-'*48)
print('Read time:     ',round(read_time,3))
print('DiffCorr time: ',round(diffcorr_time,3))
tot_time = read_time + diffcorr_time
print('Total time:    ',round(tot_time,3))
print('-'*48)
