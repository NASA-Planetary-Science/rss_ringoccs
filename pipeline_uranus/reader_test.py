import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.path.append('../')
import rss_ringoccs.reader.vgr_uranus_reader as vur
sys.path.remove('../')


rawpath = '/Volumes/sflury001/Research/TC2017/data/VoyagerUranus/data/raw/'#'../data/'

for pdir in ['I','E']:
    start = time.time()
    read_inst = vur.VGRUranusReader(rawpath,'E')
    end = time.time()
    print(round(end-start))
    print(read_inst.spm.shape,read_inst.I.shape,read_inst.Q.shape)
    pow = read_inst.I**2+read_inst.Q**2
    psmooth = np.convolve(pow,np.ones(1000),'same')
    plt.plot(read_inst.spm,psmooth,'-k',lw=1)
    plt.savefig('reader_test_'+pdir.lower()+'.png')
    plt.close()
