import numpy as np
### read in binary spectrogram file
def read_spectro(filename):
    # read in file
    input = open(filename,'r').readlines()
    # read in values for frequency and time arrays
    f_arr = [float(input[i].strip().split(':')[1]) for i in [0,1,2]]
    t_arr = [float(input[i].strip().split(':')[1]) for i in [3,4,5]]
    # construct frequency and time arrays
    freq = f_arr[0]+np.arange(int(f_arr[1]))*f_arr[2]
    time = t_arr[0]+np.arange(int(t_arr[1]))*t_arr[2]
    # read in remaining data, ring radius + spectrogram
    rho = []
    spectrogram = []
    for i in range(6,int(t_arr[1]+6)):
        # split on delimiter
        row = input[i].split(',')
        # radius is first element
        rho += [float(row[0])]
        # spectrogram is remaining elements
        spectrogram += [[float(s) for s in row[1:int(f_arr[1]+1)]]]
    return np.array(time),np.array(rho),np.array(freq),np.array(spectrogram).T
