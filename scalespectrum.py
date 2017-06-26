#! /usr/bin/env python
#

import numpy as np
import sys

def makespectfile(afile):
    """ reads in a file and returns np.arrays containing values for frequency and amplitude """
    x = []
    y = []
    with open(afile) as f:
        for line in f:
            if line.startswith('#'): continue
            (freq,flux) = line.split()
            x.append(float(freq))
            y.append(float(flux))
        return (np.asarray(x),np.asarray(y))

if __name__ == "__main__":

    if len(sys.argv) > 1:
        spectrum_file = sys.argv[1]
        freq,amp = makespectfile(spectrum_file)
        freq = freq / 1e9
        for i in range(len(freq)):
            print("%f %g" % (freq[i],amp[i]))
    
