#!/usr/bin/env python3
import os
import sys
import numpy as np

#read-in
laserfile = sys.argv[1]
laserend = float(sys.argv[2])

#read laser file
time, *rest, freq = np.loadtxt(laserfile, unpack=True)

#change laser file: set frequency to infinity (high value) to avoid any resonant laser hop
for nr,i in enumerate(time):
    if i > laserend:
        freq[nr] = 100.0

#write new laser file
np.savetxt(laserfile + '_new', np.transpose([time, *rest, freq]), fmt='%15.8E', delimiter=' ')

