#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib

def replace_nans(data):
    #load
    try:
        time, norm, *datas = np.loadtxt(data, unpack=True)
        print(np.shape(datas)[0]//2 , ' states detected')
    except:
        print('Error encountered in reading data of', data)
        return
    
    #renaming
    os.system("cp " + data + " " + data[:-4]+"_old.out")
    
    #get header
    with open(data, 'r') as f:
        header = f.readlines()[0:3]
    
    #repace nan's
    norm = np.nan_to_num(norm, nan=1.0)
    datas[0] = np.nan_to_num(datas[0] , nan=1.0)
    for i in range(0,np.shape(datas)[0]):
        datas[i] = np.nan_to_num(datas[i] , nan=0.0)

    #save
    np.savetxt(data, np.transpose([time, norm, *datas]), fmt='%23.12E', delimiter='', header=str(header[0])+str(header[1])+str(header[2][:-1]))

folder = sys.argv[1]
trajs = os.listdir(folder)
print('Looking for data in ', folder)

for traj in trajs:
    if not traj.startswith('TRAJ'):
        continue

    print(traj)
    
    replace_nans(folder + traj+"/output_data/coeff_MCH.out")
    replace_nans(folder + traj+"/output_data/coeff_diag.out")
    replace_nans(folder + traj+"/output_data/coeff_diab.out")
