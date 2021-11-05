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

    if not os.path.isfile(folder + traj + "/output.dat"):
        continue

    print(traj)

    #with open(folder + traj+"/output_data/coeff_MCH.out", "r") as myfile:
    #    lineList = myfile.readlines()
    #    print( lineList[len(lineList)-1] )

#    sys.exit(0)

    datafile = folder + traj+"/output_data/coeff_MCH.out"
    datafile2 = folder + traj+"/output_data/coeff_MCH_new.out"

    #os.system('head -n-1 ' + datafile + ' > ' + datafile2)
    #os.system('mv '+ datafile2 + ' ' + datafile)
    #sys.exit(0)


    with open(datafile, "r") as myfile:
        lineList = myfile.readlines()
        #print( lineList[len(lineList)-1] )
        #print( lineList[len(lineList)-1].split() )
        lastline =  lineList[len(lineList)-1].split()[0] + ' 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    
    os.system('head -n-1 ' + datafile + ' > ' + datafile2)
    os.system('mv '+ datafile2 + ' ' + datafile)

    with open(datafile, "a") as myfile:
        myfile.write(lastline)  

    #sys.exit(0)
#replace_nans(folder + traj+"/output_data/coeff_MCH.out")
    #replace_nans(folder + traj+"/output_data/coeff_diag.out")
    #replace_nans(folder + traj+"/output_data/coeff_diab.out")
