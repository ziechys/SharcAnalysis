#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#parser: folder nstates Geo.inp

folder = sys.argv[1]
nstates = int(sys.argv[2])
geo = sys.argv[3]
labels = sys.argv[4].split(',')
trajs = os.listdir(folder)
curr = os.getcwd()

print('Folder of trajectories: ', folder)
print('Current trajectory: ', curr)
print('Number of states: ', nstates)
print('Geometry analysis file: ', geo)

print('')

try:
    os.mkdir(curr + '/' + "graphs")
except Exception:
    pass

for traj in trajs:
    if not os.path.isfile(folder + traj + "/output.dat"):
        continue
    print(traj)
    os.chdir(folder + traj)

    os.system("$SHARC/data_extractor.x output.dat")
    os.system("$SHARC/make_gnuscript.py "+ str(nstates)  +" > plot.gp")
    
    with open('plot2.gp','w') as dataout:
        dataout.write('set terminal pdf\n')
        dataout.write("set output 'energies.pdf'\n")
        dataout.write("set pointsize 0.2\n")
    
        name = ["'Diag_energies.pdf'","'MCH_amplitudes.pdf'", "'Diag_amplitudes.pdf'", "'hop_prop.pdf'"]
        index = 0
    
        with open('plot.gp', 'r') as datain:
            for line in datain:
                if "pause" in line and index <4:
                    dataout.write('set output ' + name[index] )
                    index = index + 1
                elif  "pause" in line and index ==4:
                    break
                else:
                    dataout.write(line)
    os.system("gnuplot plot2.gp")       
    
    os.system("cp "+curr+'/'+geo+" ./")
    #with open('Geo.inp','w') as dataout:
    #    dataout.write("r 5 1\n")
    #    dataout.write("r 1 2\n")
    #    dataout.write("r 2 3\n")
    #    dataout.write("r 3 4\n")
    #    dataout.write("r 4 5")
    
    os.system("$SHARC/geo.py -g output.xyz -t 0.5 < Geo.inp > Geo.out")

    #with open('Geo.gp','w') as dataout:
    #    dataout.write("set terminal pdf\nset output 'Geo.pdf'\nset xlabel 't / fs'\nset ylabel 'r / A'\n")
    #    dataout.write("p 'Geo.out' w l t 'O5-C1', '' u 1:3 w l t 'C1-C2', '' u 1:4 w l t 'C2-C3', '' u 1:5 w l t 'C3-C4', '' u 1:6 w l t 'C4-O5'")
    
    #os.system("gnuplot Geo.gp")
    time, *geomdata = np.loadtxt('Geo.out', unpack=True)
    for nr,i in enumerate(geomdata):
        plt.plot(time, i, label = labels[nr])
    plt.xlabel('t / fs')
    plt.legend(loc='upper left')
    plt.savefig('Geo.pdf')
    plt.close()
    #plt.show()
    #sys.exit(0)

    os.system("cp energies.pdf "+curr+"/graphs/" + traj + "_energies.pdf")
    os.system("cp Diag_energies.pdf "+curr+"/graphs/" + traj + "_Diag_energies.pdf")
    os.system("cp Geo.pdf "+curr+"/graphs/" + traj + "_Geo.pdf")
    os.system("cp hop_prop.pdf "+curr+"/graphs/" + traj + "_hop_prop.pdf")    
    os.system("cp MCH_amplitudes.pdf "+curr+"/graphs/" + traj + "_MCH_amplitudes.pdf")
    os.system("cp Diag_amplitudes.pdf "+curr+"/graphs/" + traj + "_Diag_amplitudes.pdf")

    os.chdir(curr)

