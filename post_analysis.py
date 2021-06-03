#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re
import pandas as pd

def get_grad(data,nr_traj):
	#load
	try:
		step, time, diag, state, kin, pot, energy, gradient, dm, s, runtime = np.loadtxt(data+"/output.lis", unpack=True)	
	except ValueError:
		print(data[:-11] +  '   Error in dataset')
		return

	#Get hop step
	hop = 'none'
	with open(data+"/output.lis", 'r') as f:
		for nr, line in enumerate(f):
			if 'Forced jump to ground state' in line:
				hop = f.readlines()[1].split()[0]
				break

	#Check for after hop hop
	time2, cof, s1, *others = np.loadtxt(data+"/output_data/coeff_class_MCH.out", unpack=True)
	if not hop == 'none':
		switch = np.all(s1[int(hop):])
	else:
		switch = ''
		nones.append(data+"/output.lis"[:-11])

	#Check for doubled gradients
	prev = gradient[0]
	for i in range(1,len(gradient)):
		if gradient[i] == prev:
        #print('Equal gradient of %8' + str(gradient[i]) + ' found at step ' + str(step[i]) + ', i.e. ' + str(time[i]) + ' fs')
        #print("%8" % str(gradient[i]) + '   ' +  str(step[i]) +  '   ' +str(time[i]) + ' fs')
			if not hop == 'none':
				if float(int(hop)/2) <  errordata[nr_traj][1]:
					usable.append(data+"/output.lis"[:-11])
				print(data+"/output.lis"[:-11] + '   '  + "%8f" % gradient[i] + '   '  + '%4i' % step[i] + '   '  + '%10f' % time[i]  + ' fs   ' + '%7.2f' % float(int(hop)/2) + '      ' + str(switch).rjust(5)  + '              ' + '%7.2f' % errordata[nr_traj][1] )
				return
			else:
				print(data+"/output.lis"[:-11] + '   '  + "%8f" % gradient[i] + '   '  + '%4i' % step[i] + '   '  + '%10f' % time[i]  + ' fs   ' + hop + '      ' + str(switch).rjust(5) + '              ' + '%7.2f' % errordata[nr_traj][1])
				return
		prev = gradient[i]

	if not hop == 'none':
		if float(int(hop)/2) <  errordata[nr_traj][1]:
                    usable.append(data+"/output.lis"[:-11])
		print(data+"/output.lis"[:-11] +  '   No doubled gradients              '  + '%7.2f' %  float(int(hop)/2) + '      ' + str(switch).rjust(5) + '              ' + '%7.2f' % errordata[nr_traj][1])
	else:
		print(data+"/output.lis"[:-11] +  '   No doubled gradients              '  + hop.rjust(7) + '      ' + str(switch).rjust(5) + '              ' + '%7.2f' % errordata[nr_traj][1])

traj = []
thresh = []

start = False
with open(sys.argv[2], 'r') as datain:
    for line in datain:
        if 'Trajectory Files?' in line:
            start = True
        if 'TRAJ' in line and start:
            #print(line)
            datasearch = re.search('TRAJ_(\d+)\D+\d+.\d+\s+(\d+[.,]?\d+)',line)
            #print(datasearch)
            traj.append(int(datasearch.group(1)))
            thresh.append(float(datasearch.group(2)))

data = []
data.append(traj)
data.append(thresh)

data = np.transpose(data)
errordata = data[data[:, 0].argsort()]

usable = []
nones = []

folder = sys.argv[1]
trajs = os.listdir(folder)

print('Trajectory | gradient | step | time       |   GS hop   | no switch after GShop | T_use')
print()

traj_nr = 0

for traj in trajs:
	if not traj.startswith('TRAJ'):
		continue

	#print(traj)
	get_grad(folder+traj,traj_nr)
	traj_nr = traj_nr + 1

print(str(len(usable)) +' usable Trajectories:')
for i in usable:
	print(i)
	os.system('rm ' + i + '/DONT_ANALYZE')
print(str(len(nones)) + ' trajectories without GShop:')
for i in nones:
	print(i)
	os.system('rm ' + i + '/DONT_ANALYZE')
