#!/usr/bin/env python3
import os.path
import sys
import numpy as np
import matplotlib

folder = sys.argv[1]
trajs = os.listdir(folder)

end_struc = []

for traj in trajs:
    if not traj.startswith('TRAJ') or os.path.isfile(folder + traj + '/DONT_ANALYZE'):
        continue

    print(traj)

    with open(folder + traj+'/output.xyz','r') as f:
        im = f.readlines()[-11:]
        im[1] = im[1].replace('\r', '').replace('\n', '') + ' ' + traj +'\n'
        end_struc.append(im)


with open('end_structures.xyz','w') as f:
    for i in end_struc:
        for j in range(len(i)):
            f.write(i[j])

