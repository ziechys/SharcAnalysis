#imports
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import py3Dmol
import pymatgen
from pymatgen.core import Lattice, Structure, Molecule, IStructure, IMolecule
from pymatgen.analysis import molecule_matcher
from pymatgen.core import operations
from pymatgen.io.gaussian import GaussianInput
import copy 
import shutil
from scipy.fft import fft, ifft, fftfreq, rfft, fftshift
import math

def read_trajectory(filename,as_pymatgen=False):
    """
    Read a trajectory file in cartesian coordinates and return string or pymatgen molecule objects for each time step
    """
    #Read in last structure as pymatgen Molecule to extract Natoms
    ref = Molecule.from_file(filename)
    Natoms = len(ref.species)
    Nt     = sum(1 for line in open(filename))//(Natoms+2)
    print(Nt , ' structures found')
    trajs = np.zeros((Nt,3,Natoms))

    symbol = np.loadtxt(filename, skiprows=2,usecols=(0),max_rows=Natoms,unpack=True,dtype=str)
    #read in structures at all times steps
    for i in range(Nt):
        trajs[i] = np.loadtxt(filename, skiprows=(i*11)+2,usecols=(1,2,3),max_rows=Natoms,unpack=True)
    trajs = np.transpose(trajs, (0,2,1))
    if as_pymatgen: #convert to pymatgen molecule object
        trajs_pymatgen = [0]*Nt
        for i in range(Nt):
            trajs_pymatgen[i] = ref.copy()
            for nr in range(Natoms):
                trajs_pymatgen[i][nr].coords = trajs[i][nr]
        return trajs_pymatgen
    else:
        return trajs

def visualize(structure):
    """
    Plot interactively a pymatgen molecule object using py3Dmol
    """
    xyzview = py3Dmol.view(width=400,height=400)
    xyzview.addModel(structure.to('xyz'),'xyz')
    xyzview.setStyle({'stick':{}})
    xyzview.addPropertyLabels("index", {},{'fontColor':'black','font': 'sans-serif', 'fontSize': 28, 'showBackground':False,'alignment':'center'})
    xyzview.show()

def visualize_trajectory(trajectory,speed=10,labels=False):
    """
    Plot interactively and animated a full trajectory consisting of several pymatgen objects using py3Dmol
    """
    traj_string = ''
    for i in trajectory:
        traj_string = traj_string + i.to('xyz') + '\n'
    
    xyzview = py3Dmol.view(width=400,height=400)
    xyzview.addModelsAsFrames(traj_string,'xyz')
    #xyzview.addModel(trajectory[1].to('xyz'),'xyz')
    xyzview.animate({'loop':'forward','reps':1,'interval':speed})
    xyzview.setStyle({'stick':{}})
    if labels:
        xyzview.addPropertyLabels("index", {},{'fontColor':'black','font': 'sans-serif', 'fontSize': 28, 'showBackground':False,'alignment':'center'})
    xyzview.zommTo()
    xyzview.render()
    xyzview.show()

def combine_molecules(mol1,mol2,matching=True):
    """Combines coordinates of two pymatgen molecules"""
    if matching:
        #Perform Hungarian Matching
        matcher = molecule_matcher.HungarianOrderMatcher(mol1)
        mol2 = matcher.fit(mol2)[0]
    new_coords = (mol1.cart_coords + mol2.cart_coords)/2
    new_mol = mol1.copy()
    for nr in range(len(new_mol)):
        new_mol[nr].coords = new_coords[nr]
    return new_mol

def operate_combine_molecules(ref, *fargs_input):
    """Perform symmetry operation and combine coordinates of several pymatgen molecules"""
    
    fargs = copy.deepcopy(fargs_input) #to not alter input lists
    
    #symmetry operation
    for args in fargs: 
        if not args[1] is None:
            args[0].apply_operation(args[1]) 
        
    #Perform Hungarian Matching
    matcher = molecule_matcher.HungarianOrderMatcher(ref)
    for args in fargs:
        args[0] = matcher.fit(args[0])[0] 
        
    #return fargs[0][0]
    
    #Change coordinates
    new_coords = ref.cart_coords
    for args in fargs:
        new_coords = new_coords + args[0].cart_coords
    new_coords = new_coords/(len(fargs)+1)
    new_mol = ref.copy()
    for nr in range(len(new_mol)):
        new_mol[nr].coords = new_coords[nr]
    return new_mol

def operate_combine_trajectories(ref, *fargs_input, weights=[1.0]):
    """Perform symmetry operation and combine coordinates for several trajectories"""
    
    fargs = copy.deepcopy(fargs_input) #to not alter input lists
    
    #weights
    if len(weights) > 1:
        print('Weighting of trajectories was selected')
        if len(weights) != len(fargs)+1:
            print('ERROR: Number of trajectories does not match weighting definition')
            sys.exit(0)
        weights = weights/np.sum(weights)
    else:
        weights = [1/(1+len(fargs))]*(1+len(fargs))
    
    #get shortest trajectory 
    n_steps = len(ref)
    for args in fargs:
        if len(args[0]) < n_steps:
            n_steps = len(args[0])
    
    #symmetry operation
    for args in fargs: 
        if not args[1] is None:
            for molecule in args[0]:
                molecule.apply_operation(args[1])
                
    #loop over time until shortes trajectory ends
    new_traj = [0]*n_steps
    for t in range(n_steps):
        
        #Perform Hungarian Matching
        matcher = molecule_matcher.HungarianOrderMatcher(ref[t])
        for nr_args in range(len(fargs)):
            #for nr_mol in range(len(fargs[nr_args][0])): 
            fargs[nr_args][0][t] = matcher.fit(fargs[nr_args][0][t])[0]
            
        #Combine coordinates to new trajectory file
        new_coords = ref[t].cart_coords*weights[0]
        for nr_args,args in enumerate(fargs):
            new_coords = new_coords + args[0][t].cart_coords*weights[nr_args+1]
        #new_coords = new_coords/(len(fargs)+1)
        new_mol = ref[t].copy()
        for nr in range(len(new_mol)):
            new_mol[nr].coords = new_coords[nr]
        new_traj[t] = new_mol
    
    return new_traj

def traj_to_gaussian(trajectory, folder, title = '', functional = 'PBE1PBE', basis_set = '6-311++g(d,p)', route_parameters={'td': "(singlets,nstates=10)", "nosym":""}, dieze_tag = '#',link0_parameters='',submit_path='',submit_file=False,step=1):
    """
    Convert each time step of a trajectory consisting of pymatgen molecules into a seperate Gaussian input file in separate folder. 
    Submission script can be incorporates as well. 
    """
    if submit_file:
        if os.path.isfile('submit.sh'):
            os.remove('submit.sh')
    for nr in range(0,len(trajectory),step): 
        os.mkdir(folder + '/step_'+ str(nr).zfill(4))
        ffolder = folder + '/' + 'step_'+ str(nr).zfill(4) + '/'
        if title == '':
            ctitle = 'Step ' + str(nr)
        else:
            ctitle = title
        if link0_parameters == '':
            clink0_parameters = {'%chk': 'step_'+ str(nr).zfill(4) + '.chk'}
        else:
            clink0_parameters = link0_parameters
        gau = GaussianInput(trajectory[nr], title = ctitle, functional = functional, basis_set = basis_set, route_parameters=route_parameters, dieze_tag = dieze_tag , link0_parameters=clink0_parameters)
        gau.write_file(filename=ffolder + 'step_'+ str(nr).zfill(4) + '.gjf',cart_coords=True)
       
        #copy submit script in each directory
        if submit_path != '':
            shutil.copy(submit_path , ffolder)
            
        #Create submit.sh file
        if submit_file: 
            with open('submit.sh', 'a') as submitfile:
                submitfile.write('cd ' + ffolder + '\n')
                submitfile.write('sbatch g16_submit ' + 'step_'+ str(nr).zfill(4) + '.gjf'  + '\n')
                submitfile.write('cd - \n')

def readgaussian(file):
    """
    Read in GS and excited state energies from a Gaussian log file. 
    """
    excitations = []
    print('Reading: ', file)
    with open(file, 'r') as datain:
        for line in datain:
            if 'Excited State' in line:
                excitations.append(float(line.split()[4]))
            if 'SCF Done' in line:
                gs = float(line.split()[4])
    return np.array(gs), np.array(excitations)

def PES_1D(GS,Excitations,use_internal_x = True, x_coord={'coord':'bond','index':[0,1]}):
    """
    Read in all Gaussian calculations from one folder and make 1D PES from it. 
    Automatically detect x coordinate based on internal coordinate requested
    """
    #Read in GS, i.e. reference point for PES
    gsenergy, gs_excitations = readgaussian(GS)
    
    #Read in Excitation files
    nsteps  = len(Excitations)
    nstates = len(gs_excitations)
    data = np.zeros((nsteps,nstates))
    energy = np.zeros(nsteps)
    for nr,i in enumerate(Excitations):
        energy[nr] , data[nr] = readgaussian(i)
        
    #Reference energies to GS and in eV
    energy = (energy - gsenergy)*27.2114
    #Shift excited state energies by their GS energy
    for j in range(nsteps):
        data[j] = data[j] + energy[j]
    
    #Get x coordinate spacing by analysing structures as pymatgen Molecule objects
    if use_internal_x:
        if x_coord['coord'] == 'bond':
            print('x coordinate is analysed by bond length between index ', x_coord['index'])
            x = np.zeros(nsteps+1)
            #Get GS (ref) bond length.
            x[0] = Molecule.from_file(GS).get_distance(x_coord['index'][0],x_coord['index'][1])
            for nr,j in enumerate(Excitations):
                x[nr+1] = Molecule.from_file(j).get_distance(x_coord['index'][0],x_coord['index'][1])
        
    state0 = np.concatenate(([0.0],energy))
    states = [0]*nstates
    for i in range(nstates):
        states[i] = np.concatenate(([gs_excitations[i]],np.transpose(data)[i]))
    
    if use_internal_x:
        return x, state0, states
    else:
        return state0, states

def bond_from_traj(trajectory,index1,index2):
    """
    Get specific bond length from trajectory file
    """
    length = np.zeros(len(trajectory))
    for nr,traj in enumerate(trajectory):
        length[nr] = traj.get_distance(index1,index2)
    return length

def dihedral_from_traj(trajectory,index1,index2,index3,index4):
    """
    Get specific bond length from trajectory file
    """
    length = np.zeros(len(trajectory))
    for nr,traj in enumerate(trajectory):
        length[nr] = traj.get_dihedral(index1,index2,index3,index4)
    return length

def traj_to_xyz(trajectory,filename='output.xyz'):
    """
    Write out a trajectory file into a xyz file
    """
    traj_string = ''
    for i in trajectory:
        traj_string = traj_string + i.to('xyz') + '\n'
    
    with open(filename, 'w') as f:
        f.write(traj_string)

def do_FT(time, data,plottime=False,damping=False,padding=0.0,mean=False,corrdamp=False):
    """
    FT anaylsis of data. 
    Possible FT options are: splines, damping, padding, mean value substraction, post-damping correction.
    """
    #Some useful conversion
    cm2au = 4.5563352529120*10**(-6)
    au2fs = 0.024188843265857

    #Process data: Mean offset and damping
    if mean:
        offset = np.mean(data[:   round(math.floor(len(data)/np.pi/2) * np.pi * 2 )   ])
        #offset = np.mean(data)
        print('Data offset found of ', offset)
        data = data - offset
    if damping:
        #Damping function
        damping_func = np.array([np.cos( ( i /(time[-1]) ) *  (np.pi/2) )**2 for i in time  ])
        data = data * damping_func

    if plottime:
        #Plot input data in time
        plt.plot(time*au2fs, data,label=legend)
        plt.xlabel('t / fs')
        plt.ylabel('E / eV')
        plt.show()

    #FFT to energy/freq
    if padding > 0.0:
        gridpoints = padding
    else:
        gridpoints = len(data)
    dt = abs(time[1] - time[0])
    dE = (2*np.pi)/(dt)
    FFT_energy = rfft(data,n=gridpoints)
    energy = np.linspace(0, dE/2, num=gridpoints//2+1, endpoint=True)
    #Plot FFT results
    if corrdamp:
        FFT_damp = rfft(damping_func,n=gridpoints)
        energy_out = energy/cm2au
        data_out   = abs(FFT_energy/FFT_damp)/np.max( abs(FFT_energy/FFT_damp))
    else:
        energy_out = energy[1:]/cm2au #remove first point cause heavily influenced by damping
        data_out   = abs(FFT_energy[1:])/np.max( abs(FFT_energy[1:]))
    
    data_out = data_out / (np.sum(data_out)*dE)
    #data_out = data_out/np.max(data_out)

    return np.array(energy_out),np.array(data_out)

def ft_bond_from_traj(trajectory,time,index1,index2):
    au2fs = 0.024188843265857
    return do_FT(time/au2fs,bond_from_traj(trajectory,index1,index2),mean=True,padding=10000,damping=True,corrdamp=True)


def ft_dihedral_from_traj(trajectory,time,index1,index2,index3,index4):
    au2fs = 0.024188843265857
    return do_FT(time/au2fs,dihedral_from_traj(trajectory,index1,index2,index3,index4),mean=True,padding=10000,damping=True,corrdamp=True)
