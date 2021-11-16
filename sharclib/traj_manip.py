"""
author: Felix Plasser
version: 1.0
description: package for rotating and superimposing NewtonX trajectories
"""

import os, sys
try:
    import numpy
except ImportError:
    print('numpy package not installed')
    sys.exit()
    
try:
    import openbabel
except ImportError:
    print('openbabel.py package not installed')
    sys.exit()

try:
    import file_handler, vib_molden, struc_linalg
except ImportError:
    print('file_handler, vib_molden or struc_linalg not found. They should be part of this package. Check the installation or change the PYTHONPATH environment variable.')

class trajectory:
    """
    A NewtonX trajectory. Data is taken from the dyn.mld file.
    """

    def __init__(self, path, ref_struc=None, dt=0.5,GSstop=False):
        self.path = path
        self.dt = dt # timestep length in fs
        self.ref_struc = ref_struc # structure that the other structures are superimposed on
        self.hop = None
        
        self.structures = []
        t = 0.
        
        ## read in the molecules with openbabel
        obconversion = openbabel.OBConversion()
        obconversion.SetInFormat('xyz')
        mol = openbabel.OBMol()        
        
        #print 'timestep:'
        # the first structure is read in
        # return a boolean: True if read in successfull  
        notatend = obconversion.ReadFile(mol, path+'/output.xyz')
        self.num_at = mol.NumAtoms()
        #obconversion.WriteFile(mol,'teststructure.xyz')
        #sys.exit(0)               

        if GSstop:
            self.get_GShop()
            print('GS hop found at step: ',self.hop)
 
        if GSstop and not self.hop is None:
            while notatend and t < self.hop * self.dt:
                self.structures += [struc_linalg.structure(str(t) + ' fs')] # define the structure object with name based on time
                self.structures[-1].get_mol(mol, path) # read in the data from ObMol object
                t = t + self.dt

                mol = openbabel.OBMol()               #new ObMol object
                notatend = obconversion.Read(mol)     #next structure assigned to mol

        else:
            # the structures are read in as struc_linalg objects as long as structures are present
            while notatend:
                self.structures += [struc_linalg.structure(str(t) + ' fs')] # define the structure object with name based on time
                self.structures[-1].get_mol(mol, path) # read in the data from ObMol object
                t = t + self.dt
                
                mol = openbabel.OBMol()               #new ObMol object
                notatend = obconversion.Read(mol)	  #next structure assigned to mol
            
            
        self.num_tsteps = len(self.structures)
        #del obconversion

    def write_file(self,nr,fileout='.xyz'):
        """
        Write xyz file of trajectory number nr
        """
        obconversion = openbabel.OBConversion()
        obconversion.SetInFormat('xyz')

        obconversion.WriteFile(self.structures[nr].mol,self.path+'/'+str(nr)+fileout)


    def get_GShop(self):
        """
        get time / step of GS hop in a forced to GS hop trajectory
        """
        with open(self.path+"/output.lis", 'r') as f:
            for nr, line in enumerate(f):
                if 'Forced jump to ground state' in line:
                    self.hop = int(f.readlines()[1].split()[0])
                    break

            

    def clean_structures(self, mass_wt_pw=1):
        """
        Clean the structures.
        Currently just superposition, sub can be overwritten if more has to be done.
        """
        self.superimpose_structures(mass_wt_pw=mass_wt_pw) 

    def superimpose_structures(self, mass_wt_pw=1):
        """
        Superposition of the structures using a quaternion fit.
        """
        ref_structures = [] # aligned structures
        for structure in self.structures:
            ref_structures += [structure.ret_superimposed_structure(self.ref_struc, mass_wt_pw=mass_wt_pw)]
        self.structures = ref_structures

    def print_file(self, out_path, out_file='dyn.mld'):
        """
        Print the output file with the aligned structures.
        """
        try:
            os.makedirs(os.path.join(out_path, 'RESULTS'))
        except OSError:
            pass
        
        obconversion = openbabel.OBConversion()
        obconversion.SetOutFormat('xyz')
        
        obconversion.WriteFile(self.structures[0].mol, os.path.join(out_path, 'RESULTS', out_file))
        
        for structure in self.structures[1:]:
            obconversion.Write(structure.mol)

        #obconversion.CloseOutFile()
        
    def ret_coor_matrix(self, relative=False,filltoend=None):
        """
        Returns a 3N x T matrix with all the coordinates of the timesteps.
        If <relative=True> coordinates relative to <ref_struc> are returned.
        """
        if relative:
            ref_vect = self.ref_struc.ret_vector()
            coor_list = [structure.ret_vector() - ref_vect for structure in self.structures]
        else:
            coor_list = [structure.ret_vector() for structure in self.structures]
      
        if not filltoend is None:
            coor_list = coor_list + [self.structures[-1].ret_vector()] * filltoend

        return numpy.array(coor_list)
    
    def autocorr(self, mass_mat=None, out_file='RESULTS/autocorr.txt'):
        """
        Compute the autocorrelation function.
        """
        # +++ 
        tm = file_handler.table_maker(2 * [20])
        
        ac_list = []
        for disp in range(self.num_tsteps):
            print(disp)
            ac_list += [self.ret_autocorr_disp(disp=disp, mass_mat=mass_mat)]
            
        for disp in range(self.num_tsteps):
            tm.write_line([float(disp)*self.dt, ac_list[disp]])
            
        tm.write_to_file(out_file)
    
    def ret_autocorr_disp(self, disp, mass_mat=None):
        """
        Return the value of the autocorrelation function at <disp>.
        """
        coor_mat = self.ret_coor_matrix(relative=True)
        
        if mass_mat == None:
            mass_mat = numpy.identity(3 * self.num_at)
            
        m_coor_mat = numpy.dot(coor_mat, mass_mat)
        
        ret_num = 0.
        #print len(m_coor_mat)
        for ind, vec in enumerate(m_coor_mat):
            #print ind, (ind+disp)%self.num_tsteps
            ret_num += numpy.dot(vec, coor_mat[(ind+disp)%self.num_tsteps])
            
        ret_num = ret_num / self.num_tsteps
            
        return ret_num

    def normal_mode_analysis(self, nma_mat, def_struc, header=None, out_file='nma.txt', abs_list=[],timestep=1.0):
        """
        Perform a normal mode analysis and print the result to <out_file>.
        Normal modes are specified by numpy.array <nma_mat> (in cartesian coordinates), <nma_mat> is the inverse of the normal mode matrix.
        The <def_struc> should be the same structure that trajectories were superimposed onto.
        A matrix with all the nma vectors is returned.
        <abs_list> constains a list of normal mode numbers where the absolute value is taken (the first mode's index is 1)
        """
        # +++ include weighting for every mode
            # for example according to zero point vibrations
        # abs_list would better only be considered with the averaging

        def_vect = def_struc.ret_vector()
        
        tm = file_handler.table_maker([35]+len(nma_mat[0])*[20])
        if header == None: header = [i+1 for i in range(len(nma_mat[0]))]  #header as numbered modes 
        tm.write_header_line(header[0])  # column number
        tm.write_header_line(header[1])  # time, frequencies
        tm.write_header_line(header[2])  # blank, T = 1/freq
        nma_list = []
        for istruct,structure in enumerate(self.structures):
            diff = structure.ret_vector() - def_vect
            nma_vect = numpy.dot(diff, nma_mat) # one time step in the normal mode basis
            #Absolute values or not ?
            #for i in abs_list:
            #    nma_vect[i-1] = abs(nma_vect[i-1])
            nma_vect = abs(nma_vect)
            tm.write_line([timestep*istruct]+[coor for coor in nma_vect])
            nma_list += [nma_vect]
        tm.write_to_file(out_file) #Normal mode decomposition for each traj for each time step
        
        return numpy.array(nma_list)
            
if __name__=='__main__':
    ref_struc = struc_linalg.structure('ref_struc') # define the structure that all the time step structures are superimposed onto
    ref_struc.read_file('/home2/plasserf/calcs/BIP/opt/es/cc2/SV_P/DK/coord', 'tmol')
    
    traj = trajectory(path='.',ref_struc=ref_struc)
    traj.autocorr()
