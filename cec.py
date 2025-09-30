import numpy as np
from ase.io import read, write
from ase import neighborlist
import matplotlib.pyplot as plt
import os
import sys


class CEC():
    def __init__(self, working_dir, dumptraj, dumptraj_tag, noAtoms, xdatcar_prepend, timestep=0.25, rcut=1.3):
        '''
        1.
        2.
        3.
        4.
        '''

        self.working_dir = working_dir
        self.dumptraj = dumptraj
        self.dumptraj_tag = dumptraj_tag #"LAMMPS" Or "XDATCAR"
        self.noAtoms = noAtoms
        self.ts = timestep
        self.rcut = rcut
        self.xdatcar_path = f'{self.working_dir}/XDATCAR_CEC'
        self.xdatcar_prepend = xdatcar_prepend

    def get_atom_index(self, atoms, atom_char):
        '''
        Takes in Atoms object and returns the indices of where a particular atom is
        '''
        Element_idx_bool = atoms.symbols==atom_char
        Element_idx = [i for i, x in enumerate(Element_idx_bool) if x]
        return Element_idx

    def xdatcar2array(self):
        '''
        Uses the XDATCAR and stores the coordinates into an array after sorting based on the index.
        Add this to the main CEC code and add a condition that asks for XDATCAR or lammpstrj
        '''
        dump = open(self.dumptraj)
        count = 0
        direct_count = 0
        start = 0
        tag = 0
        frame = []
        atm = []
        for l in dump:

            if tag == 1:
                coords = [float(l.split()[0]),float(l.split()[1]),float(l.split()[2])]
                atm.append(coords)
                count += 1

            if count == self.noAtoms:
                tag = 0
                frame.append(atm)
                atm = []
                count = 0
            if l.split()[0] == 'Direct':
                direct_count += 1
                tag = 1
        self.frame = np.array(frame)
        
    def dump2array(self):
        '''
        Uses the dump traj file from lammps and stores the coordinates into an array after sorting based on the index.
        '''
        dump = open(self.dumptraj)
        count = 0
        direct_count = 0
        start = 0
        tag = 0
        frame = []
        atm = []
        for l in dump:
            if l.split()[0] == 'ITEM:':
                if l.split()[1] == 'TIMESTEP':
                    count = 0
            if tag == 1:
                atm.append([int(l.split()[0]),float(l.split()[2]),float(l.split()[3]),float(l.split()[4])])
                count += 1
            if count == self.noAtoms:
                tag = 0
                frame.append(atm)
                atm = []
            if l.split()[0] == 'ITEM:':
                if l.split()[1] == 'ATOMS':
                    direct_count += 1
                    tag = 1
                    count = 0
        frame = np.array(frame)
        frame_sort = []
        for i in range (len(frame)):
            fr_i = frame[i]
            fr_is = fr_i[fr_i[:,0].argsort()]
            frame_sort.append(fr_is)
        frame_sort = np.array(frame_sort)
        self.frame = frame_sort

    def read_first_atoms_object(self):
        self.atomsObj = read(self.dumptraj, index=0)
        return self.atomsObj
    
    def convert_to_adjacency(matrix):
        start = 0
        res = []
        lst = []
        n = len(matrix)

        for i in range(n):
            res.append(lst*n)
        while start < n:
            y = matrix[start]
            for i in range(len(y)):
                if y[i] == True:
                    res[start].append(i)
            start += 1
        return res
    
    

    def CECmain(self,write_cec_path, write_xdatcar=True, id0_path=None):
        '''
        Uses the dump traj file from lammps and stores the coordinates into an array after sorting based on the index.
        '''
        if write_xdatcar:
            print('DUMP to CEC in XDATCAR')
            print('------------------------')
            self.xdatcar_path = f'{self.working_dir}/XDATCAR_CEC'
            print(id0_path)
            xdatcar = open(self.xdatcar_path,'w+')
            
        if id0_path is None:
            print('READING AND WRITING FIRST ATOMS OBJECT')
            print('------------------------')
            id0 = self.read_first_atoms_object()
            # only run this for the run. 
            N_idx_bool = id0.symbols=='He'
            N_idx = [i for i, x in enumerate(N_idx_bool) if x]
            id0.symbols[N_idx] = 'N'

            C_idx_bool = id0.symbols=='H'
            C_idx = [i for i, x in enumerate(C_idx_bool) if x]
            id0.symbols[C_idx] = 'C'

            H_idx_bool = id0.symbols=='Li'
            H_idx = [i for i, x in enumerate(H_idx_bool) if x]
            id0.symbols[H_idx] = 'H'

            id0_path = f'{self.working_dir}/POSCAR'
            print('id0_path = ', id0_path)
            write(id0_path, id0) # Assumes POSCAR
            self.id0_path = id0_path # Save to use for next iterations
        else:
            id0 = read(id0_path)
            self.id0_path = id0_path
            print('FOUND BASE ATOMS OBJECT')
            
        H_idx = self.get_atom_index(id0, 'H')
        N_idx = self.get_atom_index(id0, 'N')
        
        # Create a neighborlist for all N atoms
        conn_bool_N = id0.get_all_distances(mic=True)[N_idx][:,N_idx] < 3 # assuming N-N distance is just under 3 Angs
        conn_bool_N
        N_adj_list = {}
        conn_all=conn_bool_N[0]*N_idx
        conn_fin=conn_all[conn_all!=0]

        for i_n,n in enumerate(N_idx):
            conn_all=conn_bool_N[i_n]*N_idx
            conn_fin=conn_all[conn_all!=0]
            N_adj_list[n] = conn_fin
        
        if self.dumptraj_tag == "LAMMPS":
            self.dump2array() # generates self.frame
        elif self.dumptraj_tag == "XDATCAR":
            self.xdatcar2array()
        
        # perform a global scan to generate the first CEC. Since the proton is manually added, 
        where_CEC=0 
        CEC_N_idx = N_idx[where_CEC]
        prevCEC = CEC_N_idx
        CEC_track = []
        cec_list = []
        pot_list = []
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        # TODO: Add in a constraint where the next CEC is less than X distance away. X is calculated as the 'min' possible N-N distance in GNH2. 
        # If there are instances with more than 1 potential CECs then save those instances (just the index) 

        for i,f in enumerate(self.frame):
            if len(f.T) == 4:
                pos_f = f[:,1:]
            elif len(f.T) == 3:
                pos_f = f

            id0.set_scaled_positions(pos_f)
            
            if i == 0:
                # Finds the N that has the two closest H atoms. 
                conn_NH = id0.get_all_distances(mic=True)[N_idx][:,H_idx]
                sorted_conn_NH = np.sort(conn_NH)
                where_CEC = np.argmin(np.sum(sorted_conn_NH[:,0:3], axis=1))
                CEC_N_idx = N_idx[where_CEC]
                count4 +=1
                
            else:
                pot_CEC = []
                for n in N_adj_list[prevCEC]:
                    conn_bool = id0.get_distances(n, H_idx, mic=True) < self.rcut
                    conn_H = np.sum(conn_bool) 
                    if conn_H >= 3:
                        pot_CEC.append(n)
                        count3 +=1
                           
    
                if len(pot_CEC)>=2:
                    # iterate over all potential CECs
                    min_dist = np.inf
                    for pc in pot_CEC:
                        N_N_dist = id0.get_distance(prevCEC, pc,mic=True) #added mic=True - laa
                        if N_N_dist < min_dist:
                            min_dist = N_N_dist
                            CEC_N_idx = pc
                            count2 += 1
   
                            
                            
                elif len(pot_CEC) == 1:
                    if pot_CEC[0] == 0:
                        CEC_N_idx = prevCEC
                    else:
                        CEC_N_idx = pot_CEC[0]

                    
            scaled_pos = id0.get_scaled_positions()[CEC_N_idx]
            prevCEC = CEC_N_idx
            cec_list.append(prevCEC)
            if write_xdatcar:
                xdatcar.write(f"{CEC_N_idx}\n")    #added by LAA147 to print out the N atomic index to determine proton transport mechanism for graphamine
                xdatcar.write(f'Direct configuration= {i+1}\n')
                xdatcar.write(f' {scaled_pos[0]:1.4f} {scaled_pos[1]:1.4f} {scaled_pos[2]:1.4f}\n')
        xdatcar.close()        
        with open(self.xdatcar_path, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(self.xdatcar_prepend.rstrip('\r\n') + '\n' + content)
        self.cec_list = np.array(cec_list)
        np.savetxt(f'{write_cec_path}', self.cec_list)

    def xdatcar2unwrapped(self):
        print('XDATCAR to UNWRAPPED xyz')
        print('------------------------')
        xdatcar2unwrappedxyz_path = '/ihome/kjohnson/laa147/DP_Graphamine/CEC/source_code/xdatcar2unwrappedxyz.py'
        os.chdir(self.working_dir)
        os.system(f'python {xdatcar2unwrappedxyz_path} -i {self.xdatcar_path}')

    def calc_msd(self):
        print('UNWRAPPED xyz to MSD')
        print('------------------------')
        calc_msd_path = '/ihome/kjohnson/laa147/DP_Graphamine/CEC/source_code/xyztomsd_un.py'
        os.chdir(self.working_dir)
        os.system(f'python {calc_msd_path} -i coordunwrapped_CEC.xyz')
