import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase.io import read
from ase.visualize import view
import os
import shutil
from ase.io.vasp import write_vasp
import sys
import random

rotate_no = 10   #specify the number of rotated structures you need
path = os.getcwd()
path1 = r'/ihome/kjohnson/laa147/DP_Graphamine/perturb/rotate-perturb'

# #################  Obtaining atoms information from the original POSCAR file  #######################

gamine = read('POSCAR') 
def get_atom_index(atoms, atom_char):
    '''
    Takes in Atoms object and returns the indices of where a particular atom is
    '''
    Element_idx_bool = atoms.symbols==atom_char
    Element_idx = [i for i, x in enumerate(Element_idx_bool) if x]
    return Element_idx

# Finding the indices for C, N, H
H_idx = get_atom_index(gamine, 'H')
N_idx = get_atom_index(gamine, 'N')
C_idx = get_atom_index(gamine, 'C')

pos = gamine.get_positions()

C_zsurf = np.mean(pos[C_idx][:,2])


H_all= pos[H_idx]
H_top = H_all[np.argwhere([H_all[:,-1]>C_zsurf])[:,1]]
H_top_idx = [0 for element in range(len(H_top))]

z=0
for i in range(len(H_top)):
    var = H_top[i]
    for j in range(len(H_all)):
        var1 = pos[H_idx][j]
        if (var1 == var).all():
            H_top_idx[z] = H_idx[j]
            z+=1

# #################  Obtaining atoms distances from the original POSCAR file  #######################
dis_list = gamine.get_all_distances(mic=True, vector=True)
N_all = pos[N_idx]
total_len = len(H_top) + len(N_all)
flist = [[0]*len(pos)]*total_len
flist_index = np.concatenate((H_top_idx, N_idx))

for i in range(total_len):
    for j in range(len(pos)):
        k = flist_index[i] 
        flist[i][j] = dis_list[k][j]




distance = [0 for element in range(len(H_idx))]
temp = [0 for element in range(len(H_idx))]
Neighbor_list = [[0 for i in range(3)] for j in range(len(N_all))]
A = [0 for element in range(len(N_all))]
B = [0 for element in range(len(N_all))]

# #################  Obtaining the nearest neighbors of NH2 groups  #######################
for i in range(len(N_all)):
    j = N_idx[i]
    n =0
    for k in range(len(H_idx)):
        l = H_idx[k]
        x = dis_list[j][l][0] - dis_list[j][j][0]
        y = dis_list[j][l][1] - dis_list[j][j][1]
        z = dis_list[j][l][2] - dis_list[j][j][2]
        distance[k] = np.sqrt(x**2 + y**2 + z**2)
    A[i], B[i] = np.partition(distance,1)[0:2]
    Neighbor_list[i][n] = N_idx[i]
    for m in range(len(H_idx)):
        if distance[m] == A[i]:
            n+=1
            Neighbor_list[i][n] = H_idx[m]
        if distance[m] == B[i]:
            n+=1
            Neighbor_list[i][n] = H_idx[m]

# ################# Moving the H atoms to follow the minimum image convention  #######################
poscar1=np.loadtxt('POSCAR',skiprows=8)

for p in range(len(N_all)):
    idx1 = Neighbor_list[p][0] 
    idx2 = Neighbor_list[p][1]  
    idx3 = Neighbor_list[p][2] 
    var1 = poscar1[idx2][1] - poscar1[idx1][1]
    var2 = poscar1[idx3][1] - poscar1[idx1][1]
    if var1 < -0.5:
        poscar1[idx2][1] +=1
    if var2 < -0.5:
        poscar1[idx3][1] +=1
    if var1 > 0.5:
        poscar1[idx2][1] -=1
    if var2 > 0.5:
        poscar1[idx3][1] -=1
np.savetxt('POSCAR-temp', poscar1)    

def copy_line_by_line(input_file1, output_file1, start_line1, total_no_lines1):
    # Iterate over each line in the file
    for line in input_file1.readlines()[start_line1:(start_line1 + total_no_lines1)]:
        # Write to the file if conditions are met
        output_file1.write(str(line))
        

    # Release used resources
    input_file1.close()
    output_file1.close()


input_file1 = open("POSCAR", "r")
output_file1 = open("POSCAR-temp2", "a")


start_line1 = 0
total_lines1 = 8

copy_line_by_line(input_file1, output_file1, start_line1, total_lines1)

def copy_(input_file2, output_file2, start_line2, total_no_lines2):
    
    # Iterate over each line in the file
    for line in input_file2.readlines()[start_line2:(start_line2 + total_no_lines2)]:
        # Write to the file if conditions are met
        output_file2.write(str(line))
        

    # Release used resources
    input_file2.close()
    output_file2.close()
   


input_file2 = open("POSCAR-temp", "r")
output_file2 = open("POSCAR-temp2", "a")


start_line2 = 0
total_lines2 = len(gamine.get_positions())

copy_(input_file2, output_file2, start_line2, total_lines2)

value = 0

for file in os.listdir(path):
    name = os.path.basename('POSCAR')
    
    for i in range(value, value+rotate_no):
        newname = f"{name}-{i}"
        shutil.copyfile('POSCAR-temp2', newname)
    
    

# #################  Rotating the H-atoms to create new rotated POSCAR files  #######################

def moveH(H_atom1,H_atom2, N_center,angle):
    x1=N_center[0]-H_atom1[0]
    y1=N_center[1]-H_atom1[1]
    R1=np.sqrt(x1**2 + y1**2)
    x2=N_center[0]-H_atom2[0]
    y2=N_center[1]-H_atom2[1]
    R2=np.sqrt(x2**2 + y2**2)
    theta=random.uniform(0,1)*2*np.pi
    theta2 = theta + ((angle/180)*np.pi)
    x1_new=N_center[0]+R1*np.cos(theta)
    y1_new=N_center[1]+R1*np.sin(theta)
    new_H1_pos=np.array([x1_new, y1_new, H_atom1[2]])
    x2_new=N_center[0]+R2*np.cos(theta2)
    y2_new=N_center[1]+R2*np.sin(theta2)
    new_H2_pos=np.array([x2_new, y2_new, H_atom2[2]])
    return new_H1_pos, new_H2_pos

for m in range(rotate_no):
    obj=read(f'POSCAR-{m}')
    for n in range(len(N_all)):
        idx1 = Neighbor_list[n][0]
        idx2 = Neighbor_list[n][1]  
        idx3 = Neighbor_list[n][2] 
        old_H1 = obj.get_positions()[idx2]
        old_H2 = obj.get_positions()[idx3]
        N_cent = obj.get_positions()[idx1]
        angle_b = obj.get_angle(idx2,idx1,idx3)
        new_H1, new_H2 = moveH(old_H1, old_H2, N_cent,angle_b)
        p1 = idx2
        p2 = idx3
        obj.positions[p1] = new_H1
        obj.positions[p2] = new_H2
    write(f'POSCAR-{m}', obj, format='vasp')
    
# for o in range(rotate_no):
#     atoms = read(f'POSCAR-{o}')
#     write_vasp(f'POSCAR-{o}', atoms, direct = True, vasp5 = True)
               
print("Created new POSCAR files with random rotations")
