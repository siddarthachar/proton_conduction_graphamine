import numpy as np
import argparse
from ase.io import read
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Name of XDATCAR file to be read", default="XDATCAR")
args = parser.parse_args()
print( "Opening file {}".format(args.input))

xdatcar = open(args.input, 'r')
outputunwrapped = open('coordunwrapped_CEC.xyz', 'w')
poscar = read('/ihome/kjohnson/laa147/DP_Graphamine/CEC/POSCAR_prod')
cart_positions = poscar.get_scaled_positions()

# Read in the system name:
system = xdatcar.readline()
# Read in the scale factor:
scale = float(xdatcar.readline().rstrip('\n'))
# Read the unit cell vectors, aa, bb, cc:
vec = xdatcar.readline().rstrip('\n').split()
aa = np.array([float(vec[0]), float(vec[1]),float(vec[2])])
vec = xdatcar.readline().rstrip('\n').split()
bb = np.array([float(vec[0]), float(vec[1]),float(vec[2])])
vec = xdatcar.readline().rstrip('\n').split()
cc = np.array([float(vec[0]), float(vec[1]),float(vec[2])])
unitcell = np.array([aa,bb,cc])
# Read in the element names from the XDATCAR file
element_names = xdatcar.readline().rstrip('\n').split()

element_dict = {}
# Read the number of atoms for each element type:
line = xdatcar.readline().rstrip('\n').split()
i = 0
Natoms = 0
for el in element_names:
    element_dict[el] = int(line[i])
    Natoms += int(line[i])
    print(Natoms)
    i += 1

prev_step = np.zeros((Natoms,3))
curr_step = np.zeros((Natoms,3))
first_step = np.zeros((Natoms,3))
unwrapped = np.zeros((Natoms,3))
old_step = np.zeros((Natoms,3))

print("New file writing now")
    
Ncoord_list = [0]*12
for l in range(12):
    m = 24+l
    Ncoord_list[l] += cart_positions[m]
    
Nstep = 0
while True:
    Nstep += 1
    if Nstep == 1:
# Read in the first configuration to start with:
        line = xdatcar.readline()
        if len(line) == 0:
            break
        outputunwrapped.write(str(Natoms)+"\n"+line)
        for i in range(Natoms):
            coord0 = xdatcar.readline().rstrip('\n').split() #modified by laa147 to account for the nitrogen atom index in the coordunwrapped.xyz
            coords0 = xdatcar.readline().rstrip('\n').split()
            index = int(line) - 24
            for j in range(3):
                first_step[i][j] = float(coords0[j])
# Start the unwrapped coordinates at the origin
        unwrapped = np.copy(first_step)
        prev_step = np.copy(first_step)

# Write out the first step in Cartesian coordinates:
# Multiply the vectors
        xyzcoords = np.multiply(np.dot(first_step,unitcell),scale)
        i = 0
        for el in element_names:
            for j in range(element_dict[el]):
                outputunwrapped.write(el + " " + str(xyzcoords[i][0]) + " " + str(xyzcoords[i][1]) + " " + str(xyzcoords[i][2]) +"\n")
                i += 1

# Read in the next configuration and compute the distance
    line = xdatcar.readline()
    if len(line) == 0:
        break
# Write the number of atoms and the comment line
    outputunwrapped.write(str(Natoms)+"\n"+line)
    for i in range(Natoms):
        coord = xdatcar.readline().rstrip('\n').split() #modified by laa147 to account for the nitrogen atom index in the coordunwrapped.xyz
        coords = xdatcar.readline().rstrip('\n').split()
        index = int(line) - 24
        j = 0
        for j in range(3):
            curr_step[i][j] = float(Ncoord_list[int(index)][j])
            dd = prev_step[i][j] - curr_step[i][j]
            # print(unwrapped)
            # print(dd)
            if np.fabs(dd) > 0.5:
# The distance traveled was larger than 1/2 the box length, so unwrap the coordinate
                unwrapped[i][j] += curr_step[i][j] + np.sign(dd)*1.0 - prev_step[i][j]
            else:
                unwrapped[i][j] += curr_step[i][j] - prev_step[i][j]

# Set the previous step to the current step:
    prev_step = np.copy(curr_step)
# Multiply the vectors
    xyzcoords = np.multiply(np.dot(unwrapped,unitcell),scale)
    wrappedcoords = np.multiply(np.dot(curr_step,unitcell),scale)

# Write the coordinates 
    i = 0
    for el in element_names:
        for j in range(element_dict[el]):
            outputunwrapped.write(el + " " + str(xyzcoords[i][0]) + " " + str(xyzcoords[i][1]) + " " + str(xyzcoords[i][2]) +"\n")
            i += 1
            
xdatcar.close()
outputunwrapped.close()