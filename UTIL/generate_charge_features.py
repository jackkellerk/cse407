import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skspatial.objects import Line, Points, Plane

# Important hyperparameters
search_radius = 200
line_potential_resolution = 1

def open_line(file_name):
    '''
    Opens one of the text files and parses the line into a skspatial object
    '''
    
    with open(file_name) as f:
        lines = f.readlines()

    point = [float(dim) for dim in lines[0].split()[1:]]
    line = [float(dim) for dim in lines[1].split()[1:]]
    
    return Line(point, line)

def charged_AA_loc(file_name):    

    with open(file_name) as pdb_structure:
        lines = pdb_structure.readlines()

    axis = open_line('./line_7DTC.txt')

    charged_AA = []
    charges = []

    for line in lines:
        line = line[:22] + ' ' + line[22:] # correct .pdb stupidity        
        line = line.split()

        if len(line) > 5:
            pass
        else:
            continue

        if line[2] == 'CA' and line[3] in list(chargedAAs.keys()):
            x = float(line[6])
            y = float(line[7])
            z = float(line[8])
            point = [x,y,z]
            if axis.distance_point(point) < search_radius:
                charged_AA.append(point)
                charges.append(chargedAAs[line[3]])


       

    return (np.array(charged_AA), charges)


if __name__ == "__main__":
    df = pd.DataFrame(columns=["Structure", "Channel", "Min_potential", "Max_potential", "Max_potential_grad", "Min_potential_grad"])
    potentials = pd.DataFrame(columns=['Structure', 'Channel', 'Line_potential'])

    # create amino acid charge dictionary
    chargedAAs_arr = pd.read_csv('./chargedAAs.csv',header=None).to_numpy()
    chargedAAs = dict()
    for pair in chargedAAs_arr:
        chargedAAs[pair[0]] = pair[1]

    with open('../CONFIG/config_alignment.yaml', 'r') as f:
        pdb_ids = yaml.load(f, Loader=yaml.FullLoader)

    for channel in pdb_ids.keys():
        for aligned_structure in pdb_ids[channel]:
            for structure in list(aligned_structure.values())[0]:
                print(f"{channel} / {list(aligned_structure.keys())[0]} / {structure}")

                # Open line
                line = open_line(f"line_{list(aligned_structure.keys())[0]}.txt")
                
                # get locations and charges of AA within radius
                locs_charges = charged_AA_loc(f"../DATA/ALIGNED/{channel}/{structure}.pdb")


                t = list(range(-100,100,line_potential_resolution))
                plane_points = np.array([line.point + line.vector * i for i in t])
                line_potential = []
                for point in plane_points:
                    potential = 0
                    for i in range(len(locs_charges[1])):
                        dx = locs_charges[0][i][0] - point[0]
                        dy = locs_charges[0][i][1] - point[1]
                        dz = locs_charges[0][i][2] - point[2]
                        r = np.sqrt(dx**2 + dy**2 + dz**2)
                        potential += locs_charges[1][i]/r
                    line_potential.append(potential)

                line_potential = np.array(line_potential)
                line_potential_grad = np.gradient(line_potential)

                min_potential = min(line_potential)
                max_potential = max(line_potential)
                min_potential_grad = min(line_potential_grad)
                max_potential_grad = max(line_potential_grad)

                # Append to df
                df.loc[len(df.index)] = [structure, channel, min_potential, max_potential, min_potential_grad, max_potential_grad] 
                potentials.loc[len(potentials)] = [structure, channel, line_potential]
    
    # Save df as csv
    df.to_csv("./charges.csv", index=False)
    potentials.to_csv('./line_potentials.csv', index=False)
