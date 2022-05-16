import yaml
import numpy as np
import pandas as pd
from skspatial.objects import Line, Plane

def open_line(file_name):
    '''
    Opens one of the text files and parses the line into a skspatial object
    '''
    
    with open(file_name) as f:
        lines = f.readlines()

    point = [float(dim) for dim in lines[0].split()[1:]]
    line = [float(dim) for dim in lines[1].split()[1:]]
    
    return Line(point, line)

def open_carbon_atoms(file_name):    

    with open(file_name) as pdb_structure:
        lines = pdb_structure.readlines()

    first_segment = list()
    second_segment = list()
    third_segment = list()
    fourth_segment = list()

    allCA_7EJ1 = list()

    for line in lines:
        line = line[:22] + ' ' + line[22:] # correct .pdb stupidity        

        if len(line.split()) > 5:
            residue_num = int(line.split()[5])
            chain = line.split()[4]
        else:
            continue

        if 'CA' in line:
            x = float(line.split()[6])
            y = float(line.split()[7])
            z = float(line.split()[8])
            data = [x,y,z]
            allCA_7EJ1.append(data)

        if 455 <= residue_num <= 488 and ('CA' in line):
                
            x = float(line.split()[6])
            y = float(line.split()[7])
            z = float(line.split()[8])
            data = (x,y,z)
            if (chain == 'B'):
                first_segment.append(data)
            elif chain == 'D':
                second_segment.append(data)
            elif chain == 'F':
                third_segment.append(data)
            elif chain == 'H':
                fourth_segment.append(data)    

    return np.array(allCA_7EJ1)

if __name__ == "__main__":
    df = pd.DataFrame(columns=["Structure", "Channel", "Average_radius", "Max_radius", "Min_radius"])
    df_radii = pd.DataFrame(columns=['Structure', 'Channel', 'radii'])


    with open('../CONFIG/config_alignment.yaml', 'r') as f:
        pdb_ids = yaml.load(f, Loader=yaml.FullLoader)

    for channel in pdb_ids.keys():
        for aligned_structure in pdb_ids[channel]:
            for structure in list(aligned_structure.values())[0]:
                print(f"{channel} / {list(aligned_structure.keys())[0]} / {structure}")
                
                # Important hyperparameters
                distance_between_planes = 10
                radius_to_check_for_atoms = 20

                # Open line and carbon atoms for 6V01
                line = open_line(f"line_{list(aligned_structure.keys())[0]}.txt")
                points = open_carbon_atoms(f"../DATA/ALIGNED/{channel}/{structure}.pdb")

                # Generate planes
                origin = [0, 0, 0]
                origin_projection = line.project_point(origin)
                perpendicular_line = Line(origin_projection, origin - origin_projection)
                normal_vector = perpendicular_line.vector.cross(line.vector)
                first_plane = Plane(origin_projection, normal_vector)
                second_plane = Plane(origin_projection, line.vector.cross(normal_vector))                
                plane_points = np.array([line.point + line.vector * i for i in range(-100, 100, distance_between_planes)])
                orthogonal_planes = np.array([Plane(plane_points[i], line.vector) for i in range(plane_points.shape[0])])

                # Separate points into quadrants
                quadrants = np.ndarray((orthogonal_planes.shape[0] - 1, 4), object)
                for i in range(orthogonal_planes.shape[0] - 1):
                    quadrant_1, quadrant_2, quadrant_3, quadrant_4 = [], [], [], []
                    
                    for atom in points:
                        if orthogonal_planes[i].side_point(atom) > 0 and orthogonal_planes[i + 1].side_point(atom) < 0:
                            if line.distance_point(atom) < radius_to_check_for_atoms:
                                if first_plane.side_point(atom) < 0 and second_plane.side_point(atom) > 0:
                                    quadrant_1.append(atom)                    
                                elif first_plane.side_point(atom) > 0 and second_plane.side_point(atom) > 0:
                                    quadrant_2.append(atom)                    
                                elif first_plane.side_point(atom) > 0 and second_plane.side_point(atom) < 0:
                                    quadrant_3.append(atom)                    
                                elif first_plane.side_point(atom) < 0 and second_plane.side_point(atom) < 0:
                                    quadrant_4.append(atom)
                    
                    quadrants[i][0] = quadrant_1
                    quadrants[i][1] = quadrant_2
                    quadrants[i][2] = quadrant_3
                    quadrants[i][3] = quadrant_4
                
                # For each line segment,
                # For each quadrant, find the atom closest to the line. 
                # Using the four atoms identified, calculate the channel's average radius length, maximum radius length, and minimum radius length.
                radii = []
                for segment in quadrants:

                    avg_radius = 0
                    for quadrant in segment:
                        if len(quadrant) < 1:
                            continue

                        smallest_radius = np.min([line.distance_point(atom) for atom in quadrant])
                        avg_radius += smallest_radius

                    avg_radius = avg_radius / 4
                    radii.append(avg_radius)

                # These are features from the channel
                max_radius = np.max([value for value in radii if value > 0])
                min_radius = np.min([value for value in radii if value > 0])
                avg_radius = np.mean([value for value in radii if value > 0])

                # Append to df
                df.loc[len(df.index)] = [structure, channel, avg_radius, max_radius, min_radius] 
                radii = np.array(radii)
                df_radii.loc[len(df_radii)] = [structure, channel, radii]

    
    # Save df as csv
    df.to_csv("./radius.csv", index=False)
    df_radii.to_csv('./radii.csv', index=False)
