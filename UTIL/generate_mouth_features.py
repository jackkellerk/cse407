import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skspatial.objects import Line, Points, Plane
import pandas as pd
import yaml
import urllib.request
import os
import math
import warnings
warnings.filterwarnings('ignore')

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

    allCA = list()

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
            allCA.append(data)   

    return np.array(allCA)

#to get the inner and outer shape of the mouth
def channel_mouth(line, quadrants, first_plane, second_plane):
    segment_atoms= []

 
    for segment in quadrants:
        atoms = []
        for quadrant in segment:
            if len(quadrant) < 1: 
                continue
            #get all the x and y values of each point
            #then get the polar coordinates from each
            for atom in quadrant:
                x = first_plane.distance_point(atom) * first_plane.side_point(atom)
                y = second_plane.distance_point(atom) * second_plane.side_point(atom)
                theta = np.arctan2(x,y)
                r = line.distance_point(atom)
                atoms.append([r,theta])

        #get only the segments that has the points
        if len(atoms) > 1:
            atoms = sorted(atoms, key = lambda x:x[1]) 
            atoms = np.array(atoms)
            segment_atoms.append(atoms)

    segment_atoms = np.array(segment_atoms)
    figure, ax = plt.subplots()
    x_value, y_value = [], []

    #plot the raw shape of the mouth only 
    # for atom in segment_atoms[3]:
    #     x_value.append(atom[0]*math.cos(atom[1]))
    #     y_value.append(atom[0]*math.sin(atom[1]))
    # ax.plot(x_value, y_value, color = "blue")

    #do the theta slicing (we're using pi/6 as delta theta)
    thetas = np.linspace(0, 2*math.pi, 14)
    inner_shape = []
    outer_shape = []

    #get all points with their radius between each theta slices
    for i in range(0, len(thetas)):
        theta_slice = []
        for atom in segment_atoms[3]:
            if thetas[i] < abs(atom[1]) < thetas[i+1]:
                theta_slice.append(atom)
        #to get the smallest and largest radius
        theta_slice = sorted(theta_slice, key=lambda x:x[0])

        if theta_slice:
            inner_shape.append(theta_slice[0]) #the first array has the smallest radius for each slice
            outer_shape.append(theta_slice[-1]) #the last one has the biggest radius for each slice

    inner_shape = np.array(inner_shape)
    avg_smallest_radius = 0
    avg_largest_radius = 0

    #get the average radius for the inner and outer mouth of the channel
    #this where we make the circles for both inner and outer shapes using the average radius
    for x,y in zip(inner_shape, outer_shape):
        avg_smallest_radius += x[0]
        avg_largest_radius += y[0]
    avg_smallest_radius /= len(inner_shape)
    avg_largest_radius /= len(outer_shape)

    #plot the inner mouth
    # x_value, y_value = [], []
    # for i in range(0, len(thetas)):
    #     x_value.append(avg_smallest_radius*math.cos(thetas[i]))
    #     y_value.append(avg_smallest_radius*math.sin(thetas[i]))
    # ax.plot(x_value, y_value, color = "red")

    #plot theouter mouth
    x_value, y_value = [], []
    # for i in range(0, len(thetas)):
    #     x_value.append(avg_largest_radius*math.cos(thetas[i]))
    #     y_value.append(avg_largest_radius*math.sin(thetas[i]))
    # ax.plot(x_value, y_value, color = "red")

    return (avg_smallest_radius, avg_largest_radius)

def create_segments(line, points, first_plane, second_plane):
    # Separate the line/channel into arbitrary size length. Create orthogonal planes to the two perpendicular planes.

    plane_points = np.array([line.point + line.vector * i for i in range(-100, 100, distance_between_planes)])
    orthogonal_planes = np.array([Plane(plane_points[i], line.vector) for i in range(plane_points.shape[0])])

    # Determine what quadrant carbon atoms are in

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

    # print("Number of points per quadrant:")
    # print("N.B. shape = (number of planes x number of quadrants)")
    # for plane in quadrants:
    #     print([len(plane[i]) for i in range(4)])
    return channel_mouth(line, quadrants, first_plane, second_plane)



def do_everything(line, points):
    #this where everything is done, including the plane segments and others
    #create the two perpedniduclular planes
    origin = [0, 0, 0]
    origin_projection = line.project_point(origin)
    perpendicular_line = Line(origin_projection, origin - origin_projection)
    normal_vector = perpendicular_line.vector.cross(line.vector)
    first_plane = Plane(origin_projection, normal_vector)
    second_plane = Plane(origin_projection, line.vector.cross(normal_vector))
 
    return create_segments(line, points, first_plane, second_plane)


if __name__ == "__main__":
    df = pd.DataFrame(columns=["Structure", "Channel", "Inner_mouth_radius", "Outer_mouth_radius"])

    with open('../CONFIG/config_alignment.yaml', 'r') as f:
        pdb_ids = yaml.load(f, Loader=yaml.FullLoader)
        #get all the structures and append to lines and points_arrays
    lines, points_array = [], []
    for channel in pdb_ids.keys():
        for ref_struct in pdb_ids[channel]:
            for struct in list(ref_struct.values())[0]:
                line = open_line(f"line_{list(ref_struct.keys())[0]}.txt")
                points = open_carbon_atoms(f"../DATA/ALIGNED/{channel}/{struct}.pdb")

                # This is the parameter that determines how far apart the orthogonal planes are from each other.
                distance_between_planes = 10

                # This is the parameter that determines the radius from the middle line that is used to check 
                # for atoms to use for calculations for measurement of channel radius.
                radius_to_check_for_atoms = 20

                # # lines = np.array(lines)
                # points_array = np.array(points_array)                
                inner_radius, outer_radius = do_everything(line, points)
                
                #append to df
                df.loc[len(df.index)] = [struct, channel, inner_radius, outer_radius] 

    # Save df as csv
    df.to_csv("./opening_mouth.csv", index=False)
                                    