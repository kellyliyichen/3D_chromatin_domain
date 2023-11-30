#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import sys
import os


def get_surface_area(hull):
    # Compute the surface area of the convex hull
    area = 0
    for simplex in hull.simplices:
        vertices = hull.points[simplex]
        area += np.linalg.norm(np.cross(vertices[1]-vertices[0], vertices[2]-vertices[0]))/2
    return area


def get_volume(hull):
    # Compute the volume of the convex hull
    centroid = np.mean(hull.points[hull.vertices], axis=0)
    volume = 0
    for simplex in hull.simplices:
        triangle = hull.points[simplex]
        # Compute the normal vector of the triangle
        normal = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
        # Compute the area of the triangle
        area = np.linalg.norm(normal)/2
        # Compute the distance of the centroid to triangle
        D = -np.dot(normal, triangle[0])
        dist = np.abs(np.dot(normal,centroid)+D)/np.linalg.norm(normal)
        # Update the volume using the distance and area of the triangle
        volume += dist * area/3
    return volume


def cal_convex(points):
    # Compute the convex hull of the points
    hull = ConvexHull(points)
    A = get_surface_area(hull)
    V = get_volume(hull)
    sphericity = np.pi**(1/3) * (6 * V) ** (2/3) / A
    return hull.vertices, sphericity


def read_data(sc_loc_file, rep):
    data = pd.read_csv(sc_loc_file, sep = '\t', names = ['bin_id', 'bin_tad', 'z', 'x', 'y', 'cell_id'])
    data['cell_id'] = data['cell_id'].astype(str) + rep
    data['tad_id'] = data['bin_tad'].str.split('&', 1).str[1]
    data['tad_cell'] = data['tad_id'] + '&' + data['cell_id']
    data['all_id'] = data['bin_tad'] + '&' + data['cell_id']
    data.dropna(axis=0, inplace=True)
    return data




if __name__ == '__main__':
    outdir = './convex_surface/data/'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok = True)

    rep = sys.argv[1]
    sc_loc_file = sys.argv[2]
    file1 = outdir + 'sphericity_' + rep + '.txt'
    file2 = outdir + 'surface_' + rep + '.txt'

    f1 = open(file1, 'w')

    data = read_data(sc_loc_file, rep)
    all_domain = data['tad_cell'].unique()
    res_df = pd.DataFrame(columns=['tad_cell', 'bin_id'])

    for dom in all_domain:
        #print(cell)
        sub_data = data[data['tad_cell'] == dom]
        n = len(sub_data.index)
        if n < 20:
            continue
        df = sub_data.loc[:, ['z', 'x', 'y']]
        df.index = sub_data['all_id'].values
        
        surface_idx, sphericity = cal_convex(df.values)
        f1.write(dom + '\t' + str(sphericity) + '\n') 
        surface_points = sub_data.iloc[surface_idx, :]
        res_df = pd.concat([res_df, surface_points.loc[:, ['tad_cell', 'bin_id']]], ignore_index=False, sort=False) 

    f1.close()
    res_df.to_csv(file2, sep = '\t', index=False, header=True)
