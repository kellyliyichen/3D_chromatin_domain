#!/usr/bin/env python


import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
import sys
import os

def pdist(data, dist_cutoff=500):
    dist = pairwise_distances(data.values, metric='euclidean')
    dist = dist <= dist_cutoff
    return pd.DataFrame(dist)


def pdist_convert(data, slope, intercept, diag, log = True):
    dist = pairwise_distances(data.values, metric='euclidean')
    #n = len(dist)
    #dist[range(n), range(n)] = 1 ## to avoid log0 error
    zeros = np.where(dist==0) ## for some reason, sometimes two different regions have the same coordinates
    dist[zeros] = 1 ## to avoid log0 error
    contact = np.exp(slope * np.log(dist) + intercept)
    #contact[range(n), range(n)] = diag
    contact[zeros] = diag
    if log == True:
        contact = np.log(contact)
    return pd.DataFrame(contact)


def get_tad_bins_intra(df, chrom='21'):
    tad_bin_dict = {}
    all_tad_bins = {}
    for index, row in df.iterrows():
        bin_index = index
        binID, tadID = row['bin_id'], row['tad_id']
        ch = binID.split(':')[0][3:]
        if ch != str(chrom):
            continue

        if bin_index not in all_tad_bins:
            all_tad_bins[bin_index] = (binID, tadID)
        else:
            raise ValueError("Repeated binID index in the TAD file!")

        if tadID not in tad_bin_dict:
            tad_bin_dict[tadID] = [bin_index]
        else:
            tad_bin_dict[tadID].append(bin_index)
    return all_tad_bins, tad_bin_dict



def intraTAD_ratio(all_tad_bins, tad_bin_dict, contact_mat, intraTAD_ratio_file, cell, rep):
    fout = open(intraTAD_ratio_file, 'a')
    for bin_index, ID in all_tad_bins.items():
        binID = ID[0]
        tadID = ID[1]
        all_interactions = contact_mat.iloc[bin_index, :].sum()
        if all_interactions > 0:
            intraTAD_list = tad_bin_dict[tadID]
            intraTAD_interactions = 0
            for i in intraTAD_list:
                intraTAD_interactions += contact_mat.iloc[bin_index, i]
            fout.write(binID + '\t'+ str(intraTAD_interactions/all_interactions) + '\t' + tadID + '\t' + str(cell) + rep + '\n')





if __name__ == '__main__':
    outdir = sys.argv[1]
    rep = sys.argv[2]
    sc_loc_file = sys.argv[3]
    intraTAD_ratio_file = sys.argv[4]

    res = '50000'
    chrom = '21'


    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    data = pd.read_csv(sc_loc_file, sep = '\t', names = ['bin_id', 'bin_tad', 'z', 'x', 'y', 'cell_id'])
    data['tad_id'] = data['bin_tad'].str.split('&', 1).str[1]
    data.dropna(axis=0, inplace=True)

    all_cell = data['cell_id'].unique()
    for cell in all_cell:
        #print(cell)
        sub_data = data[data['cell_id'] == cell]
        sub_data.reset_index(inplace = True)
        n = len(sub_data.index)
        df1 = sub_data.loc[:, ['bin_id', 'tad_id']]
        df2 = sub_data.loc[:, ['z', 'x', 'y']]

        contact_mat = pdist(df2, dist_cutoff=500)

        #print("extract TAD information...")
        all_tad_bins, tad_bin_dict = get_tad_bins_intra(df1, chrom)
        #print("calculate intra-TAD ratio...")
        intraTAD_ratio(all_tad_bins, tad_bin_dict, contact_mat, intraTAD_ratio_file, cell, rep)
        #print("intra-TAD ratio is done!")



