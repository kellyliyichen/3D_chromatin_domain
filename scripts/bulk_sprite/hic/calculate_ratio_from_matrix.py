#!/usr/bin/env python

import sys
import re
import numpy as np
import pandas as pd


def get_chrom_sizes(chrom_size_file, res):
    hg38_sizes = {}
    with open(chrom_size_file, "r") as f:
        for line in f:
            if line.strip() == "":
                continue
            line=line.strip().split("\t")
            hg38_sizes[line[0]] = int(line[1])
    size_ls = [0]
    for chrom in range(1,23):
        size_ls.append(hg38_sizes["chr" + str(chrom)] // int(res) + 1)
    cusize = np.cumsum(size_ls)
    print(cusize)
    return cusize



def get_contact_matrix(res, cusize, contact_pre):
    m = cusize[-1]
    contact_mat = np.zeros((m,m), dtype='float16')
    for i in range(1,23):
        for j in range(i,23):
            print(i, j)
            fileIn = contact_pre + str(i) + "_" + str(j) + ".tab"
            print(fileIn)

            with open(fileIn, 'r') as f:
                for line in f:
                    b1, b2, conf = line.strip().split("\t")
                    bin1 = int(b1) // int(res)
                    bin2 = int(b2) // int(res)
                    bini = bin1 + cusize[i-1]
                    binj = bin2 + cusize[j-1]

                    if conf != "NaN":
                        contact_mat[bini, binj] = float(conf)
                        contact_mat[binj, bini] = float(conf)
    return pd.DataFrame(data=contact_mat)



def get_tad_bins(res, cusize, tad_file):
    tad_bin_dict = {}
    all_tad_bins = {}
    with open(tad_file, 'r') as f:
        for line in f:
            binID, tadID = line.strip().split()
            bin_index = int(re.search('bin(.+?)\|', binID).group(1)) - 1
            ch = tadID.split(':')[0][3:]
            if ch not in list(map(str, range(1,23))):
                continue
            chrom = int(ch)
            bin_index_chrom = bin_index + cusize[chrom-1]
            if bin_index_chrom not in all_tad_bins:
                all_tad_bins[bin_index_chrom] = (binID, tadID)
            else:
                raise ValueError("Repeated binID index in the TAD file!")

            if tadID not in tad_bin_dict:
                tad_bin_dict[tadID] = [bin_index_chrom]
            else:
                tad_bin_dict[tadID].append(bin_index_chrom)
    return all_tad_bins, tad_bin_dict



def intraTAD_ratio(all_tad_bins, tad_bin_dict, contact_mat, intraTAD_ratio_file):
    fout = open(intraTAD_ratio_file, 'w')
    for bin_index, ID in all_tad_bins.items():
        binID = ID[0]
        tadID = ID[1]
        all_interactions = contact_mat.iloc[bin_index, :].sum()
        if all_interactions > 0:
            intraTAD_list = tad_bin_dict[tadID]
            intraTAD_interactions = 0
            for i in intraTAD_list:
                intraTAD_interactions += contact_mat.iloc[bin_index, i]
            fout.write(binID + '\t'+ str(intraTAD_interactions/all_interactions) + '\n')




if __name__ == '__main__':
    chrom_size_file = sys.argv[1]
    contact_pre = sys.argv[2]
    tad_file = sys.argv[3]
    res = sys.argv[4]
    intraTAD_ratio_file = sys.argv[5]

    print("get chromosome size...")
    cusize = get_chrom_sizes(chrom_size_file, res)
    print("extract inter-chromosomal contact matrix...")
    contact_mat = get_contact_matrix(res, cusize, contact_pre)
    print("extract TAD information...")
    all_tad_bins, tad_bin_dict = get_tad_bins(res, cusize, tad_file)
    print("calculate intra-TAD ratio...")
    intraTAD_ratio(all_tad_bins, tad_bin_dict, contact_mat, intraTAD_ratio_file)
    print("intra-TAD ratio is done!")


