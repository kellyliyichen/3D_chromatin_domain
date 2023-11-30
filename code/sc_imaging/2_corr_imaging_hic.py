#!/usr/bin/env python

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import os

from scipy.stats import pearsonr
from scipy.stats import spearmanr

import pickle
import sys



def average(data, k):
    all_cells = data['cellID'].unique()
    sample_cells = np.random.choice(all_cells, k, replace=False)
    sub_data = data[data['cellID'].isin(sample_cells)]
    res = sub_data.groupby(['binNo']).mean()['ratio'].to_frame()
    res = res.rename(columns={"ratio": str(k)})
    return res


def all_ave(data, k_arr, outfile):
    all_regs = data['binNo'].unique()
    res_df = pd.DataFrame(index=all_regs)
    res_df.index.name = 'binNo'
    for k in k_arr:
        res = average(data, k)
        res_df = pd.merge(res_df, res, left_index=True, right_index=True, how='outer')
    #res_df.to_csv(outfile, sep='\t')
    return res_df


def get_data(hic_ratio_file, image_ratio, k_arr):
    hic_ratio = pd.read_csv(hic_ratio_file, sep='\t',  names=['bin_name', 'hic_ratio'])
    hic_ratio['binNo'] = hic_ratio['bin_name'].str.split('|',expand=True).iloc[:,1]
    join_ratio = pd.merge(image_ratio, hic_ratio, on='binNo', how='outer')
    join_ratio.drop(['bin_name'], axis=1, inplace=True)

    m = len(k_arr)
    
    pcc_arr = np.zeros(m)
    scc_arr = np.zeros(m)
    num_regs = np.zeros(m)
    for i in range(m):
        x = join_ratio.loc[:, str(k_arr[i])]
        y = join_ratio.loc[:, "hic_ratio"]
        nas = np.logical_or(np.isnan(x), np.isnan(y))
        n_reg = len(x[~nas])
        num_regs[i] = n_reg
        if n_reg<3:
            continue
        pcc = pearsonr(x[~nas], y[~nas])[0]
        scc = spearmanr(x[~nas], y[~nas])[0] 
        pcc_arr[i] = pcc
        scc_arr[i] = scc

    #print(num_regs)
    #print(pcc_arr)
    #print(scc_arr)
    
    return pcc_arr, scc_arr, num_regs



if __name__ == "__main__":
    hic_ratio_file = sys.argv[1]
    image_ratio_file = sys.argv[2]
    sample = sys.argv[3]
    outdir = sys.argv[4]

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok = True)

    k_arr = list(range(0,1000,10))[1:]
    #print(len(k_arr))

    print("Read data...")
    data = pd.read_csv(image_ratio_file, sep='\t', names=['binNo', 'ratio', 'tadID', 'cellID'])

    print("Calculate average...")
    res_df = all_ave(data, k_arr, None)

    print("Correlate with Hi-C...")
    pcc_arr, scc_arr, num_regs = get_data(hic_ratio_file, res_df, k_arr)

    pcc_file = outdir + 'pcc_' + sample + '.pkl'
    with open(pcc_file, 'wb') as f:
        pickle.dump(pcc_arr, f) 

    scc_file = outdir + 'scc_' + sample + '.pkl'
    with open(scc_file, 'wb') as f:
        pickle.dump(scc_arr, f) 


