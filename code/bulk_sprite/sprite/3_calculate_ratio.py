#!/usr/bin/env python

## calculate intra-TAD ratio

import sys


def count_intra_inter_TAD(fileIn):
	"""
	For each region, count the number of intra- and inter-TAD interactions
	fileIn: input file with the following columns
	chrom	start	end	binID	clusterID	clusterSize	TADID
	output file format:
	binID	clusterID	clusterSize	intra inter
	"""
	count_dict = {} #key: cluster, value: [size, num_bins, domain_dict {key: domain, value: [binNo]}]

	with open(fileIn, 'r') as f:
	    for line in f:
	        _, _, _, binNo, cluster, size, domain = line.strip().split()
	        if int(size) > 1000:
	            continue

	        if cluster not in count_dict:
	            count_dict[cluster] = [size, 1, {}]
	        else:
	            count_dict[cluster][1] += 1

	        if domain not in count_dict[cluster][2]:
	            count_dict[cluster][2][domain] = []

	        count_dict[cluster][2][domain].append(binNo)


	fileOut = open(fileIn + ".count", 'w')

	for cluster in count_dict:
	    size = count_dict[cluster][0]
	    num_bins = count_dict[cluster][1]
	    domain_dict = count_dict[cluster][2]
	    for domain, bin_list in domain_dict.items():
	        if domain == ".":
				#intra = 0
				#inter = num_bins - 1
	            continue
	        else:
	            intra = len(bin_list) - 1
	            inter = num_bins - len(bin_list)
	        for binNo in bin_list:
	            #print(i, cluster, size, intra, inter)
	            fileOut.write('\t'.join([binNo, cluster, str(size), str(intra), str(inter)]) + '\n')

	fileOut.close()

	del count_dict



def average_cluster_wsize(fileIn, min_size, max_size):
	"""
	For each bin in each cluster, calculate intra-TAD ratio. 
	Then for each bin, calculate averaged intra-TAD ratio among all clusters involving this bin.
	
	"""
	bin_dict = {}
	with open(fileIn, 'r') as f:
	    for line in f:
	        binNo, _, size, intra, inter = line.strip().split()
	        # only keep cluster size within [min_size, max_size]
	        if int(size) < int(min_size) or int(size) > int(max_size): 
	            continue 
	        norm = int(intra) + int(inter)
	        if norm == 0:
	            continue
	        intra_norm = int(intra) / norm
	        inter_norm = int(inter) / norm
	        if binNo not in bin_dict:
	            bin_dict[binNo] = [1, intra_norm, inter_norm, int(size)]
	        else:
	            bin_dict[binNo][0] += 1
	            bin_dict[binNo][1] += intra_norm
	            bin_dict[binNo][2] += inter_norm
	            bin_dict[binNo][3] += int(size)

	fileOut = open(fileIn + ".average." + min_size + "_" + max_size + "_wsize", 'w')

	for binNo, count in bin_dict.items():
		#if count[0] < 3:
			#continue
	    ave_intra = count[1] / count[0]
	    ave_inter = count[2] / count[0]
	    ave_size = count[3] / count[0]
	    #intra_pro = ave_intra /(ave_intra + ave_inter)
	    #inter_pro = ave_inter /(ave_intra + ave_inter)
	    
	    fileOut.write('\t'.join([binNo, str(ave_intra), str(ave_inter), str(ave_size), str(count[0])]) + '\n')    


	fileOut.close()


if __name__ == '__main__':
	fileIn = sys.argv[1]
	min_size = sys.argv[2]
	max_size = sys.argv[3]
	
	count_intra_inter_TAD(fileIn)
	average_cluster_wsize(fileIn + '.count', min_size, max_size)
