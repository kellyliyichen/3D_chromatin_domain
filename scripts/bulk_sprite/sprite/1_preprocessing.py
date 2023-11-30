#!/usr/bin/env python

## Preprocess the sprite data

import sys
import gzip

def cluster_to_bed(fileIn, fileOut):
	"""
	Convert SPRITE cluster file into bed file
	fileIn: input SPRITE cluster file, each row is a cluster containing clusterID, regions in this cluster
	e.g., cluster1	chr12:124313728	chr2:166750477
	fileOut: output bed file, each row is a region with clusterID and clusterSize
	e.g. chr12	124313728	124313729	cluster1	2

	"""

	fout = open(fileOut, 'w')
	with gzip.open(fileIn, 'rt') as f:
	    for line in f:
	        line = line.strip().split()
	        if len(line) > 1001:
	            continue
	        for i in range(1, len(line)):
	            chrom = line[i].split(":")[0]
	            pos = line[i].split(":")[1]
	            size = len(line) - 1
	            fout.write(chrom + "\t" + pos + "\t" + str(int(pos)+1) + "\t" + line[0] + "\t" + str(size) + "\n")
	fout.close()


if __name__ == "__main__":
	fileIn = sys.argv[1]
	fileOut = sys.argv[2]
	cluster_to_bed(fileIn, fileOut)