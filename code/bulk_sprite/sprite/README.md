Calculate intra-TAD ratio from SPRITE data.

1. Convert cluster file to bed file
Input: SPRITE cluster file
Output: SPRITE data in bed format
Usage:
```
python 1_preprocessing.py $sprite_cluster_file $sprite_bed_file
```

2. Intersect SPRITE data with equal-sized bins and TADs
Input: SPRITE data in bed format, genomic bins in bed format, TAD in bed format
Output: SPRITE data with bin and TAD information
Usage:
```
bash 2_intresect_TAD.sh $sprite_bed_file $bins_file $TAD_file
```

3. Calculate intra-TAD ratio
Input: SPRITE data with bin and TAD information, minimum cluster size, maximum cluster size
Output: intra-TAD ratio with average cluster size and number of cluster
Usage:
```
python 3_calculate_ratio.py $sprite_tad_file 2 1000
```