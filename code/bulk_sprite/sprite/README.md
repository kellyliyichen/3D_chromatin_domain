Calculate intra-TAD ratio from SPRITE data.

1. Convert cluster file to bed file
```
python 1_preprocessing.py $sprite_cluster_file $sprite_bed_file
```

2. Intersect SPRITE data with equal-sized bins and TADs
```
bash 2_intresect_TAD.sh $sprite_bed_file $bins_file $TAD_file
```

3. Calculate intra-TAD ratio
```
python 3_calculate_ratio.py $sprite_tad_file 2 1000
```