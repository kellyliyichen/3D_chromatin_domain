#!/usr/bin/bash

## Intersect SPRITE regions with equal-sized genomic bins and then intersect bins with TAD

sprite_bed_file=$1
bins_file=$2
TAD_file=$3

for res in 5000 10000 25000 50000
do
    sprite_bin_out=${sprite_bed_file}.${res}.bed
    bedtools intersect -a $bins_file -b $sprite_bed_file -sorted -wo | \
    cut -f1-4,8,9 | sort -T ./tmp -u | sort -T ./tmp -k1,1 -k2,2n  > ${sprite_bin_out} &

    bedtools intersect -a $sprite_bin_out -b $TAD_file -f 0.51 -sorted -wo | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}'|\
    #sort -u -T ./tmp > $dir$p.$i.hg38.${res}.bed.domain &    
    sort -u -T ./tmp > ${sprite_bin_out}.domain &
done


