#!/usr/bin/env bash

mkdir -p result/enrich_to_random/histone

cat script/IS_sample_list.txt | while read f_name ; do
  #histone
  python script/calcEnrichmentForRandomExpectation.py \
    IS_data/corrected/${f_name}.bed \
    script/cd4t_histone_list.txt \
    chip_data/histone/ \
    > result/enrich_to_random/histone/${f_name}.txt

done
