#!/usr/bin/env bash

mkdir -p result/pairwise_comparision/histone

cat script/IS_sample_list.txt | while read f1_name ; do
  cat script/IS_sample_list.txt | while read f2_name ; do
    if [ "${f1_name}" != "${f2_name}" ] ; then
      #histone
      python2.7 script/pairwiseComparision_revise.py \
        IS_data/corrected/${f1_name}.bed \
        IS_data/corrected/${f2_name}.bed \
        script/cd4t_histone_list.txt \
        "chip_data/histone/" \
        > result/pairwise_comparision/histone/${f1_name}_vs_${f2_name}.txt

    fi
  done
done
