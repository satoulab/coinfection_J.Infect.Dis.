#!/usr/bin/env bash

mkdir -p chip_data/histone
 
histone mark from Roadmap
cat "script/cd4t_histone_list.txt" | while read data ; do
  wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/${data}.gz -P chip_data/histone
  gzip -d chip_data/histone/${data}.gz
done

#correct bed
mkdir IS_data/corrected
ls IS_data | grep '.bed' | while read f ; do
  sed -E "s/\s+/\t/g" IS_data/${f} | sed -E "s/\t$//g" > IS_data/corrected/${f}
done

ls IS_data/corrected | sed 's/.bed//g' > script/IS_sample_list.txt

#install scipy
pip install scipy --user
