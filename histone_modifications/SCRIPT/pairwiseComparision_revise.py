#!/usr/bin/env python2.7

import scipy.stats as stats
import sys,math

argvs = sys.argv

IS_A_f_name = argvs[1]
IS_B_f_name = argvs[2]
chip_list_f = open(argvs[3])
chip_dir = argvs[4]

######parameters########
chrom_list_f = "script/chrom_list.txt"
distance = 2000
########################

#countOverlapsInObservedData()
def countOverlapsInObservedData(F1_name,F2_name,Chrom_list_f,Distance):
  import subprocess
  Cmd = 'bedtools slop -i %(f2_name)s -g %(chrom_list_f)s -b %(distance)d | \
       bedtools intersect -a "stdin" -b %(f1_name)s | cut -f 1-3 | sort | uniq | wc -l | cut -d " " -f 1' % {'f1_name':F1_name,'f2_name':F2_name,'chrom_list_f':Chrom_list_f,'distance':Distance}
  Count = subprocess.check_output(Cmd, shell=True).strip()
  Count = int(Count)
  return Count

#countTotalLine()
def countTotalLine(IS_f_name):
  IS_f = open(IS_f_name)
  Count = 0
  for line in IS_f:
    Count += 1
  return Count

chip_f_name_d = {}
for line in chip_list_f:
  chip_f_name = line.strip()
  chip_name = chip_f_name.split('.')[0]
  chip_f_name = chip_dir + chip_f_name
  chip_f_name_d[chip_name] = chip_f_name
num_of_test = len(chip_f_name_d)

print "#Chip_name\tA_hit_count\tB_hit_count\tA_nonhit_count\tB_nonhit_count\tOddsRatio\tPval_log10\tFER_log10"
A_count = countTotalLine(IS_A_f_name)
B_count = countTotalLine(IS_B_f_name)
for chip_name in chip_f_name_d:
  chip_f_name = chip_f_name_d[chip_name]
  A_hit_count = countOverlapsInObservedData(IS_A_f_name,chip_f_name,chrom_list_f,distance)
  B_hit_count = countOverlapsInObservedData(IS_B_f_name,chip_f_name,chrom_list_f,distance)
  A_nonhit_count = A_count - A_hit_count
  B_nonhit_count = B_count - B_hit_count
  oddsRatio, pval = stats.fisher_exact([[A_hit_count, B_hit_count], [A_nonhit_count, B_nonhit_count]])
  FER = pval * num_of_test
  pval_log10 = math.log(pval,10)
  FER_log10 = math.log(FER,10)
  l = [chip_name] + [str(round(i,2)) for i in [A_hit_count,B_hit_count,A_nonhit_count,B_nonhit_count,oddsRatio, pval_log10,FER_log10]]
  print "\t".join(l)
