#!/usr/bin/env python

import sys,math
argvs = sys.argv

IS_f_name = argvs[1]
chip_list_f = open(argvs[2])
chip_dir = argvs[3]

######parameters########
chrom_list_f = "script/chrom_list.txt"
distance = 2000
num_of_randomaization = 100
########################


#countOverlapsInObservedData()
def countOverlapsInObservedData(F1_name,F2_name,Chrom_list_f,Distance):
  import subprocess

  Cmd = 'bedtools slop -i %(f2_name)s -g %(chrom_list_f)s -b %(distance)d | \
       bedtools intersect -a "stdin" -b %(f1_name)s | cut -f 1-3 | sort | uniq | wc -l | cut -d " " -f 1' % {'f1_name':F1_name,'f2_name':F2_name,'chrom_list_f':Chrom_list_f,'distance':Distance}
  Count = subprocess.check_output(Cmd, shell=True).strip()
  Count = int(Count)
  return Count

#countOverlapsInRandomizedData()
def countOverlapsInRandomizedData(F1_name,F2_name,Chrom_list_f,Distance,NumRandomization):
  import subprocess
  import numpy as np

  Cmd = 'bedtools shuffle -i %(f2_name)s -g %(chrom_list_f)s -noOverlapping | \
         bedtools slop -i "stdin" -g %(chrom_list_f)s -b %(distance)d | \
         bedtools intersect -a "stdin" -b %(f1_name)s | cut -f 1-3 | sort | uniq | wc -l | cut -d " " -f 1' % {'f1_name':F1_name,'f2_name':F2_name,'chrom_list_f':Chrom_list_f,'distance':Distance}

  Count_l = []
  for i in range(NumRandomization):
    l = []
    Count = subprocess.check_output(Cmd, shell=True).strip()
    Count = int(Count)
    Count_l.append(Count)
  Count_a = np.array(Count_l)
  Mean = np.mean(Count_a)
  Var = np.var(Count_a)
  Std = np.std(Count_a)
  return Mean,Var,Std

#calcPvalUponNormdist()
def calcPvalUponNormdist(Count,Mean,Std):
  from scipy.stats import norm
  if Count >= Mean:
    Pval = norm.sf(x=Count,loc=Mean,scale=Std)*2
  else:
    Pval = norm.cdf(x=Count,loc=Mean,scale=Std)*2
  return Pval

#calcPvalUponPoissondist()
def calcPvalUponPoissondist(Count,Mean):
  from scipy.stats import poisson
  if Count >= Mean:
    Pval = poisson.sf(Count,Mean)*2
  else:
    Pval = poisson.cdf(Count,Mean)*2
  return Pval

chip_f_name_d = {}
for line in chip_list_f:
  chip_f_name = line.strip()
  chip_name = chip_f_name.split('.')[0]
  chip_f_name = chip_dir + chip_f_name
  chip_f_name_d[chip_name] = chip_f_name

num_of_test = len(chip_f_name_d)
print "#Chip_name\tObserved_count\tEnrichment\tPval_norm\tPval_poisson\tFER_norm\tFER_poisson\tMean_in_random\tVar_in_random\tStd_in_random"
for chip_name in chip_f_name_d:
  chip_f_name = chip_f_name_d[chip_name]
  count = countOverlapsInObservedData(IS_f_name,chip_f_name,chrom_list_f,distance)
  mean,var,std = countOverlapsInRandomizedData(IS_f_name,chip_f_name,chrom_list_f,distance,num_of_randomaization)
  if mean > 0:
    enrichment = count/mean
    pval_norm = calcPvalUponNormdist(count,mean,std)
    pval_poisson = calcPvalUponPoissondist(count,mean)
    FER_norm = pval_norm * num_of_test
    FER_poisson = pval_poisson * num_of_test
    l = [chip_name,count,enrichment,pval_norm,pval_poisson,FER_norm,FER_poisson,mean,var,std]
    print "\t".join([str(i) for i in l])
