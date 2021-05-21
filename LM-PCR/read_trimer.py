#!/usr/bin/env python

"""
This program trims (discards) sequences after the END_SEQUENCE (e.g., TAGCA for HIV-1) in a fastq file.

Usage:
python read_trimer.py <fastq_file> <END_SEQUENCE> <output_file>

"""

import sys
argvs = sys.argv

fastq_f = open(argvs[1])
viral_end_seq = argvs[2]
out_f_name =argvs[3]

line_l = []
for line in fastq_f:
  line = line.strip()
  line_l.append(line)

out_f = open(out_f_name,"w")
for i in range(0,len(line_l),4):
  name,seq,info,score = line_l[i:i+4]
  if viral_end_seq in seq:
    viral_end_index = seq.find(viral_end_seq) + len(viral_end_seq)
    seq_genome = seq[viral_end_index:]
    score_genome = score[viral_end_index:]
    if len(seq_genome) > 0:
      res_l = [name,seq_genome,info,score_genome]
      print >> out_f, "\n".join(res_l)

