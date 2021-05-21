#!/bin/bash

####################################################
#                                                  # 
#            LM-PCR Analysis for HTLV-1            #
#                     by Yorifumi                  #  
#                                                  #
####################################################

# ************ Main Script Begins Here ************ #

	echo "Begin analysis"
    echo "Task started at"
    	date
    begin=$(date +%s)

# Make analysis directories
    
    echo "Create analysis directories"
    	mkdir 1_qc
        mkdir 2_map
    echo "Directories 1_qc & 2_map created"
    echo "Begin step 1 in 1_qc directory"
        cd 1_qc

# Remove adaptor sequence

    echo "Remove adaptor from read 1"
    	cutadapt -b GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -O 5 ../*R1.fastq -o R1_step1.fq
    
# Remove HIV-1 LTR sequence
    
    echo "Extract reads containing HTLV-1 LTR sequence, ACACA"
    	read_skipper.pl R1_step1.fq ACACA 
	
	echo "Trim HIV-1 LTR sequence"       
		fastx_trimmer -f 6 -i R1_step1_skip.fq -o R1_noVirus.fq -Q33
	
# Read QC
	
	echo "Reads quality check"
    	qcleaner_renew_v3.1.pl --i1 R1_noVirus.fq --i2 ../*R2.fastq --o1 R1_clean.fastq --o2 R2_clean.fastq --log qclog.txt --trim skip
    
    echo "Step 1 completed"
    echo "Proceed to step 2 in 2_map directory"
    	cd ../2_map

# Mapping to reference genome
    
    echo "Mapping using bwa"
    echo "Mapping against HTLV-1+hg19"
		bwa mem -t 4 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /home/genome/HTLV-1+hg19/genome.fa ../1_qc/R1_clean.fastq ../1_qc/R2_clean.fastq > virus+host.sam
    
    echo "Mapping complete - proceed to extract reads"
    
# Reads extraction and selection
    
    echo "Begin reads extraction and selection"
    	samtools view -Sb virus+host.sam > virus+host.bam
    
    echo "Filter and accept only first mapped reads"
        samtools view -bh -F 256 -o virus+host_uniq.bam virus+host.bam
    
    echo "Sort reads"
        samtools sort virus+host_uniq.bam virus+host_uniq_sort

	echo "Remove unmapped reads"
		samtools view -bh -F 12 -o virus+host_map.bam virus+host_uniq_sort.bam
    
    echo "Sort reads"
    	samtools sort virus+host_map.bam virus+host_map_sort
    
    echo "Convert bam file to sam file"
        samtools view -h virus+host_map_sort.bam > virus+host_map_sort.sam
    
    echo "Remove read-pairs which are more than 1kb apart"
    	awk '((-1001 < $9 && $9 < 1001) || $1~/^@/) {print}' virus+host_map_sort.sam > virus+host_map_sort_ins.sam
    
    echo "Remove read-pairs which are mapped to different chromosomes"
    	awk '($7=="=" || $1~/^@/) {print}' virus+host_map_sort_ins.sam > virus+host_map_sort_ins_chr.sam
    
    echo "Remove noise"
		seq_finder.pl virus+host_map_sort_ins_chr.sam ACACA /home/genome/HTLV-1+hg19/genome.fa > virus+host_map_sort_ins_chr_rm.sam
    
    echo "Extract 3’LTR junction"
		grep -v htlv-1 virus+host_map_sort_ins_chr_rm.sam > host_map_sort_ins_chr_rm.sam
    
    echo "Remove reads with mapping quality 0"
		awk -F"\t" '($5>0 || $1~/^@/) {print}' host_map_sort_ins_chr_rm.sam > host_map_sort_ins_chr_rm_uniq.sam
    
    echo "Generate export file"
		../../../host_export_maker5.pl  host_map_sort_ins_chr_rm_uniq.sam 20 2 > host_export20_2.csv
	
	echo "Calculate copy number"
		../../../copy_analysis_fromExport.pl host_export20_2.csv > copy_number20_2.txt
    
    echo "Data processing complete! Task completed on"
       	date
        end=$(date +%s)
        duration=$(($end-$begin))
        
	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
    
# ************ Main Script Ends Here ************ #
