#!/usr/bin/perl
##
## read_skipper.pl
## 
## Amelieff Corporation. All Rights Reserved, 2014
## by dedachi kenichi
##
## リードの先頭にウイルスのLTR3'配列5塩基が含まれていないリードを除去する
## USAGE: read_skipper.pl input.fastq TAGCA
##

use strict;
use warnings;

if($ARGV[0] eq ''){
        print "ERROR! Please set input file.\n";
        exit(1);
}

if($ARGV[1] eq ''){
        print "ERROR! Please set virus seq.\n";
        exit(1);
}

open IN, "< $ARGV[0]";
my $input = $ARGV[0];

if($input =~ /\.fastq$/){
        $input =~ s/\.fastq//;
        open OUT, "> $input\_skip.fastq";
}
elsif($input =~ /\.fq$/){
        $input =~ s/\.fq//;
        open OUT, "> $input\_skip.fq";
}
else{
        print " *.fq or *.fastq file is required\n";
        exit(1);
}

my $count=0;
my $skip='yes';
my @arry;

my $allcount = 0;
my $delcount = 0;

while(<IN>){
        s/\r//g;
        chomp;
        $arry[$count]=$_;
        if($count==1){
                if( /^$ARGV[1]/ ){
                        $skip='no';
                }
        }
        $count++;
        if($count==4){
                if($skip eq 'no'){
                        for(my $i=0 ; $i<4 ; $i++){
                                print OUT "$arry[$i]\n";
                        }
                }
                else{
                        $delcount++;
                }
                $allcount++;
                $count = 0;
                $skip='yes';
        }
}

print "Number of all reads = $allcount\n";
print "Number of delete reads = $delcount\n";
