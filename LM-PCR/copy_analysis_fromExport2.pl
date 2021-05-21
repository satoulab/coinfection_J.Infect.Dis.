#!/usr/bin/perl
##
## copy_analysis_fromExport.pl
## 
## Amelieff Corporation. All Rights Reserved, 2014
##
## exportファイルからcopy_numを算出する
##

use strict;
use warnings;

if(@ARGV != 1){
  print "Usage: perl $0 export-file\n";
  exit 1;
}

my $file = shift @ARGV; # 配列の先頭から要素を取り除きファイル名を設定する

my %key2cnt; 

open my $fh, '<', $file or die;
while(<$fh>){
  chomp;
  next if /^Chr/;
  my @arr = split(/\t/, $_);
  my $key = join("_", $arr[0], $arr[1]);
  $key = join("\t", $key, $arr[3]);
  $key2cnt{$key} ++;
}
close $fh;

print join("\t", "#chr_pos", "Ins_Strand", "copy_num"), "\n";
foreach my $key(sort { $key2cnt{$b} <=>$key2cnt{$a} } keys %key2cnt){
  print join("\t", $key, $key2cnt{$key}), "\n";
}
