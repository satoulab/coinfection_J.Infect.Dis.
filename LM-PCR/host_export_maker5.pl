#!/usr/bin/perl
##
## host_export_maker.pl
## 
## Amelieff Corporation. All Rights Reserved, 2014
## by dedachi kenichi, Emi Hattori
## update: 2015/01/06
##
## 宿主ゲノムへのマッピング結果からexportファイルを出力する。
##
## USAGE: host_export_maker.pl sam
##
##

use strict;
use warnings;

die "Usage: perl $0 sam score num\n" if @ARGV != 3;

my %pQ = (
'!'=>0, '"'=>1, '#'=>2, '$'=>3, '%'=>4, '&'=>5, "'"=>6, '('=>7, ')'=>8, '*'=>9, '+'=>10,
','=>11, '-'=>12, '.'=>13, '/'=>14, '0'=>15, '1'=>16, '2'=>17, '3'=>18, '4'=>19, '5'=>20,
'6'=>21, '7'=>22, '8'=>23, '9'=>24, ':'=>25, ';'=>26, '<'=>27, '='=>28, '>'=>29, '?'=>30,
'@'=>31, 'A'=>32, 'B'=>33, 'C'=>34, 'D'=>35, 'E'=>36, 'F'=>37, 'G'=>38, 'H'=>39, 'I'=>40, 'J'=>41
);

my %first;
my %other;
my %id2cnt;

my ($sam, $sco, $num) = @ARGV;

open my $SAM, "< $sam" or die;
while(<$SAM>){
  s/\r//g;
  chomp;
  next if /^@/;

  my @splt = split "\t", $_;
  my ($id, $flag, $chr, $pos, $cig, $qua) = ($splt[0], $splt[1], $splt[2], $splt[3], $splt[5], $splt[10]);
$chr =~s/^chr//;
  $id2cnt{$id} ++;

  my $flagBin = sprintf "%b", $flag;

  my $str = $flagBin =~ /1\d{4}$/ ? '-' : '+';
  my $len = &smoke($cig);

  my $lowSco = 0;
  #foreach my $q(split("", substr($qua, 0, 5))){
  foreach my $q(split("", $qua)){
    $lowSco ++ if $pQ{$q} < $sco;
#print "$q\t$pQ{$q}\n";
  }
#print "$qua\n$lowSco\n\n";

  if($flagBin =~ /1\d{6}$/){  # first segment in the template
    $first{$id} = {'chr'=>$chr, 'pos'=>$pos, 'len'=>$len, 'str'=>$str, 'lowSco'=>$lowSco};
  }
  elsif($flagBin =~ /1\d{7}$/){  # last segment in the template
    my $tmp = exists $other{$id} ? $other{$id} : [];
    push @$tmp, {'chr'=>$chr, 'pos'=>$pos, 'len'=>$len, 'str'=>$str, 'lowSco'=>$lowSco};
    $other{$id} = $tmp;
  }
}
close $SAM;

my %pos2cnt;
foreach my $id(keys %first){
  my ($fchr, $fpos, $flen, $fstr, $flow) = ($first{$id}{'chr'}, $first{$id}{'pos'}, $first{$id}{'len'}, $first{$id}{'str'}, $first{$id}{'lowSco'});

  next if $flow >= $num;
  next if $id2cnt{$id} > 2;

  my $frag_len = 0;
  my $frag_sta = $fstr eq '+' ? $fpos : $fpos + $flen;

  next if ! exists $other{$id};

  if(exists $other{$id}){
    my $others = $other{$id};
    foreach my $o(@$others){
      next if $fchr ne $$o{'chr'};
      next if $$o{'lowSco'} >= $num;

      my $tmp = 0;
      if($fstr eq '+' and $$o{'str'} eq '-'){
        $tmp = $$o{'pos'} + $$o{'len'} - $frag_sta;
      }
      elsif($fstr eq '-' and $$o{'str'} eq '+'){
        $tmp = $frag_sta - $$o{'pos'};
      }
      $frag_len = $tmp if $frag_len < $tmp;
    }
  }

  next if $frag_len == 0;

  my $frag_pos = $fchr;
  if($fstr eq '+'){
    $frag_pos .= join("", ":$frag_sta-", $frag_sta + $frag_len - 1);
  }
  else{
    $frag_sta --;
    $frag_pos .= join("", ":", $frag_sta - $frag_len + 1, "-", $frag_sta);
  }
  $pos2cnt{join("\t", $fchr, $frag_sta, $frag_len, $fstr, $frag_pos)} ++;
}

my $TMP = "tmp". int(rand 100);
while(-s $TMP){
  $TMP = "tmp". int(rand 100);
}

open my $fTMP, '>', $TMP or die;
foreach my $pos(keys %pos2cnt){
  print $fTMP join("\t", $pos, $pos2cnt{$pos}), "\n";
}
close $fTMP;

print join("\t", qw/Chr Ins_Start Ins_Len Ins_Strand Position Reads/), "\n";
system "sort -k1,1 -k2,2n -k3,3n $TMP";

unlink $TMP;


sub smoke {
  my $cig = shift @_;
  my $len = 0;
  while($cig =~ /^(\d+)(\D+)(.*)$/){
    $cig = $3;
    next if ($2 ne 'M' and $2 ne 'D');
    $len += $1;
  }
  return $len;
}

