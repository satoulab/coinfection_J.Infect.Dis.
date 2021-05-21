#!/usr/bin/perl
#######################################################################
# filename : qcleaner-renew.pl
# Date     : 2012/07/25
# Updated  : 2013/02/27
# Version  : 3.1
# Brief    : FASTQファイルのクリーニングを行う。
#            【参考】
#
# Author   : Kenichi Dedachi, Amelieff Co. Ltd.
#
# Copyright(C) Amelieff Co. Ltd. 2012

use strict;
use warnings;
use Getopt::Long;
use feature 'state';
use File::Basename;
## ファイルシステムの設定 file::basename用 ##
fileparse_set_fstype('Unix');


#--------Error Message---------#
my $SOFT = "qcleaner_renew.pl";
my $VERSION = "3.1";
my $errorMsg = <<ERROR;
Name: $SOFT
Version: $VERSION
usage for single-end: $SOFT --i <input file> --o <output file> --log <logfile name> --qp <qvalue,percent> -n <number> --trim <qvalue> --length <base num> --fqc_nogroup --skip_fqc --skip_cmp --write
usage for paired-end: $SOFT --i1 <input file 1> --i2 <input file 2> --o1 <output file 1> --o2 <output file 2> --log <logfile name> --qp <qvalue,percent> -n <number> --trim <qvalue> --length <base num>  --fqc_nogroup --skip_fqc --skip_cmp --write

  --i <file>              	Input fastq file (single-end)
  --o <file>              	Output fastq file (single-end)
  --i1 <file>             	Input fastq_1 file (paired-end)
  --i2 <file>             	Input fastq_2 file (paired-end)
  --o1 <file>             	Output fastq_1 file (paired-end)
  --o2 <file>             	Output fastq_2 file (paired-end)
  --log <file>            	Output log file
  --qp <[int,int] or skip>	Check a low quality read, <[quality value, percent]>
  --n  <int or skip>      	Allowable number of N content
  --trim <int or skip>    	Both term bases are trimmed based on this qvalue
  --length <int>          	Check a short read, can not specify 0 or skip
  --fqc_nogroup           	No grouping for FastQC result
  --skip_fqc              	Do not perform FastQC phase
  --skip_cmp              	Do not perform comparing fastq files (paired-end)
  --write                 	Output excluded reads in each phase
  
  Copyright(C) Amelieff Co. Ltd. 2012
ERROR


#---------Define the Global Parameter------#
#--------Checking Argument---------#
## file name ##
my $IN_FILE = '';
my $OUT_FILE = '';
my $IN_1_FILE = '';
my $IN_2_FILE = '';
my $OUT_1_FILE = '';
my $OUT_2_FILE = '';
my $OUT_LOG_FILE = '';
## cleaning option ##
my $qp = '';
my $N_content = '';
my $trim = '';
my $delete_length = '';
## function option ##
my $fqc_nogroup = 0;
my $skip_fqc = 0;
my $skip_cmp = 0;
my $p64 = 0;
my $write = 0;

GetOptions (
	'i=s'=>\$IN_FILE,
	'o=s'=>\$OUT_FILE,
	'i1=s'=>\$IN_1_FILE,
	'i2=s'=>\$IN_2_FILE,
	'o1=s'=>\$OUT_1_FILE,
	'o2=s'=>\$OUT_2_FILE,
	'log=s'=>\$OUT_LOG_FILE,
	'qp=s'=>\$qp,
	'n=s'=>\$N_content,
	'trim=s'=>\$trim,
	'length=s'=>\$delete_length,
	'fqc_nogroup'=>\$fqc_nogroup,
	'skip_fqc'=>\$skip_fqc,
	'skip_cmp'=>\$skip_cmp,
	'write'=>\$write,
) or die $errorMsg;


##sigleかpairか判定##
my ($SINGLE,$PAIRED) = (0, 0);

if ( -f $IN_FILE && ! -f $IN_1_FILE && ! -f $IN_2_FILE  ){
	$SINGLE = 1;
}
elsif ( -f $IN_1_FILE && -f $IN_2_FILE && ! -f $IN_FILE ){
	$PAIRED = 1;
}
else{
	die "helpは -h オプションを指定してください。\n入力ファイル名(シングルorペアエンド)を正しく指定してください\n";
}


#----------Define default thresholds----------#
my $Q_THRE_qp = 20;
my $P_THRE_qp = 80;
my $N_THRE = 6;
my $TRIM_THRE = 20;
my $LEN_THRE = 20;

my $trim_bup;

#----------Option argument treat-----------#
##クリーニング条件オプションに正しい値が与えられたか判定##
die "argument of --qp option is strange\nPlease see USAGE" if ($qp !~ /[0-9]{1,3},[0-9]{1,3}/ && $qp ne 'skip' && $qp ne 'SKIP' && $qp ne '' );
die "argument of --n option is strange\nPlease see USAGE" if ($N_content !~ /[0-9]{1,4}/ && $N_content ne 'skip' && $N_content ne 'SKIP' && $N_content ne '' );
die "argument of --trim option is strange\nPlease see USAGE" if ($trim !~ /[0-9]{1,4}/ && $trim ne 'skip' && $trim ne 'SKIP' && $trim ne '' );
die "Can not skip length option or specify 0\nPlease see USAGE" if ($delete_length eq 0 || $delete_length eq 'skip' || $delete_length eq 'SKIP');
die "argument of --length option is strange\nPlease see USAGE" if ($delete_length !~ /[0-9]{1,4}/ && $delete_length ne '' );


## 出力ファイル名が正常に与えられたか判定##
die "Output file name start with -\nPlease see USAGE" if ( $OUT_FILE =~ /^-/ || $OUT_1_FILE =~ /^-/ || $OUT_2_FILE =~ /^-/ || $OUT_LOG_FILE =~ /^-/ ); #ファイル名冒頭が - だったらエラー
die "Please specify output file name\nPlease see USAGE" if ( $SINGLE && ($OUT_FILE eq '' || $OUT_LOG_FILE eq '') );
die "Please specify output file name\nPlease see USAGE" if ( $PAIRED && ($OUT_1_FILE eq '' || $OUT_2_FILE eq '' || $OUT_LOG_FILE eq '') );



#---------Execute Subroutine & Others------#
&qcleaner;


#------------Define Subroutine-------------#

sub qcleaner{

	## 時間計測のためのtime関数読み込み ##
	my $starttime_recognizable = localtime();
	my $starttime = time();
	my ($in, $in_1, $in_2);		#入力ファイルハンドラ
	my ($input_gz_flag, $input_gz_flag_1, $input_gz_flag_2) = (0, 0, 0);
	
	
	## FastQC結果や中間ファイルなど出力するtempディレクトリを作成、除去リード出力ファイル名取得##
	my $TEMPD_NAME;																	#カレントに作成されるtempディレクトリ名用
	
	my ($fname_save, $fname_save_1, $fname_save_2);
	
	if($SINGLE){
		my $temp = $IN_FILE;
		$TEMPD_NAME = basename($temp);
		$TEMPD_NAME =~ s/\..*/_QCleaner/;

		$temp = $IN_FILE;
		$temp =~ /(.*)\..*/;
		$fname_save = basename($1);
	}
	elsif($PAIRED){
		my $temp = $IN_1_FILE;
		$TEMPD_NAME = basename($temp);
		$TEMPD_NAME =~ s/\..*/_QCleaner/;
		
		$temp = $IN_1_FILE;
		$temp =~ /(.*)\..*/;
		$fname_save_1 = basename($1);
		
		#以下3行は拡張子前ファイル名を取り出すための処理
		$temp = $IN_2_FILE;
		$temp =~ /(.*)\..*/;
		$fname_save_2 = basename($1);
	}	
	
	if( ! -d $TEMPD_NAME){		
		system ( "mkdir $TEMPD_NAME" );
	}


	
	## 入力ファイルopen処理 ##
	if ($SINGLE){
		if($IN_FILE =~ /\.gz$/){		#圧縮ファイルなら展開する
			$input_gz_flag = 1;
#			my $decomp_in_fname = $IN_FILE;
#			$decomp_in_fname =~ s/(.*)\.gz$/$1/;
#			$decomp_in_fname = basename($decomp_in_fname);
#			system( "gzip -c -d $IN_FILE > $TEMPD_NAME/$decomp_in_fname" );		# カレントのtempディレクトリに解凍する
#			$IN_FILE = "$TEMPD_NAME/" . $decomp_in_fname;			#入力ファイル名を解凍後のものに上書き
			open $in, "zcat $IN_FILE |" or die "Cannot open $IN_FILE file\n";
		}
		else{
			open $in, '<', $IN_FILE or die "Cannot open $IN_FILE file\n";
		}
	}
	elsif ($PAIRED){
		if($IN_1_FILE =~ /\.gz$/){		#圧縮ファイルなら展開する
			$input_gz_flag_1 = 1;
#			my $decomp_in_fname = $IN_1_FILE;
#			$decomp_in_fname =~ s/(.*)\.gz$/$1/;
#			$decomp_in_fname = basename($decomp_in_fname);
#			system( "gzip -c -d $IN_1_FILE > $TEMPD_NAME/$decomp_in_fname" );		# カレントのtempディレクトリに解凍する
#			$IN_1_FILE = "$TEMPD_NAME/" . $decomp_in_fname;			#入力ファイル名を解凍後のものに上書き
			open $in_1, "zcat $IN_1_FILE |" or die "1Cannot open $IN_1_FILE file\n";
		}
		else{
			open $in_1, '<', $IN_1_FILE or die "Cannot open $IN_1_FILE file\n";
		}
		
		if($IN_2_FILE =~ /\.gz$/){		#圧縮ファイルなら展開する
			$input_gz_flag_2 = 1;
#			my $decomp_in_fname = $IN_2_FILE;
#			$decomp_in_fname =~ s/(.*)\.gz$/$1/;
#			$decomp_in_fname = basename($decomp_in_fname);
#			system( "gzip -c -d $IN_2_FILE > $TEMPD_NAME/$decomp_in_fname" );		# カレントのtempディレクトリに解凍する
#			$IN_2_FILE = "$TEMPD_NAME/" . $decomp_in_fname;			#入力ファイル名を解凍後のものに上書き
			open $in_2, "zcat $IN_2_FILE |" or die "2Cannot open $IN_2_FILE file\n";
		}
		else{
			open $in_2, '<', $IN_2_FILE or die "Cannot open $IN_2_FILE file\n";
		}
	}
	else{
		die "入力ファイルが正しく開けません\n";
	}
	
	
	## 出力ファイルopen処理 ##
	my ($out, $out_1, $out_2, $out_log, $out_log_valid);			#出力ファイルハンドラ用
	my $OUT_LOG_FILE_VALID;											#validationログ出力用
	my ($dpath_out, $dpath_out_1, $dpath_out_2, $dpath_log);		#ディレクトリまでのパス用
	
	if ($SINGLE){
		$dpath_out = dirname($OUT_FILE);			#dirまでのパス取得
		
		if( ! -d $dpath_out){						#ディレクトリが存在していなければ作成
			system ( "mkdir -p $dpath_out" );
		}
		
		if($OUT_FILE =~ /\.gz$/){					#gz形式だったら
			my $remove_gz_filetail = $OUT_FILE;
			$remove_gz_filetail =~ s/(.*)\.gz$/$1/;
			open $out, '>', "$remove_gz_filetail"; 	#gzを取って、出力ファイルopen
		}
		else{	
			open $out, '>', "$OUT_FILE"; 			#出力ファイルopen
		}	
		
			
	}
	elsif($PAIRED){
		$dpath_out_1 = dirname($OUT_1_FILE);			#dirまでのパス取得
		$dpath_out_2 = dirname($OUT_2_FILE);
	
		if( ! -d $dpath_out_1){						#ディレクトリが存在していなければ作成
			system ( "mkdir -p $dpath_out_1" );
		}
		if( ! -d $dpath_out_2){
			system ( "mkdir -p $dpath_out_2" );
		}
		
		if($OUT_1_FILE =~ /\.gz$/){					#pair_1がgz形式だったら
			my $remove_gz_filetail = $OUT_1_FILE;
			$remove_gz_filetail =~ s/(.*)\.gz$/$1/;
			open $out_1, '>', "$remove_gz_filetail"; 	#gzを取って、出力ファイルopen
		}
		else{	
			open $out_1, '>', "$OUT_1_FILE"; 			#出力ファイルopen
		}	
		
		if($OUT_2_FILE =~ /\.gz$/){					#pair_2がgz形式だったら
			my $remove_gz_filetail = $OUT_2_FILE;
			$remove_gz_filetail =~ s/(.*)\.gz$/$1/;
			open $out_2, '>', "$remove_gz_filetail"; 	#gzを取って、出力ファイルopen
		}
		else{	
			open $out_2, '>', "$OUT_2_FILE"; 			#出力ファイルopen
		}	
		
	}	
	
	## ログファイル出力関係 ##
	$dpath_log = dirname($OUT_LOG_FILE);
	if( ! -d $dpath_log){
		system ( "mkdir -p $dpath_log" );
	}
	open $out_log, '>', "$OUT_LOG_FILE"; 				#ログ出力ファイルopen
	$OUT_LOG_FILE_VALID = $OUT_LOG_FILE . ".valid";		#validation ログ出力ファイル名作成
	open $out_log_valid, '>', "$OUT_LOG_FILE_VALID";	#validlog open
	
	
	

	my ($REMOVED_FNAME, $removed_file);																#除去リード出力ファイル名用
	my ($ID_FNAME, $ID_1_FNAME, $ID_2_FNAME, $TEMPID_1_FNAME, $TEMPID_2_FNAME, $SINGLEDIFF_FNAME, $SIN_REMOVE_FNAME);		#ID出力用ファイル名
	my ($TEMP_FNAME, $TEMP_1_FNAME, $TEMP_2_FNAME);									#TEMP出力用ファイル名
	

	
	
	## 除去リード出力ファイルopen ##
	if( $write == 1){
		if($SINGLE){
			$REMOVED_FNAME = $fname_save . '.removed.fastq';	
			open $removed_file, '>', "$TEMPD_NAME/$REMOVED_FNAME"; 
		}
		elsif($PAIRED){
			$REMOVED_FNAME = $fname_save_1 . '.removed.fastq';	
			open $removed_file, '>', "$TEMPD_NAME/$REMOVED_FNAME"; 
		}
	}
	
	
	## 除去されなかったID出力ファイルopen ##
	my ($id, $id_1, $id_2);
	if($SINGLE){
		$ID_FNAME = $fname_save . '.id';	
		open $id, '>', "$TEMPD_NAME/$ID_FNAME"; 
	}
	elsif($PAIRED){
		$ID_1_FNAME = $fname_save_1 . '.id_1';	
		$ID_2_FNAME = $fname_save_2 . '.id_2';	
		open $id_1, '>', "$TEMPD_NAME/$ID_1_FNAME"; 
		open $id_2, '>', "$TEMPD_NAME/$ID_2_FNAME"; 
	}
	
	
	## 中間処理fastq出力ファイルopen ##
	my ($temp_out, $temp_out_1, $temp_out_2);
		if($SINGLE){
		$TEMP_FNAME = $fname_save . '.temp_fastq';	
		open $temp_out, '>', "$TEMPD_NAME/$TEMP_FNAME"; 
	}
	elsif($PAIRED){
		$TEMP_1_FNAME = $fname_save_1 . '.temp_fastq_1';	
		$TEMP_2_FNAME = $fname_save_2 . '.temp_fastq_2';	
		open $temp_out_1, '+>', "$TEMPD_NAME/$TEMP_1_FNAME"; 
		open $temp_out_2, '+>', "$TEMPD_NAME/$TEMP_2_FNAME"; 
	}
	
	
	## 中間処理ID出力ファイル・シングルfastq出力ファイルopen ##
	my ($tempid_1, $tempid_2, $single_diff, $single_removed);
	if($PAIRED && ! $skip_cmp){
		$TEMPID_1_FNAME = $fname_save_1 . '.temp_id_1';	
	    $TEMPID_2_FNAME = $fname_save_2 . '.temp_id_2';	
	    $SINGLEDIFF_FNAME = $fname_save_1 . '.singlediff';	
	    $SIN_REMOVE_FNAME = $fname_save_1 . '.single_removed';	
	    open $tempid_1, '>', "$TEMPD_NAME/$TEMPID_1_FNAME"; 
	    open $tempid_2, '>', "$TEMPD_NAME/$TEMPID_2_FNAME"; 
		open $single_diff, '+>', "$TEMPD_NAME/$SINGLEDIFF_FNAME"; 
		open $single_removed, '>', "$TEMPD_NAME/$SIN_REMOVE_FNAME"; 
	}
	
	
	
	
	
	### 処理１．ファイルの書式チェック ###
	my @phase_1_dataset = &phase_1_stylecheck($out_log, $in, $in_1, $in_2, $out_log_valid);		#戻り値は配列


	## log 出力 ##
	print $out_log "######################################\n### QCleaner-RENEW.pl  version $VERSION ###\n######################################\n";
	if($SINGLE){
		print $out_log "\nINPUT file name\n[$IN_FILE]\n";
	}
	elsif($PAIRED){
		print $out_log "\nINPUT file name\n[$IN_1_FILE] & [$IN_2_FILE]\n";
	}
		print $out_log "\n処理開始時間：$starttime_recognizable\n\n表記の意味： 処理後／処理前\n";

	
	### 処理２．FastQCの実行 ###
	if($skip_fqc == 0){
		&FastQC($TEMPD_NAME, $IN_FILE, $IN_1_FILE, $IN_2_FILE );		#FastQC実行サブルーチン呼び出し
	}
	
	### 処理３から７を実行するサブルーチン ###
	my ($casava_flag, $del_base_phase6, $del_r);
	my ($casava_flag_1, $del_base_phase6_1, $del_r_1);
	my ($casava_flag_2, $del_base_phase6_2, $del_r_2);
	my ($q, $ratio_p);
	my ($skip_cmp_len, $skip_cmp_len_1, $skip_cmp_len_2);
	
	if($SINGLE){
		($q, $ratio_p, $casava_flag, $del_base_phase6, $del_r, $skip_cmp_len) = &phase_3to7_cleaning($removed_file, $in, $id, $temp_out, $out);		# $del_r は配列のリファレンス
	}
	elsif($PAIRED){
		($q, $ratio_p, $casava_flag_1, $del_base_phase6_1, $del_r_1, $casava_flag_2, $del_base_phase6_2, $del_r_2, $skip_cmp_len_1, $skip_cmp_len_2) = &phase_3to7_cleaning($removed_file, $in, $id, $temp_out, $out, $in_1, $in_2, $id_1, $id_2, $temp_out_1, $temp_out_2, $out_1, $out_2);	# $del_r_1 _2 は配列のリファレンス
	}

	
	### 処理８.ペアを揃える ###
	my (@allbase_p8, @del_r_p8);
	if($PAIRED && ! $skip_cmp){
		($allbase_p8[0], $allbase_p8[1], $del_r_p8[0], $del_r_p8[1]) = &phase_8_single_truncated($TEMPD_NAME, $ID_1_FNAME, $ID_2_FNAME, $TEMPID_1_FNAME, $TEMPID_2_FNAME, $SINGLEDIFF_FNAME, $single_diff, $out_1, $out_2, $temp_out_1, $temp_out_2, $single_removed);
	}
	elsif($skip_cmp != 0 && $PAIRED){
		$allbase_p8[0] = $skip_cmp_len_1;
		$allbase_p8[1] = $skip_cmp_len_2;
	}
	elsif($SINGLE){
		$allbase_p8[0] = $skip_cmp_len;
	}
	
	
	### 出力の処理 ###
	
	# 出力ファイル名末尾が gz であれば、圧縮
	if($SINGLE){
		if($OUT_FILE =~ /\.gz$/){
			my $remove_gz_filetail = $OUT_FILE;
			$remove_gz_filetail =~ s/(.*)\.gz$/$1/;
			system( "gzip -f $remove_gz_filetail" );
		}
	}
	if($PAIRED){
		if($OUT_1_FILE =~ /\.gz$/){
			my $remove_gz_filetail = $OUT_1_FILE;
			$remove_gz_filetail =~ s/(.*)\.gz$/$1/;
			system( "gzip -f $remove_gz_filetail" );
		}
		if($OUT_2_FILE =~ /\.gz$/){
			my $remove_gz_filetail = $OUT_2_FILE;
			$remove_gz_filetail =~ s/(.*)\.gz$/$1/;
			system( "gzip -f $remove_gz_filetail" );
		}
	}
	
	
	##出力ファイルサイズチェック、0ならエラーでdie
	if($SINGLE){
		my $fsize_outfile = -s $OUT_FILE;
		if($fsize_outfile == 0){
			die "出力ファイルサイズが0のためQCleanerを強制終了します。\nPhred score の検出に失敗した可能性があります。\n$OUT_LOG_FILE_VALID をご確認ください。\n";
		}
	}
	if($PAIRED){
		my $fsize_outfile_1 = -s $OUT_1_FILE;
		my $fsize_outfile_2 = -s $OUT_2_FILE;
		if($fsize_outfile_1 == 0 || $fsize_outfile_2 == 0){
			die "出力ファイルサイズが0のためQCleanerを強制終了します。\nPhred score の検出に失敗した可能性があります。\n$OUT_LOG_FILE_VALID をご確認ください。\n";
		}
	}
	
		
	### 処理９．FastQCの実行 ###
		
	if($skip_fqc == 0){
		&FastQC($TEMPD_NAME, $OUT_FILE, $OUT_1_FILE, $OUT_2_FILE);		#FastQC実行サブルーチン呼び出し
	}
	
	
	# 入力ファイルが gz でそれを解凍した場合は、解凍後のファイルを消す
#	if($input_gz_flag == 1){
#		system ( "rm -f $IN_FILE" );
#	}
#	elsif($input_gz_flag_1 == 1){
#		system ( "rm -f $IN_1_FILE" );
#	}
#	elsif($input_gz_flag_2 == 1){
#		system ( "rm -f $IN_2_FILE" );
#	}
	
	# FastQCの結果を取ってくる
	my($fqc_in, $fqc_in_1, $fqc_in_2, $fqc_out, $fqc_out_1, $fqc_out_2);
	
	if(! $skip_fqc){
		if($SINGLE){
			my $SUM_IN = basename($IN_FILE);
			my $SUM_OUT = basename($OUT_FILE);
			
			if($IN_FILE =~ /\.gz$/){
				if($IN_FILE =~ /\.fastq.gz$/){
					$SUM_IN =~ s/(.*)\..*\.gz$/$1_fastqc/;
				}
				else{
					$SUM_IN =~ s/(.*)\.gz$/$1_fastqc/;
				}
			}
			else{
				if($IN_FILE =~ /\.fastq$/){
					$SUM_IN =~ s/(.*)\..*$/$1_fastqc/;
				}
				else{
					$SUM_IN = $SUM_IN . "_fastqc";
				}	
			}
			
			if($OUT_FILE =~ /\.gz$/){
				if($OUT_FILE =~ /\.fastq.gz$/){
					$SUM_OUT =~ s/(.*)\..*\.gz$/$1_fastqc/;
				}
				else{
					$SUM_OUT =~ s/(.*)\.gz$/$1_fastqc/;
				}
			}
			else{
				if($OUT_FILE =~ /\.fastq$/){
					$SUM_OUT =~ s/(.*)\..*$/$1_fastqc/;
				}
				else{
					$SUM_OUT = $SUM_OUT . "_fastqc";
				}	
			}
			
			open $fqc_in, '<', "$TEMPD_NAME/$SUM_IN/summary.txt" or die "Cannot open $TEMPD_NAME/$SUM_IN/summary.txt\n";
			open $fqc_out, '<', "$TEMPD_NAME/$SUM_OUT/summary.txt" or die "Cannot open $TEMPD_NAME/$SUM_OUT/summary.txt\n";
		}
		if($PAIRED){
			my $SUM_1_IN = basename($IN_1_FILE);
			my $SUM_1_OUT = basename($OUT_1_FILE);
			my $SUM_2_IN = basename($IN_2_FILE);
			my $SUM_2_OUT = basename($OUT_2_FILE);
			
			### $IN_1_FILE ###			
			if($IN_1_FILE =~ /\.gz$/){
				if($IN_1_FILE =~ /\.fastq.gz$/){
					$SUM_1_IN =~ s/(.*)\..*\.gz$/$1_fastqc/;
				}
				else{
					$SUM_1_IN =~ s/(.*)\.gz$/$1_fastqc/;
				}
			}
			else{
				if($IN_1_FILE =~ /\.fastq$/){
					$SUM_1_IN =~ s/(.*)\..*$/$1_fastqc/;
				}
				else{
					$SUM_1_IN = $SUM_1_IN . "_fastqc";
				}	
			}
			
			if($OUT_1_FILE =~ /\.gz$/){
				if($OUT_1_FILE =~ /\.fastq.gz$/){
					$SUM_1_OUT =~ s/(.*)\..*\.gz$/$1_fastqc/;
				}
				else{
					$SUM_1_OUT =~ s/(.*)\.gz$/$1_fastqc/;
				}
			}
			else{
				if($OUT_1_FILE =~ /\.fastq$/){
					$SUM_1_OUT =~ s/(.*)\..*$/$1_fastqc/;
				}
				else{
					$SUM_1_OUT = $SUM_1_OUT . "_fastqc";
				}	
			}
			
			
			### IN_2_FILE ###
			if($IN_2_FILE =~ /\.gz$/){
				if($IN_2_FILE =~ /\.fastq.gz$/){
					$SUM_2_IN =~ s/(.*)\..*\.gz$/$1_fastqc/;
				}
				else{
					$SUM_2_IN =~ s/(.*)\.gz$/$1_fastqc/;
				}
			}
			else{
				if($IN_2_FILE =~ /\.fastq$/){
					$SUM_2_IN =~ s/(.*)\..*$/$1_fastqc/;
				}
				else{
					$SUM_2_IN = $SUM_2_IN . "_fastqc";
				}	
			}
			
			if($OUT_2_FILE =~ /\.gz$/){
				if($OUT_2_FILE =~ /\.fastq.gz$/){
					$SUM_2_OUT =~ s/(.*)\..*\.gz$/$1_fastqc/;
				}
				else{
					$SUM_2_OUT =~ s/(.*)\.gz$/$1_fastqc/;
				}
			}
			else{
				if($OUT_2_FILE =~ /\.fastq$/){
					$SUM_2_OUT =~ s/(.*)\..*$/$1_fastqc/;
				}
				else{
					$SUM_2_OUT = $SUM_2_OUT . "_fastqc";
				}	
			}

			open $fqc_in_1, '<', "$TEMPD_NAME/$SUM_1_IN/summary.txt" or die "Cannot open $TEMPD_NAME/$SUM_1_IN/summary.txt\n";
			open $fqc_out_1, '<', "$TEMPD_NAME/$SUM_1_OUT/summary.txt" or die "Cannot open $TEMPD_NAME/$SUM_1_OUT/summary.txt\n";
			open $fqc_in_2, '<', "$TEMPD_NAME/$SUM_2_IN/summary.txt" or die "Cannot open $TEMPD_NAME/$SUM_2_IN/summary.txt\n";
			open $fqc_out_2, '<', "$TEMPD_NAME/$SUM_2_OUT/summary.txt" or die "Cannot open $TEMPD_NAME/$SUM_2_OUT/summary.txt\n";
		}
	}
	
	
	
	
	## ログへの出力
	#phase 2
	
	if(! $skip_fqc){
		if($SINGLE){
			print $out_log "\n----------\nPhase_2（処理前のFastQC実行結果）\n$IN_FILE\n";
			while(<$fqc_in>){
				print $out_log "$_";
			}
		}
		elsif($PAIRED){
			print $out_log "\n----------\nPhase_2（処理前のFastQC実行結果）\n$IN_1_FILE\n";
			while(<$fqc_in_1>){
				print $out_log "$_";
			}
			print $out_log "\n$IN_2_FILE\n";
			while(<$fqc_in_2>){
				print $out_log "$_";
			}
		}
	}	
	elsif($skip_fqc){
		print $out_log "\n----------\nPhase_2（処理前のFastQC実行結果）\n phase_2 はスキップされました\n";
	}
	
	
	my ($rem_r, $rem_r_1, $rem_r_2);
	my ($rem_b, $rem_b_1, $rem_b_2);
	
	#phase 3
	print $out_log "\n----------\nPhase_3（CASAVA filter [:Y:] のリード除去）\n";
	
	if($SINGLE){
		if($casava_flag == 1){
			my $temp_r = $phase_1_dataset[1] - @$del_r[3];
			print $out_log "$IN_FILE\n[リード数変化] $temp_r／$phase_1_dataset[1], @$del_r[3] リードを除去\n";
			$rem_r = $temp_r;
		}
		elsif($casava_flag == 0){
			$rem_r = $phase_1_dataset[1];
			print $out_log "$IN_FILE\n CASAVAフィルタは適応されていません\n";
		}
		else{
			print "CASAVA関連処理に異常があります\n";
			die;
		}
	}
	elsif($PAIRED){
		if($casava_flag_1 == 1){
			my $temp_r = $phase_1_dataset[1] - @$del_r_1[3];
			print $out_log "$IN_1_FILE\n[リード数変化] $temp_r／$phase_1_dataset[1], @$del_r_1[3] リードを除去\n";
			$rem_r_1 = $temp_r;
		}
		elsif($casava_flag_1 == 0){
			$rem_r_1 = $phase_1_dataset[1];
			print $out_log "$IN_1_FILE\n CASAVAフィルタは適応されていません\n";
		}
		else{
			print "CASAVA関連処理に異常があります\n";
			die;
		}
		
		if($casava_flag_1 == 1){
			my $temp_r = $phase_1_dataset[3] - @$del_r_2[3];
			print $out_log "\n$IN_2_FILE\n[リード数変化] $temp_r／$phase_1_dataset[1], @$del_r_2[3] リードを除去\n";
			$rem_r_2 = $temp_r;
		}
		elsif($casava_flag_1 == 0){
			$rem_r_2 = $phase_1_dataset[3];
			print $out_log "$IN_2_FILE\n CASAVAフィルタは適応されていません\n";
		}
	}
	
	
	#phase 4
	my $hyaku_p = $ratio_p * 100;
	
	print $out_log "\n----------\nPhase_4（クオリティが$q未満の塩基が$hyaku_p%以上存在するリードを除去)\n";
	if($qp eq 'skip' || $qp eq 'SKIP'){
		print $out_log "phase_4 はスキップされました。\n"
	}
	else{
		if($SINGLE){
			my $temp_r = $rem_r - @$del_r[4];
			print $out_log "$IN_FILE\n[リード数変化] $temp_r／$rem_r, @$del_r[4] リードを除去\n\n";
			$rem_r = $temp_r;
		}
		elsif($PAIRED){
			my $temp_r = $rem_r_1 - @$del_r_1[4];
			print $out_log "$IN_1_FILE\n[リード数変化] $temp_r／$rem_r_1, @$del_r_1[4] リードを除去\n\n";
			$rem_r_1 = $temp_r;
			$temp_r = $rem_r_2 - @$del_r_2[4];
			print $out_log "$IN_2_FILE\n[リード数変化] $temp_r／$rem_r_2, @$del_r_2[4] リードを除去\n\n";
			$rem_r_2 = $temp_r;
		}
	}
	
	
	#phase 5
	print $out_log "\n----------\nPhase_5（未知の塩基（N）を一定量含むリードを除去）\n";
	if($N_content eq 'skip' || $N_content eq 'SKIP'){
		print $out_log "phase_5 はスキップされました。\n"
	}
	else{
		if($SINGLE){
			my $temp_r = $rem_r - @$del_r[5];
			print $out_log "$IN_FILE\n[リード数変化] $temp_r／$rem_r, @$del_r[5] リードを除去\n\n";
			$rem_r = $temp_r;
		}
		elsif($PAIRED){
			my $temp_r = $rem_r_1 - @$del_r_1[5];
			print $out_log "$IN_1_FILE\n[リード数変化] $temp_r／$rem_r_1, @$del_r_1[5] リードを除去\n\n";
			$rem_r_1 = $temp_r;
			$temp_r = $rem_r_2 - @$del_r_2[5];
			print $out_log "$IN_2_FILE\n[リード数変化] $temp_r／$rem_r_2, @$del_r_2[5] リードを除去\n\n";
			$rem_r_2 = $temp_r;
		}
	}
	
	
	#phase 6
	print $out_log "\n----------\nPhase_6（リード両末端のクオリティが$trim_bup未満の塩基をトリミング）\n";
		if($trim eq 'skip' || $trim eq 'SKIP'){
		print $out_log "phase_6 はスキップされました。\n"
	}
	else{
		if($SINGLE){
			my $temp_b = $phase_1_dataset[0] - $del_base_phase6;
			print $out_log "$IN_FILE\n[塩基数変化] $temp_b／$phase_1_dataset[0], $del_base_phase6 塩基を除去\n\n";
			$rem_b = $temp_b;
		}
		elsif($PAIRED){
			my $temp_b = $phase_1_dataset[0] - $del_base_phase6_1;
			print $out_log "$IN_1_FILE\n[塩基数変化] $temp_b／$phase_1_dataset[0], $del_base_phase6_1 塩基を除去\n\n";
			$rem_b_1 = $temp_b;
			$temp_b = $phase_1_dataset[2] - $del_base_phase6_2;
			print $out_log "$IN_2_FILE\n[塩基数変化] $temp_b／$phase_1_dataset[2], $del_base_phase6_2 塩基を除去\n\n";
			$rem_b_2 = $temp_b;
		}
	}
	
	
	#phase 7
	print $out_log "\n----------\nPhase_7（リード長が$delete_length未満のリードを除去）\n";
	if($SINGLE){
		my $temp_r = $rem_r - @$del_r[7];
		print $out_log "$IN_FILE\n[リード数変化] $temp_r／$rem_r, @$del_r[7] リードを除去\n\n";
		$rem_r = $temp_r;
	}
	elsif($PAIRED){
		my $temp_r = $rem_r_1 - @$del_r_1[7];
		print $out_log "$IN_1_FILE\n[リード数変化] $temp_r／$rem_r_1, @$del_r_1[7] リードを除去\n\n";
		$rem_r_1 = $temp_r;
		$temp_r = $rem_r_2 - @$del_r_2[7];
		print $out_log "$IN_2_FILE\n[リード数変化] $temp_r／$rem_r_2, @$del_r_2[7] リードを除去\n\n";
		$rem_r_2 = $temp_r;
	}
	
	
	#phase 8
	if($PAIRED && ! $skip_cmp){
		print $out_log "\n----------\nPhase_8（ペアエンドを揃える）\n";
		my $temp_r = $rem_r_1 - $del_r_p8[0];
		print $out_log "$IN_1_FILE\n[リード数変化] $temp_r／$rem_r_1, $del_r_p8[0] リードを除去\n\n";
		$rem_r_1 = $temp_r;
		$temp_r = $rem_r_2 - $del_r_p8[1];
		print $out_log "$IN_2_FILE\n[リード数変化] $temp_r／$rem_r_2, $del_r_p8[1] リードを除去\n\n";
		$rem_r_2 = $temp_r;
	}
	elsif($skip_cmp && $PAIRED){
		print $out_log "\n----------\nPhase_8（ペアエンドを揃える）\n入力データはペアエンドリードですが、phase8はスキップされました\n\n";
	}
	elsif($SINGLE){
		print $out_log "\n----------\nPhase_8（ペアエンドを揃える）\n入力データがシングルリードのため、phase8はスキップされました\n\n";
	}
	
	
	#phase 9
	
	if(! $skip_fqc){
		if($SINGLE){
			print $out_log "\n----------\nPhase_9（処理後のFastQC実行結果）\n$IN_FILE\n";
			while(<$fqc_out>){
				print $out_log "$_";
			}
		}
		elsif($PAIRED){
			print $out_log "\n----------\nPhase_9（処理後のFastQC実行結果）\n$IN_1_FILE\n";
			while(<$fqc_out_1>){
				print $out_log "$_";
			}
			print $out_log "\n$IN_2_FILE\n";
			while(<$fqc_out_2>){
				print $out_log "$_";
			}
		}
	}	
	elsif($skip_fqc){
		print $out_log "\n----------\nPhase_9（処理後のFastQC実行結果）\n phase_9 はスキップされました\n";
	}
	
	
	my $endtime_recognizable = localtime();
	my $endtime = time();
	
	print $out_log "\n----------\n処理終了時間\n $endtime_recognizable\n\n";
	
	if($SINGLE){
		my $rem_r_ratio = ($rem_r / $phase_1_dataset[1]) * 100;
		my $rem_b_ratio = ($allbase_p8[0] / $phase_1_dataset[0]) * 100;
		
		print $out_log "--------------------\n最終レポート\n（処理後/処理前リード数, リード残存割合, 処理後/処理前塩基数, 塩基残存割合）\n";
		print $out_log "$IN_FILE   $rem_r/$phase_1_dataset[1] ";
		printf $out_log ("%5.2f%%   ", $rem_r_ratio);
		print $out_log "$allbase_p8[0]/$phase_1_dataset[0] ";
		printf $out_log ("%5.2f%%\n\n", $rem_b_ratio);
	}
	elsif($PAIRED){
		my $rem_r_ratio_1 = ($rem_r_1 / $phase_1_dataset[1]) * 100;
		my $rem_r_ratio_2 = ($rem_r_2 / $phase_1_dataset[3]) * 100;
		my $rem_b_ratio_1 = ($allbase_p8[0] / $phase_1_dataset[0]) * 100;
		my $rem_b_ratio_2 = ($allbase_p8[1] / $phase_1_dataset[2]) * 100;
		
		print $out_log "--------------------\n最終レポート\n（処理後/処理前リード数, リード残存割合, 処理後/処理前塩基数, 塩基残存割合）\n";
		print $out_log "$IN_1_FILE   $rem_r_1/$phase_1_dataset[1] ";
		printf $out_log ("%5.2f%%   ", $rem_r_ratio_1);
		print $out_log "$allbase_p8[0]/$phase_1_dataset[0] ";
		printf $out_log ("%5.2f%%\n\n", $rem_b_ratio_1);
		
		print $out_log "$IN_2_FILE   $rem_r_2/$phase_1_dataset[3] ";
		printf $out_log ("%5.2f%%   ", $rem_r_ratio_2);
		print $out_log "$allbase_p8[1]/$phase_1_dataset[2] ";
		printf $out_log ("%5.2f%%\n\n", $rem_b_ratio_2);
	}
	
	my $perform_time = $endtime - $starttime;
	print $out_log "（処理に要した時間）\n $perform_time sec\n";
	
	close $out_log;
	
	
	
	
	## 除去されなかったID出力ファイル削除 ##
	if($SINGLE){
		close $id;
		system ( "rm -f $TEMPD_NAME/$ID_FNAME" ); 
	}
	elsif($PAIRED){
		close $id_1;
		close $id_2;
		system ( "rm -f $TEMPD_NAME/$ID_1_FNAME"); 
		system ( "rm -f $TEMPD_NAME/$ID_2_FNAME"); 
	}
	
	
	## 中間処理fastq出力ファイルclose ##
	if($SINGLE){
		close $temp_out;
		system( "rm -f $TEMPD_NAME/$TEMP_FNAME");
	}
	elsif($PAIRED){
		close $temp_out_1;
		close $temp_out_2;
		system( "rm -f $TEMPD_NAME/$TEMP_1_FNAME");
		system( "rm -f $TEMPD_NAME/$TEMP_2_FNAME");
	}
	
	
	## 中間処理ID出力ファイル・シングルfastq出力ファイルclose ##

	if($PAIRED && ! $skip_cmp){
	    close $tempid_1;
	    close $tempid_2;
		close $single_diff;
		close $single_removed;
		system ( "rm -f $TEMPD_NAME/$TEMPID_1_FNAME"); 
		system ( "rm -f $TEMPD_NAME/$TEMPID_2_FNAME"); 
		system ( "rm -f $TEMPD_NAME/$SINGLEDIFF_FNAME"); 
	#	system ( "rm -f $TEMPD_NAME/$SIN_REMOVE_FNAME"); 
	}
	
}

	
	
sub phase_8_single_truncated {
	
	my ( $TEMPD_NAME, $ID_1_FNAME, $ID_2_FNAME, $TEMPID_1_FNAME, $TEMPID_2_FNAME, $SINGLEDIFF_FNAME, $single_diff, $out_1, $out_2, $temp_out_1, $temp_out_2, $single_removed ) = @_;
	
	system( "env TMPDIR=$TEMPD_NAME/ sort $TEMPD_NAME/$ID_1_FNAME > $TEMPD_NAME/$TEMPID_1_FNAME" );
	system( "env TMPDIR=$TEMPD_NAME/ sort $TEMPD_NAME/$ID_2_FNAME > $TEMPD_NAME/$TEMPID_2_FNAME" );
	system( "env TMPDIR=$TEMPD_NAME/ diff -u $TEMPD_NAME/$TEMPID_1_FNAME $TEMPD_NAME/$TEMPID_2_FNAME > $TEMPD_NAME/$SINGLEDIFF_FNAME" );
	
	
	<$single_diff>; # diff結果の最初の三行はヘッダなので飛ばす
	<$single_diff>;
	<$single_diff>;
	
	my %array_single_1;
	my %array_single_2;
	
	while (<$single_diff>){   ##シングルのIDのみ配列に書き出す（メモリが心配）
		
		s/\r//g;
		chomp;
		next if ( /^$/ );
		my $tmp = $_;
		
		if ($tmp =~ /^-/){
			substr($tmp,0,1) = '';
			$array_single_1{$tmp} = '1';
		}
		elsif ($tmp =~ /^\+/){
			substr($tmp,0,1) = '';
			$array_single_2{$tmp} = '1';
		}
		else{
		}
	}		
	
	my ($finalallbase_1, $del_1_p8) = &pair_check($temp_out_1, $out_1, \%array_single_1, $single_removed);
	my ($finalallbase_2, $del_2_p8) = &pair_check($temp_out_2, $out_2, \%array_single_2, $single_removed);
	
	return ($finalallbase_1, $finalallbase_2, $del_1_p8, $del_2_p8);
}



sub pair_check {
	
	my ($input, $output, $array_single, $single_removed) = @_; 
	my @array;
	my ($cnt, $del_r, $cnt_allbase) = (0, 0, 0);
	
	seek($input, 0, 0);
	
	while(<$input>){
	
		s/\r//g;
		chomp;
		next if ( /^$/ );
		$array[$cnt] = $_;
		$cnt++;
		
		next if $cnt < 4;		#ここまでで4行だけ読み込む
		
		$cnt = 0;
		
		my $q0 = $array[0];
		
		if ($q0 =~ /\s/){
			$q0 =~ /(.*\s)/;
			$q0 = $1;
		}
		elsif ($q0 =~ /\//){
			$q0 =~ /(.*\/)/;
			$q0 = $1;
		}
		else{
			if($PAIRED){
				die "ID行にスペースまたは/が含まれていません。IDの抽出ができないため終了します\n";
			}
		}
		
		if (defined($$array_single{$q0})){
			delete $$array_single{$q0};
			$del_r++;
			print $single_removed "$array[0]\n";
		    print $single_removed "$array[1]\n";
		   	print $single_removed "$array[2]\n";
			print $single_removed "$array[3]\n";
		}
		else{
			my $len = length($array[1]);
			$cnt_allbase += $len;
			print $output "$array[0]\n";
			print $output "$array[1]\n";
			print $output "$array[2]\n";
			print $output "$array[3]\n";
		}
	}
	return ($cnt_allbase, $del_r);
}



sub phase_3to7_cleaning {
	
	my ($removed_file, $in, $id, $temp_out, $out, $in_1, $in_2, $id_1, $id_2, $temp_out_1, $temp_out_2, $out_1, $out_2) = @_;
	
	my $cnt = 0;
	my ($q, $p, $ratio_p);			# --qp オプション用の変数
	
	if($qp ne 'skip' && $qp ne 'SKIP' && $qp ne '' ){		# --qpオプションをparseする
		my @temp = split /,/, $qp;
		$q = $temp[0];
		$p = $temp[1];
		$ratio_p = (100 - $p) * 0.01;
	}
	elsif($qp eq '' || $qp eq 'skip' || $qp eq 'SKIP'){				# --qp オプション指定なしでデフォルト値
		$q = $Q_THRE_qp;
	    $p = $P_THRE_qp;
	    $ratio_p = (100 - $p) * 0.01;
	}
	
	if($N_content eq ''){			# --n オプション指定なしでデフォルト値
		$N_content = $N_THRE;
	}
	
	if($trim eq ''){				# --trim オプション指定なしでデフォルト値
		$trim = $TRIM_THRE;
	}
	
	$trim_bup = $trim;
	
	if($delete_length eq ''){		# --length オプション指定なしでデフォルト値
		$delete_length = $LEN_THRE;
	}
	
	
	# $phred64のinputの場合の処理
	my $q_thre;
	
	if($p64 == 1){			
		$q_thre = $q + 64;
		if( $trim ne 'skip' && $trim ne 'SKIP' ){
			$trim = $trim + 64;
		}
	}
	else{
		$q_thre = $q + 33;
		if( $trim ne 'skip' && $trim ne 'SKIP' ){
			$trim = $trim + 33;
		}
	}	
	
	my ($casava_flag, $del_base_phase6, @del_r);
	my ($casava_flag_1, $del_base_phase6_1, @del_r_1);
	my ($casava_flag_2, $del_base_phase6_2, @del_r_2);
	
	my ($skip_cmp_len, $skip_cmp_len_1, $skip_cmp_len_2);
	### クリーニング開始 ###
	if($SINGLE){
		if($IN_FILE =~ /\.gz$/){
			close $in;
			open $in, "zcat $IN_FILE |" or die "Cannot open $IN_FILE file\n";
		}
		else{
			seek($in, 0, 0);		#ファイルポインタを先頭に
		}
		my @state_init = (1,1,1,1,1);
		($casava_flag, $del_base_phase6, $skip_cmp_len, @del_r) = &cleaning($in, $out, $temp_out, $id, $removed_file, $q_thre, $trim, $ratio_p, @state_init);
		
		return($q, $ratio_p, $casava_flag, $del_base_phase6, \@del_r, $skip_cmp_len);
	}
	elsif($PAIRED){
		if($IN_1_FILE =~ /\.gz$/){
			close $in_1;
			open $in_1, "zcat $IN_1_FILE |" or die "Cannot open $IN_1_FILE file\n";
		}
		else{
			seek($in_1, 0, 0);		#ファイルポインタを先頭に
		}
		my @state_init = (1,1,1,1,1);
		($casava_flag_1, $del_base_phase6_1, $skip_cmp_len_1, @del_r_1) = &cleaning($in_1, $out_1, $temp_out_1, $id_1, $removed_file, $q_thre, $trim, $ratio_p, @state_init);

		if($IN_2_FILE =~ /\.gz$/){
			close $in_2;
			open $in_2, "zcat $IN_2_FILE |" or die "Cannot open $IN_2_FILE file\n";
		}
		else{
			seek($in_2, 0, 0);		#ファイルポインタを先頭に
		}
		@state_init = (1,1,1,1,1);
		($casava_flag_2, $del_base_phase6_2, $skip_cmp_len_2, @del_r_2) = &cleaning($in_2, $out_2, $temp_out_2, $id_2, $removed_file, $q_thre, $trim, $ratio_p, @state_init);
		return($q, $ratio_p, $casava_flag_1, $del_base_phase6_1, \@del_r_1, $casava_flag_2, $del_base_phase6_2, \@del_r_2, $skip_cmp_len_1, $skip_cmp_len_2);
	}
}



sub cleaning {
	
	my ($input, $out_nocmp, $output, $id_output, $removed_file, $q_thre, $trim, $ratio_p, @state_init) = @_;
	my @array;
	
	my ($cnt, $del_base) = (0, 0);
	my $casava_flag = 0;
	my @del_r;
	
	my $skip_cmp_len = 0;		# skip_cmpの時だけ全塩基数計測する
	
	while(<$input>){
	
		s/\r//g;
		chomp;
		next if ( /^$/ );
		$array[$cnt] = $_;
		$cnt++;
		next if $cnt < 4;
		
		my $base_len = length ($array[1]);			#リードの長さ計測
		#削除されたらフラグ立てる。配列は各処理での削除されたリード数を保存。わかりやすいように処理番号と配列indexを揃える
		my $del_flag = 0;
					
		### 処理３．CASAVA：Yの除去 ###
		( $casava_flag, $del_flag, $del_r[3] ) = &phase_3_casava_check($state_init[0], $removed_file, \@array);
		$state_init[0] = 0;
		
		### 処理４．--qpの処理（低クオリティリードの除外）###
		if ( $del_flag == 0 && $qp ne 'skip' && $qp ne 'SKIP' ){
			( $del_flag, $del_r[4] ) = &phase_4_qp_check($state_init[1], $base_len, $ratio_p, $q_thre, $removed_file, \@array);
			$state_init[1] = 0;
		}
		
		### 処理５．--nの処理(Nが多いリードの除外) ###
		if ( $del_flag == 0 && $N_content ne 'skip' && $N_content ne 'SKIP' ){
			( $del_flag, $ del_r[5] ) = &phase_5_N_check($state_init[2], $removed_file, \@array);
			$state_init[2] = 0;
		}
		
		### 処理６．--trim の処理(トリミング) ###
		if ( $del_flag == 0 && $trim ne 'skip' && $trim ne 'SKIP' ){
			( $del_base ) = &phase_6_trim_check($state_init[3], $trim, \@array);
			$state_init[3] = 0;
		}
		
		### 処理７．--lengthの処理（短いリード除外）###
		if ( $del_flag == 0 ){
			( $del_flag, $del_r[7] ) = &phase_7_length_check($state_init[4], $removed_file, \@array);
			$state_init[4] = 0;
		}
		
		if($del_flag == 0){
		
			my $id_line = $array[0];

			if ($id_line =~ /\s/){
				$id_line =~ /(.*\s)/;
				print $id_output "$1\n";
			}
			elsif ($id_line =~ /\//){
				$id_line =~ /(.*\/)/;
				print $id_output "$1\n";
			}
			else{
				if($PAIRED){
				die "ID行にスペースまたは/が含まれていません。IDの抽出ができないため終了します\n";
				}
			}
			
			if($skip_cmp != 0 || $SINGLE){
			
				$skip_cmp_len += length($array[1]);
				
				for( my $i = 0 ; $i < 4 ; $i++ ){
					print $out_nocmp "$array[$i]\n";
				}
			}
			else{
				for( my $i = 0 ; $i < 4 ; $i++ ){
					print $output "$array[$i]\n";
				}
			}
		}
		$cnt = 0;
	}

	return ($casava_flag, $del_base, $skip_cmp_len, @del_r);	
}
			

			
sub phase_7_length_check{
	 
	my ($state_init, $removed_file, $array) = @_;
	
	state $del_r_num = 0;
	my $del_flag = 0;

    ## state initialize ##
	if ( $state_init != 0 ){
		$del_r_num = 0;
	}
	
	my $cnt_base = length($$array[1]);				# 塩基数カウント
	
	if ( $cnt_base < $delete_length ){
			$del_r_num++;
			$del_flag = 1;

		if ($write == 1){
			print $removed_file "$$array[0]\n";	#excludedに出力
			print $removed_file "$$array[1]\n";
			print $removed_file "$$array[2]\n";
	    	print $removed_file "$$array[3]\n";
		}
	}
	
	return ($del_flag, $del_r_num);
}
	
	
			
sub phase_6_trim_check{
	
	my ($state_init, $trim, $array) = @_;
	
	state $delete_base_num = 0;
	
	if($state_init != 0){
		$delete_base_num = 0;
	}
	
	my @base_items = split //, $$array[1];			#塩基を一文字ずつ配列に格納
	my @quality_items = split //, $$array[3];		#クオリティ値を一文字ずつ配列に格納

	my @quality_items_bup = @quality_items;
	
	
	foreach my $tmp (@quality_items_bup){
		my $ascii = ord $tmp;						#クオリティ値をASCIIコードに変換
		last if ($ascii >= $trim);				#クオリティの閾値の判定
		shift (@quality_items);
		shift (@base_items);
		$delete_base_num++;
					
	}
	
	@base_items = reverse(@base_items);
	@quality_items = reverse(@quality_items);
	
	@quality_items_bup = @quality_items;
	
	foreach my $tmp (@quality_items_bup){
		my $ascii = ord $tmp;						#クオリティ値をASCIIコードに変換
		last if ($ascii >= $trim);				#クオリティの閾値の判定
		shift (@quality_items);
		shift (@base_items);
		$delete_base_num++;			
	}
	
	@base_items = reverse(@base_items);
	@quality_items = reverse(@quality_items);

	$$array[3] = join ('', @quality_items);
	$$array[1] = join ('', @base_items);			#トリミングされた塩基の配列を文字列にして格納
	
	
	return($delete_base_num);
}
	
			
			
sub phase_5_N_check{

	my ($state_init, $removed_file, $array) = @_;

	state $del_r_num = 0;
	my $del_flag = 0;

    ## state initialize ##
	if ( $state_init != 0 ){
		$del_r_num = 0;
	}

	my $cnt_N = $$array[1] =~ tr/N/N/;		# Nの個数カウント

	if ( $cnt_N >= $N_content ){
		$del_r_num++;
		$del_flag = 1;

		if ($write == 1){
			print $removed_file "$$array[0]\n";	#excludedに出力
			print $removed_file "$$array[1]\n";
			print $removed_file "$$array[2]\n";
	    	print $removed_file "$$array[3]\n";
		}
	}
	
	return ($del_flag, $del_r_num);
}
	
	
	
sub phase_4_qp_check{
	
	my ($state_init, $base_len, $ratio_p, $q_thre, $removed_file, $array) = @_;
	
	state $del_r_num = 0;
	my $del_flag = 1;
	
	## state initialize ##
	if ( $state_init != 0 ){
		$del_r_num = 0;
	}
	
	my @base_items = split //, $$array[1];				#塩基を一文字ずつ配列に格納
	my @quality_items = split //, $$array[3];			#クオリティ値を一文字ずつ配列に格納
	
	my $i = $base_len * $ratio_p;
	my $ratio_thisread = sprintf("%d", $i);		#整数化（.以下切り捨て）int()は推奨されないらしい
    my $q_count = 0;
    
	foreach my $tmp (@quality_items){
		my $ascii = ord $tmp;						#クオリティ値をASCIIコードに変換
		if ($ascii >= $q_thre){			#クオリティの閾値の判定
			$q_count++;
		}
		if ($q_count >= $ratio_thisread){
			$del_flag = 0;
			last;
		}
	}
	
	if ($del_flag == 1){
			$del_r_num++;
			if ($write == 1){
				print $removed_file "$$array[0]\n";	#excludedに出力
				print $removed_file "$$array[1]\n";
				print $removed_file "$$array[2]\n";
				print $removed_file "$$array[3]\n";
			}
	}
	
	return ($del_flag, $del_r_num);
	
}



sub phase_3_casava_check {
	
	my ($state_init, $removed_file, $array) = @_;
	
	state $casava_flag = 0;
	state $del_r_num = 0;
	my $del_flag = 0;
	
	## state initialize ##
	if ( $state_init != 0 ){
		$casava_flag = 0;
		$del_r_num = 0;
	}
	
	if ($$array[0] =~ /^@.* [^:]*:N:[^:]*:/){	#CASAVA filter OK!
		$casava_flag = 1;
	}	
	elsif ($$array[0] =~ /^@.* [^:]*:Y:[^:]*:/){	#CASAVA filter NO!
		$casava_flag = 1;
		$del_r_num++;								#削除リード数カウント
		$del_flag = 1;

		if ($write == 1){
			print $removed_file "$$array[0]\n";	#excluded_phase3に出力
			print $removed_file "$$array[1]\n";
			print $removed_file "$$array[2]\n";
			print $removed_file "$$array[3]\n";
		}
	}
	else{
		$casava_flag = 0;
	}
	
	return ($casava_flag, $del_flag, $del_r_num);
}

		

sub FastQC {

	my ($TEMPD_NAME, $FNAME, $FNAME_1, $FNAME_2) = @_;
	
	## FastQC実行 ##
	if($SINGLE){
		if($fqc_nogroup == 1){
			system ( "fastqc -o $TEMPD_NAME -f fastq $FNAME -nogroup" );
		}
		else{
			system ( "fastqc -o $TEMPD_NAME -f fastq $FNAME" );
		}
	}
	elsif($PAIRED){
		if($fqc_nogroup == 1){
			system ( "fastqc -o $TEMPD_NAME -f fastq $FNAME_1 $FNAME_2 -nogroup" );
		}
		else{
			system ( "fastqc -o $TEMPD_NAME -f fastq $FNAME_1 $FNAME_2" );
		}
	}
}



sub phase_1_stylecheck {
	
	my ($out_log, $in, $in_1, $in_2, $out_log_valid) = @_;
	
	### validation log への出力 ###
	print $out_log_valid "######################################\n### QCleaner-RENEW.pl  version $VERSION ###\n######################################\n";
	print $out_log_valid "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!!!!!! START INPUT FILE VALIDATION !!!!!!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";
	
	print $out_log_valid "\n\n[ First 10 read ]\n";
	
	if($SINGLE){
		print $out_log_valid "\nINPUT file name\n[$IN_FILE]\n\n";
		for(my $count=0 ; $count<40 ; $count++){
			my $line = <$in>;
			print $out_log_valid "$line";
		}
		
		if($IN_FILE =~ /\.gz$/){
			close $in;
			open $in, "zcat $IN_FILE |" or die "Cannot open $IN_FILE file\n";
		}
		else{
			seek($in, 0, 0);		#ファイルポインタを先頭に
		}
	}
	elsif($PAIRED){
		print $out_log_valid "\nINPUT file name\n[$IN_1_FILE]\n\n";
		for(my $count=0 ; $count<40 ; $count++){
			my $line = <$in_1>;
			print $out_log_valid "$line";
		}
		print $out_log_valid "\nINPUT file name\n[$IN_2_FILE]\n\n";
		for(my $count=0 ; $count<40 ; $count++){
			my $line = <$in_2>;
			print $out_log_valid "$line";
		}
		
		if($IN_1_FILE =~ /\.gz$/){
			close $in_1;
			open $in_1, "zcat $IN_1_FILE |" or die "Cannot open $IN_1_FILE file\n";
		}
		else{
			seek($in_1, 0, 0);		#ファイルポインタを先頭に
		}
		
		if($IN_2_FILE =~ /\.gz$/){
			close $in_2;
			open $in_2, "zcat $IN_2_FILE |" or die "Cannot open $IN_2_FILE file\n";
		}
		else{
			seek($in_2, 0, 0);		#ファイルポインタを先頭に
		}
	}
	
	print $out_log_valid "\n\n[Validation phase]\n";
	
	## スタイルチェック処理 ##
	my (@dataset, @dataset_1, @dataset_2);
	if($SINGLE){
		@dataset = &stylecheck($in, $out_log, $IN_FILE, $out_log_valid);
		print $out_log_valid "\nRead format... $dataset[1] reads OK!\nID format... $dataset[1] reads OK!\n";
		
		if($p64 == 0){
			print $out_log_valid "\nPhred score is +33\n";
		}
		elsif($p64 == 1){
			print $out_log_valid "\nPhred score is +64\n";
		}
		
		return (@dataset);		#リファレンスの値をそのまま返す
	}
	elsif($PAIRED){
		@dataset_1 = &stylecheck($in_1, $out_log, $IN_1_FILE, $out_log_valid);
		@dataset_2 = &stylecheck($in_2, $out_log, $IN_2_FILE, $out_log_valid);
		print $out_log_valid "\nRead format... forward $dataset_1[1] and reverse $dataset_2[1] reads OK!\nID format... forward $dataset_1[1] and reverse $dataset_2[1] reads OK!\n";
		
		if($p64 == 0){
			print $out_log_valid "\nPhred score is +33\n";
		}
		elsif($p64 == 1){
			print $out_log_valid "\nPhred score is +64\n";
		}
		
		return (@dataset_1, @dataset_2);		#結合したハッシュのリファレンスを返す
	}


}



sub stylecheck{
	
	my ($in, $out_log, $input_fname, $out_log_valid) = @_;
	my ($cnt, $cnt_r, $cnt_l, $cnt_b) = (0, 0, 0, 0);
	my @array;
	my $phred64_flag_readid;
	
	while(<$in>){
		s/\r//g;
		chomp;
		$cnt_l++;					#行数カウント
		next if ( /^$/ );
		
		$array[$cnt] = $_;
		$cnt++;
		next if $cnt < 4;
		
		## ID行とオプション行の組み合わせが正しいか ##
		if ( $array[0] !~ /^@/ || $array[2] !~ /^\+/ ){			# リードに異常があったら終了
			print $out_log_valid "\n入力ファイル[$input_fname]の$cnt_l行目に異常がありますので、QCleanerを強制終了します。\n";
			die "\n入力ファイル[$input_fname]の$cnt_l行目に異常がありますので、QCleanerを強制終了します。\n";
		}
		
		## IDが取ってこれるかどうかの判定はpairedの時だけする。
		if($PAIRED){
			if ($array[0] !~ /\s/ && $array[0] !~ /\//){
				print $out_log_valid "入力ファイル[$input_fname]の$cnt_l行目のリードのID行にスペースまたは/が含まれていません。IDの抽出ができないためQCleanerを強制終了します。\n";
				die "入力ファイル[$input_fname]の$cnt_l行目のリードのID行にスペースまたは/が含まれていません。IDの抽出ができないためQCleanerを強制終了します。\n";
			}
		}
		
		
		## phred スコア自動判定 ##
		if($p64 == 0){
			my @quality_items = split //, $array[3];		#クオリティ値を一文字ずつ配列に格納
			my @quality_items_bup = @quality_items;
			foreach my $tmp (@quality_items_bup){
				my $ascii = ord $tmp;					#クオリティ値をASCIIコードの数値に変換
				if ($ascii >= 75){					#74以上があればphred 64
					$p64 = 1;
				}
			}
		}
		
		$cnt_r++;					#リード数カウント
		$cnt_b += length($array[1]); #塩基数カウント
		$cnt = 0;
		
	}
	
	## 総塩基数と総リード数まとめ
	return ($cnt_b, $cnt_r);
}
	
