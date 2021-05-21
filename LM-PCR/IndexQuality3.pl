### 以下の操作でマトリックス的なものを作る。 ###########################################################
### 行に個々のリードが入る。 ###########################################################################
### カラムは左から、R1のID、配列、３つ目のもの(+)、クオリティ値、R2のID、、、I1のID、、、の順に並ぶ。###

###### R1                             ######
###### コマンド引数からファイルを指定 ######

open (IN1, $ARGV[0]); #1番目の標準入力(R1)

#@dataに要素を入力 [0][1][2][3]
#横方向にリードID、シークエンス配列、３つ目のもの(+)、クオリティ値が順に並び、
#縦方向に各リードが並ぶ。

$i=0;
$j=0;

while ($line=<IN1>){

   $data[$i][$j] = $line;

   if($j<3){$j=$j+1;}else{$i=$i+1;$j=0;}

}

###### 同じ操作を R2 にも行う         ######
###### コマンド引数からファイルを指定 ######

open (IN2, $ARGV[1]); #R2

#以下の操作で R1 の横に R2 を並べる。[4][5][6][7]

$i=0;
$j=4;

while ($line=<IN2>){

   $data[$i][$j] = $line;

   if($j<7){$j=$j+1;}else{$i=$i+1;$j=4;}

}

###### 今度はインデックスのFastqファイルのデータをR2の横に並べる。######
###### コマンド引数からファイルを指定                             ######

open (IN3, $ARGV[2]); #I1

#R2の横にI1を並べる。[8][9][10][11]

$i=0;
$j=8;

while ($line=<IN3>){

   $data[$i][$j] = $line;

   if($j<11){$j=$j+1;}else{$i=$i+1;$j=8;}

}

########################################################################################################



######今度は、標準入力で、出力ファイルを指定しておく。 ######

open (OUT1, ">$ARGV[3]");   ###出力ファイル for R1
open (OUT2, ">$ARGV[4]");   ###出力ファイル for R2
open (OUT3, ">$ARGV[5]");   ###出力ファイル for I1

#############################################################



for ($i=0 ; $i<@data ; $i++){				###一行ずつ、見ていく。


		@quality=split(//, $data[$i][11]);	###インデックスのクオリティ値を1文字ずつに分解。

		$number[0]=ord($quality[0])-33;		###クオリティ値の文字を、数字に直す。
		$number[1]=ord($quality[1])-33;
		$number[2]=ord($quality[2])-33;
		$number[3]=ord($quality[3])-33;
		$number[4]=ord($quality[4])-33;
		$number[5]=ord($quality[5])-33;
		$number[6]=ord($quality[6])-33;
		$number[7]=ord($quality[7])-33;

			if($number[0]>30 and $number[1]>30 and $number[2]>30 and $number[3]>30 and $number[4]>30 and $number[5]>30 and $number[6]>30 and $number[7]>30){	###クオリティ値の設定部分！！！！！！この場合は、8塩基すべて30より高ければ。


        			print OUT1 $data[$i][0];	###R1を出力
        			print OUT1 $data[$i][1];
       				print OUT1 $data[$i][2];
       				print OUT1 $data[$i][3];

     				print OUT2 $data[$i][4];	###R2を出力
      				print OUT2 $data[$i][5];
    				print OUT2 $data[$i][6];
      				print OUT2 $data[$i][7];

        			print OUT3 $data[$i][8];	###一応、インデックスのfastqも出力。
        			print OUT3 $data[$i][9];
        			print OUT3 $data[$i][10];
         			print OUT3 $data[$i][11];


	}

}

close (IN1);
close (IN2);
close (IN3);
close (OUT1);
close (OUT2);
close (OUT3);
