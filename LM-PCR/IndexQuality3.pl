### �ȉ��̑���Ń}�g���b�N�X�I�Ȃ��̂����B ###########################################################
### �s�ɌX�̃��[�h������B ###########################################################################
### �J�����͍�����AR1��ID�A�z��A�R�ڂ̂���(+)�A�N�I���e�B�l�AR2��ID�A�A�AI1��ID�A�A�A�̏��ɕ��ԁB###

###### R1                             ######
###### �R�}���h��������t�@�C�����w�� ######

open (IN1, $ARGV[0]); #1�Ԗڂ̕W������(R1)

#@data�ɗv�f����� [0][1][2][3]
#�������Ƀ��[�hID�A�V�[�N�G���X�z��A�R�ڂ̂���(+)�A�N�I���e�B�l�����ɕ��сA
#�c�����Ɋe���[�h�����ԁB

$i=0;
$j=0;

while ($line=<IN1>){

   $data[$i][$j] = $line;

   if($j<3){$j=$j+1;}else{$i=$i+1;$j=0;}

}

###### ��������� R2 �ɂ��s��         ######
###### �R�}���h��������t�@�C�����w�� ######

open (IN2, $ARGV[1]); #R2

#�ȉ��̑���� R1 �̉��� R2 ����ׂ�B[4][5][6][7]

$i=0;
$j=4;

while ($line=<IN2>){

   $data[$i][$j] = $line;

   if($j<7){$j=$j+1;}else{$i=$i+1;$j=4;}

}

###### ���x�̓C���f�b�N�X��Fastq�t�@�C���̃f�[�^��R2�̉��ɕ��ׂ�B######
###### �R�}���h��������t�@�C�����w��                             ######

open (IN3, $ARGV[2]); #I1

#R2�̉���I1����ׂ�B[8][9][10][11]

$i=0;
$j=8;

while ($line=<IN3>){

   $data[$i][$j] = $line;

   if($j<11){$j=$j+1;}else{$i=$i+1;$j=8;}

}

########################################################################################################



######���x�́A�W�����͂ŁA�o�̓t�@�C�����w�肵�Ă����B ######

open (OUT1, ">$ARGV[3]");   ###�o�̓t�@�C�� for R1
open (OUT2, ">$ARGV[4]");   ###�o�̓t�@�C�� for R2
open (OUT3, ">$ARGV[5]");   ###�o�̓t�@�C�� for I1

#############################################################



for ($i=0 ; $i<@data ; $i++){				###��s���A���Ă����B


		@quality=split(//, $data[$i][11]);	###�C���f�b�N�X�̃N�I���e�B�l��1�������ɕ����B

		$number[0]=ord($quality[0])-33;		###�N�I���e�B�l�̕������A�����ɒ����B
		$number[1]=ord($quality[1])-33;
		$number[2]=ord($quality[2])-33;
		$number[3]=ord($quality[3])-33;
		$number[4]=ord($quality[4])-33;
		$number[5]=ord($quality[5])-33;
		$number[6]=ord($quality[6])-33;
		$number[7]=ord($quality[7])-33;

			if($number[0]>30 and $number[1]>30 and $number[2]>30 and $number[3]>30 and $number[4]>30 and $number[5]>30 and $number[6]>30 and $number[7]>30){	###�N�I���e�B�l�̐ݒ蕔���I�I�I�I�I�I���̏ꍇ�́A8����ׂ�30��荂����΁B


        			print OUT1 $data[$i][0];	###R1���o��
        			print OUT1 $data[$i][1];
       				print OUT1 $data[$i][2];
       				print OUT1 $data[$i][3];

     				print OUT2 $data[$i][4];	###R2���o��
      				print OUT2 $data[$i][5];
    				print OUT2 $data[$i][6];
      				print OUT2 $data[$i][7];

        			print OUT3 $data[$i][8];	###�ꉞ�A�C���f�b�N�X��fastq���o�́B
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
