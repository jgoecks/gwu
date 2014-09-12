#!/bin/bash

# Executes workflow for all samples in a folder
# Version 4: includes reattaching primers and mapping

exit

# input file names
fol=WellStar  # folder name containing fastq files
suf=.fastq
shopt -s nullglob
arr=(../Winship/$fol/*.fastq)
num=${#arr[@]}
for ((i=0; i<num; i+=2))
do
  echo -e "Analyzing:\n\t${arr[i]}\n\t${arr[i+1]}"
  #str1=${arr[i]%%-*}
  #str2=${arr[i]#$str1-}
  str=${arr[i]%%_*}
  dir=${str##*/}
  file1=${arr[i]}
  file2=${arr[i+1]}
  #echo $dir
  #exit

  if test -d $dir; then
    echo "Repeated: $dir"
    dir=${dir}2
  fi
  dir=${dir}/
  mkdir $dir

  # stitch together reads
  echo "fastq-join"
  overlap=20
  fastq-join -m $overlap $file1 $file2 -o %.fastq
    # produces join.fastq, un1.fastq, un2.fastq

  # remove primers
  echo "removePrimer"
  pri=../Winship/primers/primers.txt
  len=20   # length of rev primer
  e=2      # errors (subs) allowed to fwd primer
  er=2     # errors (subs) allowed to rev primer
  join=join
  nm=${join}-nm$suf
  joinc=${join}comb
  out1=$dir$joinc$suf
  #cor=$dir${join}-pr-corr$suf
  removePrimer -i $join$suf -o $out1 -p $pri -l $dir$join.log -w $nm \
    -e $e -er $er -rq -r $len #-c $cor

  # retrieve reads whose primers weren't found
  echo "getReads"
  nm1=${join}-nm1$suf
  nm2=${join}-nm2$suf
  perl getReads.pl $nm $file1 $nm1
  perl getReads.pl $nm $file2 $nm2

  # cat with unjoined reads
  un1=un1$suf
  un2=un2$suf
  cat $un1 >> $nm1
  cat $un2 >> $nm2

  # remove primers
  echo "removePrimer"
  nm1p=nm1-pr$suf
  nm2p=nm2-pr$suf
  len=8
  e=2
  er=0
  removePrimer -i $nm1 -o $nm1p -p $pri -l ${dir}nm1.log -w ${dir}nu1-nm$suf \
    -e $e -er $er -r $len
  removePrimer -i $nm2 -o $nm2p -p $pri -l ${dir}nm2.log -w ${dir}nu2-nm$suf \
    -e $e -er $er -r $len

  # get common reads
  echo "getCommon"
  nm1pg=nm1-pr-gc$suf
  nm2pg=nm2-pr-gc$suf
  unjc=unjcomb
  out2=$dir$unjc$suf
  chim=${dir}chimera$suf
  perl getCommon.pl $nm1p $nm2p $nm1pg $nm2pg $out2 $chim

  # join common reads
  echo "fastq-join again"
  fastq-join -m $overlap $nm1pg $nm2pg -o nm%.fastq
    # produces nmjoin.fastq, nmun1.fastq, nmun2.fastq

  # combine joined with previous joined
  joinf=nmjoin$suf
  cat $joinf >> $out1

  # combine unjoined with other singletons
  cat nmun1$suf nmun2$suf >> $out2

  rm $join$suf $un1 $un2 $nm $nm1 $nm2 $nm1p $nm2p $nm1pg $nm2pg $joinf \
    nmun1$suf nmun2$suf
  #mv *$suf $dir


  ##################################################################################################
  # 2 output files in dir/:
  #   -  joincomb.fastq    joined paired reads, both primers removed
  #   -  unjcomb.fastq     unjoined reads with primer(s) removed

  #perl length.pl $out1
  #perl length.pl $out2
  ##################################################################################################

  # find reads of anomalous length
  per=0.01
  len=5
  lenFile=${dir}length.txt
  perl length2.pl $out1 $lenFile $per $len

  # reattach primers for joined reads
  echo "attach primers"
  temp=temp$suf
  perl attachPrimers.pl $out1 $pri $temp

  # quality filter
  echo "quality filter"
  out1q=$dir${joinc}-qt$suf
  out2q=$dir${unjc}-qt$suf
  len=3
  wqual=20
  tqual=30
  qualTrim -i $out1 -o $out1q -l $len -q $wqual -t $tqual
  qualTrim -i $out2 -o $out2q -l $len -q $wqual -t $tqual

  # for reads with reattached primers, remove those that failed quality filtering
  re=$dir${joinc}+pr$suf
  perl getReads.pl $out1q $temp $re
  rm $temp

  # bwa
  sam=.sam
  gen=../Winship/genomicRegions3.fasta
  echo "bwa: $re"
  bwa mem $gen $re > $dir${joinc}+pr$sam
  out1qb=$dir$joinc$sam
  out2qb=$dir$unjc$sam
  echo "bwa: $out1q"
  bwa mem -L 100,100 -E 3,3 $gen $out1q > $out1qb
  echo "bwa: $out2q"
  bwa mem -L 100,100 -E 3,3 $gen $out2q > $out2qb

  # cat sam files
  out2qbh=$dir${unjc}-h$sam
  outs=${dir}comb$sam
  sed '/^@/ d' < $out2qb > $out2qbh  # remove header
  cat $out1qb $out2qbh > $outs
  rm $out2qbh

  # make pileup
  echo -e "samtools\n  convert"
  bam=.bam
  qual=20
  outb=${dir}comb$bam
  samtools view -b -S -q $qual $outs > $outb   # exclude mapping qual < 20
  echo "  sort"
  pref=${dir}combsort
  samtools sort $outb $pref  # produces combsort.bam
  echo "  pileup"
  sout=$pref$bam
  pout=$pref.pileup
  Q=0
  d=100000
  samtools mpileup -B -Q $Q -d $d -f $gen $sout > $pout

  # convert pileup, VarScan
  echo "convert pileup"
  perl pileupProc.pl $pout ${pout}2
  mv ${pout}2 $pout
  echo "VarScan"
  qual=25
  final=${dir}combsort.cns
  java -jar ~/tools/VarScan/VarScan.v2.3.7.jar pileup2cns --min-avg-qual $qual < $pout > $final
  reg=regions2-.txt
  java -jar ~/tools/VarScan/VarScan.v2.3.7.jar limit $final --regions-file $reg --output-file ${final}2
  mv ${final}2 $final
  echo "final answer: $final"

done
