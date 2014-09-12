#!/bin/bash

# Executes workflow for all samples in a folder

exit

# input file names
fol=WellStar
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
  pr=${join}-pr$suf
  removePrimer -i $join$suf -o $pr -p $pri -l $dir$join.log -w $nm -e $e -er $er -rq -r $len

  # retrieve reads whose primers weren't found
  echo "getReads"
  nm1=${join}-nm1$suf
  nm2=${join}-nm2$suf
  perl getReads.pl $nm $file1 $nm1
  perl getReads.pl $nm $file2 $nm2

  # cat with unjoined reads
  un1=un1$suf
  un2=un2$suf
  nu1=nu1$suf
  nu2=nu2$suf
  cat $nm1 $un1 > $nu1
  cat $nm2 $un2 > $nu2

  # remove primers
  echo "removePrimer"
  nu1p=nu1-pr$suf
  nu2p=nu2-pr$suf
  len=8
  e=2
  er=0
  removePrimer -i $nu1 -o $nu1p -p $pri -l ${dir}nu1.log -w ${dir}nu1-nm$suf -e $e -er $er -r $len
  removePrimer -i $nu2 -o $nu2p -p $pri -l ${dir}nu2.log -w ${dir}nu2-nm$suf -e $e -er $er -r $len

  # get common reads
  echo "getCommon"
  nu1pg=nu1-pr-gc$suf
  nu2pg=nu2-pr-gc$suf
  out3=${dir}single$suf
  chim=${dir}chimera$suf
  perl getCommon.pl $nu1p $nu2p $nu1pg $nu2pg $out3 $chim

  # join common reads
  echo "fastq-join again"
  fastq-join -m $overlap $nu1pg $nu2pg -o nu%.fastq
    # produces nujoin.fastq, nuun1.fastq, nuun2.fastq

  # combine joined with previous joined
  joinc=joincomb
  out1=$dir$joinc$suf
  joinf=nujoin$suf
  cat $pr $joinf > $out1

  # combine unjoined with other singletons
  unjc=unjcomb
  out2=$dir$unjc$suf
  cat $out3 nuun1$suf nuun2$suf > $out2

  rm $join$suf $un1 $un2 $nm $pr $nm1 $nm2 $nu1 $nu2 $nu1p $nu2p $nu1pg $nu2pg $joinf \
    $out3 nuun1$suf nuun2$suf
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

  # quality filter
  echo "quality filter"
  out1q=$dir${joinc}-qt$suf
  out2q=$dir${unjc}-qt$suf
  len=3
  wqual=20
  tqual=30
  qualTrim -i $out1 -o $out1q -l $len -q $wqual -t $tqual
  qualTrim -i $out2 -o $out2q -l $len -q $wqual -t $tqual

  # bwa
  sam=.sam
  gen=../Winship/genomicRegions3.fasta
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
