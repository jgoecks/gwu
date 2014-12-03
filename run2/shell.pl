# runs the whole friggin thing

use strict;
use warnings;

#my @fol = qw(HeadandNeck);   ###
#my @fol = qw(Johann NATCH/FirstBatch NATCH/Repeats
#  Pre-HIPAA/FirstBatch Pre-HIPAA/Repeats
#  RAD001 TSTcomparison WellStar);
my @fol = qw(NATCH/FirstBatch NATCH/Repeats
  Pre-HIPAA/FirstBatch Pre-HIPAA/Repeats);
system "mkdir NATCH";
system "mkdir Pre-HIPAA";
my $gen = "/home/john/Downloads/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
my $idx = "/home/john/Downloads/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome";
my $bed = "primers.bed";

open(LOG, ">log2.txt");
my $time = time;

# first step: retrieve primer-target sequences
#my $prim = "../primers.txt";   ###
my $prim = "primers.txt";
system "perl getPrimers.pl $bed $gen $prim";
print LOG "getPrimers.pl: ", time - $time,
  " seconds\n";
$time = time;

# loop through folders
my $pre = "../../Winship";
for (my $f = 0; $f < scalar @fol; $f++) {
  opendir(DIR, "$pre/$fol[$f]");
  my @file = grep {/\.fastq/} readdir DIR;
  closedir DIR;
  @file = sort @file;
  system "mkdir $fol[$f]";

  for (my $x = 0; $x < scalar @file; $x += 2) {
    my $file1 = "$pre/$fol[$f]/$file[$x]";
    my $file2 = "$pre/$fol[$f]/$file[$x+1]";

    # check files for errors
    my @cut = split('_', $file[$x]);
    if ($cut[$#cut-1] ne "R1") {
      print "Error with $fol[$f]/$cut[0]\n";
      next;
    }
    $cut[$#cut-1] = "R2";
    if ($file[$x+1] ne join('_', @cut)) {
      print "Error with $fol[$f]/$cut[0]\n";
      next;
    }

    # make directory for output files
    while (-d "$fol[$f]/$cut[0]") {
      $cut[0] .= '-';
    }
    my $dr = "$fol[$f]/$cut[0]";
    system "mkdir $dr";

    # join
    my $tr1 = "join.fastq";
    my $tr2 = "un1.fastq";
    my $tr3 = "un2.fastq";
    system "fastq-join -m 20 $file1 $file2 -o %.fastq";

    # removePrimer with -rq
    my $out1 = "joined.fastq";
    my $log1 = "joinlog.txt";
    my $tr4 = "join-nopr.fastq";
    my $rpParam = "-fp -1,1 -rp -1,1 -ef 2 -er 3 -rl 16 -el 1 -b $bed -bp -1,1";
    system "./removePrimer6 -i $tr1 -p $prim -o $out1 $rpParam -rq -l $log1 -w $tr4";

    # get failures
    my $tr5 = "nopr1.fastq";
    my $tr6 = "nopr2.fastq";
    system "perl getReads.pl $tr4 $file1 $tr5";
    system "perl getReads.pl $tr4 $file2 $tr6";

    # cat failures with unjoined
    system "cat $tr2 >> $tr5";
    system "cat $tr3 >> $tr6";    

    # removePrimer individually
    my $tr7 = "nopr1-pr.fastq";
    my $tr8 = "nopr2-pr.fastq";
    my $log2 = "nopr1log.txt";
    my $log3 = "nopr2log.txt";
    system "./removePrimer6 -i $tr5 -p $prim -o $tr7 $rpParam -l $log2";
    system "./removePrimer6 -i $tr6 -p $prim -o $tr8 $rpParam -l $log3";

    # get common reads and singletons
    my $tr9 = "gc1.fastq";
    my $tr10 = "gc2.fastq";
    my $out2 = "single.fastq";
    my $tr11 = "chim.fastq";
    system "perl getCommon.pl $tr7 $tr8 $tr9 $tr10 $out2 $tr11";

    # join matches
    my $tr12 = "gcjoin.fastq";
    my $tr13 = "gcun1.fastq";
    my $tr14 = "gcun2.fastq";
    system "fastq-join -m 20 $tr9 $tr10 -o gc%.fastq";

    # cat joined together, and singletons together
    system "cat $tr12 >> $out1";
    system "cat $tr13 $tr14 >> $out2";

    # cat joined and singletons
    my $out3 = "comb.fastq";
    system "cat $out1 $out2 > $out3";

    # quality trim
    my $out4 = "comb-qt.fastq";
    my $qtParam = "-l 3 -q 20 -t 30";
    system "qualTrim -i $out3 -o $out4 $qtParam";

    # map
    my $out5 = "comb-qt.sam";
    my $bwtParam = "-D 200 -N 1 -L 18 -i S,1,0.50 -k 20";
    system "bowtie2 -x $idx -U $out4 -S $out5 $bwtParam -p 4";

    # find length variants
if (0) {
    my $len1 = "len1.txt";
    #my $flvParam = "0.01 3";  # 1%, 3bp diff
    system "perl ../findLengthVars.pl $out1 $bed $len1"; # $flvParam";
    my $len2 = "len2.fastq";
    system "perl ../getLengthVars.pl $out1 $len1 $len2";
    my $len3 = "len3.txt";
    my $len3v = "len3v.txt";
    system "perl ../alignLengthVars5.pl $len2 $prim $bed $len3 $len3v";
    my $out6 = "comb-qtFil.sam";
    my $len4 = "len4.txt";
    system "perl ../filterSAMLengthVars.pl $len3 $out5 $out6 $len4";
}

    # remove extra files
    unlink $tr1, $tr2, $tr3, $tr4, $tr5, $tr6, $tr7, $tr8,
      $tr9, $tr10, $tr11, $tr12, $tr13, $tr14;

    # move files to directory
    system "mv $out1 $out2 $out3 $out4 $out5 $log1 $log2 $log3 $dr"; #$out6

    print LOG "$fol[$f]/$file[$x]: ", time - $time,
      " seconds\n";
    $time = time;
  }
}
