# comparing outputs from workflow

use strict;
use warnings;

my @fol = qw(HeadandNeck Johann Melanoma
  NATCH/FirstBatch NATCH/Repeats
  Pre-HIPAA/FirstBatch Pre-HIPAA/Repeats
  RAD001 RCC-Papillary TSTcomparison WellStar);

my %pr;
open(PLOG, ">comparePrimer5.count");
open(OUT, ">compareOutputs5.txt");
print OUT "File\tInput Reads\tJoined\tAdded\tSingletons\tChimeras",
  "\tQual.removed\tNo primer\tGood mapped\tSingletons mapped\n";
my $total = 0;
for (my $f = 0; $f < scalar @fol; $f++) {
  opendir(DIR, $fol[$f]);
  my @dir = readdir DIR;
  closedir DIR;
  @dir = sort @dir;
  #print @dir, "\n";

  opendir(D2, "../Winship/$fol[$f]");
  my @d2 = grep {/\.fastq/} readdir D2;
  closedir D2;
  @d2 = sort @d2;

  for (my $x = 2; $x < scalar @dir; $x++) {
    print OUT "$fol[$f]/$dir[$x]\t";
    #opendir(DIR, "$fol[$f]/$dir[$x]");

    # get count of original reads
    my $pos = 2 * ($x - 2);
    my $file = $d2[$pos];
    my $f2 = $d2[$pos+1];
    my $count = countFastq("../Winship/$fol[$f]/$file", 0);
    my $count2 = countFastq("../Winship/$fol[$f]/$f2", 0);
    die "Counts $count and $count2 in $file and $f2 don't match\n" if ($count != $count2);
    $file =~ s/_R1_/_R2_/;
    die "File names $file and $f2 don't match\n" if ($file ne $f2);
    $count += $count2;
    print OUT "$count\t";

    # get count of good joined
    my $log = "$fol[$f]/$dir[$x]/join.log";
    open(LOG, $log) || die "Cannot open $log\n";
    my $line;
    for (my $y = 0; $y < 4; $y++) {
      $line = <LOG>;
    }
    chomp $line;
    my @spl = split(": ", $line);
    printf OUT "%.1f\t", 200*$spl[1]/$count;
    #print OUT "$spl[1]\t";

    # get count of bonus joined
    my $c2 = countFastq("$fol[$f]/$dir[$x]/joincomb.fastq", 0);
    printf OUT "%.1f\t", 200*($c2-$spl[1])/$count;
    #print OUT $c2 - $spl[1], "\t";

    # get count of singletons
    my $c3 = countFastq("$fol[$f]/$dir[$x]/unjcomb.fastq", 0);
    printf OUT "%.1f\t", 100*$c3/$count;
    #print OUT "$c3\t";

    # count chimeras
    my $c4 = countFastq("$fol[$f]/$dir[$x]/chimera.fastq", 0);
    printf OUT "%.1f\t", 100*$c4/$count;
    #print OUT "$c4\t";

    # count good quality reads
    my $c5a = countFastq("$fol[$f]/$dir[$x]/joincomb-qt.fastq", 1);
    my $c5b = countFastq("$fol[$f]/$dir[$x]/unjcomb-qt.fastq", 1);
    my $c5 = 2*($c2 - $c5a) + $c3 - $c5b;
$total += $c5a + $c5b;

    printf OUT "%.1f\t", 100*$c5/$count;
    #print OUT "$c5\t";

    # count no primer found
    my $c6 = countFastq("$fol[$f]/$dir[$x]/nu1-nm.fastq", 0);
    $c6 += countFastq("$fol[$f]/$dir[$x]/nu2-nm.fastq", 0);
    printf OUT "%.1f\t", 100*$c6/$count;
    #print OUT "$c6\t";

    # parse SAM files
    my $map = parseSam("$fol[$f]/$dir[$x]/joincomb.sam", $c5a);
    printf OUT "%.3f\t", 100*$map/$c5a;
    my $mapu = parseSam("$fol[$f]/$dir[$x]/unjcomb.sam", $c5b);
    printf OUT "%.1f\n", 100*$mapu/$c5b;
    #print OUT "\n";
  }
}
close OUT;

print "total quality reads: $total\n";
die;

# print primer counts
open(IN, "../Winship/primers/primers.txt") || die "Cannot open primers.txt\n";
my $tot = 0;
while (my $line = <IN>) {
  my @spl = split(",", $line);
  if (exists $pr{$spl[0]}) {
    print PLOG "$spl[0]\t$pr{$spl[0]}\n";
    $tot += $pr{$spl[0]};
    delete $pr{$spl[0]};
  } else {
    print PLOG "$spl[0]\t0\n";
  }
}
print PLOG "Total: $tot\n";
close IN;
close PLOG;

foreach my $key (keys %pr) {
  print "Error! $key not found\n";
}


sub countFastq {
  my $file = shift @_;
  my $get = shift @_;
  open(FQ, $file) || die "Cannot open $file\n";
  my $count = 0;
  while (my $line = <FQ>) {
    next if (substr($line, 0, 1) ne "@");
    $count++;
    if ($get) {
      chomp $line;
      my @spl = split(" ", $line);
      die "$line has no primer\n" if (scalar @spl < 3);
      $pr{$spl[2]}++;
    }
    for (my $x = 0; $x < 3; $x++) {
      $line = <FQ>;
    }
  }
  close FQ;
  return $count;
}

sub parseSam {
  my $file = shift @_;
  my $reads = shift @_;
  open(SAM, $file) || die "Cannot open $file\n";
  my $count = 0;
  my $map = 0;
  while (my $line = <SAM>) {
    next if (substr($line, 0, 1) eq "@");
    my @spl = split("\t", $line);
    next if ($spl[1] & 0x800);
    $count++;
    $map++ if (!($spl[1] & 0x4));
  }
  die "Read count $count in $file does not match $reads\n" if ($count != $reads);
  close SAM;
  return $map;
}
