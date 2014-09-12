# comparing outputs from workflow
#   focusing on primer bias

use strict;
use warnings;

my @fol = qw(HeadandNeck Johann Melanoma
  NATCH/FirstBatch NATCH/Repeats
  Pre-HIPAA/FirstBatch Pre-HIPAA/Repeats
  RAD001 RCC-Papillary TSTcomparison WellStar);

# start output file
open(OUT, ">comparePrimer2.txt");
print OUT "Primer\tTotal\tOverall\n";
open(IN, "comparePrimer.count");
my $tot = 98715066;
while (my $line = <IN>) {
  chomp $line;
  #next if ($line =~ m/Total/);
  my @spl = split("\t", $line);
  if ($spl[0] eq "Total") {
    print OUT "Total\t\t$spl[1]\n";
  } else {
    printf OUT "$line\t%.3f\n", 100*$spl[1]/$tot;
  }
}
close IN;
close OUT;

# loop through input files
my $total = 0;
my %pr;
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
    my $file = "comparePrimer2.txt";
    my $tfile = "comparePrimerTemp.txt";
    system "cp $file $tfile";
    open(IN, $file);
    open(OUT, ">$tfile");

    #print OUT "$fol[$f]/$dir[$x]\t";
    %pr = ();
    my $count = 0;
    #opendir(DIR, "$fol[$f]/$dir[$x]");

    # count good quality reads
    my $c5a = countFastq("$fol[$f]/$dir[$x]/joincomb-qt.fastq");
    my $c5b = countFastq("$fol[$f]/$dir[$x]/unjcomb-qt.fastq");
    $count += $c5a + $c5b;
    $total += $c5a + $c5b;

    # print output
    chomp (my $line = <IN>);
    print OUT "$line\t$dir[$x]\n";
    while (my $line = <IN>) {
      chomp $line;
      my @spl = split("\t", $line);
      if ($spl[0] eq "Total") {
        print OUT "$line\t$count\n";
      } elsif (exists $pr{$spl[0]}) {
        printf OUT "$line\t%.3f\n", 100*$pr{$spl[0]}/$count;
      } else {
        print OUT "$line\t0\n";
      }
    }

    close IN;
    close OUT;
    system "mv $tfile $file";
#die;
  }
}

print "total quality reads: $total\n";


sub countFastq {
  my $file = $_[0];
  #my %pr = %{$_[1]};
  open(FQ, $file) || die "Cannot open $file\n";
  my $count = 0;
  while (my $line = <FQ>) {
    next if (substr($line, 0, 1) ne "@");
    $count++;
    chomp $line;
    my @spl = split(" ", $line);
    die "$line has no primer\n" if (scalar @spl < 3);
    $pr{$spl[2]}++;
    for (my $x = 0; $x < 3; $x++) {
      $line = <FQ>;
    }
  }
  close FQ;
  return $count;
}
