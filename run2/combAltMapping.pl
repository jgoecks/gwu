# combine altMapping results for the samples --
#   no need to repeat all that stuff

use strict;
use warnings;

my @fol = qw(HeadandNeck Johann
  NATCH/FirstBatch NATCH/Repeats
  Pre-HIPAA/FirstBatch Pre-HIPAA/Repeats
  RAD001 TSTcomparison WellStar);

# loop through folders
my %tot;
for (my $f = 0; $f < scalar @fol; $f++) {
  opendir(DIR, $fol[$f]);
  my @dir = readdir DIR;
  closedir DIR;
  @dir = sort @dir;

  for (my $d = 2; $d < scalar @dir; $d++) {

    my $dr = "$fol[$f]/$dir[$d]";
    print "Analyzing $dr\n";

    open(IN, "$dr/altMapping.txt");
    while (my $line = <IN>) {
      next if (substr($line, 0, 1) eq '#');
      chomp $line;
      my @spl = split("\t", $line);
      my @div = split('-', $spl[4]);
      for (my $x = $div[0]; $x < $div[$#div]+1; $x++) {
        if (exists $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x}) {
          if ($spl[5] != $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x}) {
            #print "Warning: $dr has $spl[5] for $pre, previously $res{$pre}{$x}\n";
            $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x} = 1;
          }
        } else {
if ($spl[0] eq "AKT1_2") { print "$dr: $line\n"; }
          $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x} = $spl[5];
        }
      }
    }
  }
}

open(OUT, ">altMapping3.txt");
print OUT "#Amplicon\tPrimersRemoved\tStrand\tChrom\tPosition(s)\tMatch?\n";
foreach my $am (sort {$a cmp $b} keys %tot) {
  foreach my $bo (sort {$b cmp $a} keys %{$tot{$am}}) {
    foreach my $st (sort {$a cmp $b} keys %{$tot{$am}{$bo}}) {
      foreach my $ch (sort {$a cmp $b} keys %{$tot{$am}{$bo}{$st}}) {

        my $res; my $min;
        my $prev = -10;  # previous position analyzed
        my @loc = sort {$a <=> $b} keys %{$tot{$am}{$bo}{$st}{$ch}};
        for (my $x = 0; $x < scalar @loc; $x++) {
          # combine results for neighboring positions --
          #   consider a match if any is a match
          if ($loc[$x] < $prev + 4) {
            $res = 1 if ($tot{$am}{$bo}{$st}{$ch}{$loc[$x]});
          } else {
            if ($x) {
              printf OUT "$am\t$bo\t$st\t$ch\t%s\t$res\n",
                ($min != $prev ? "$min-$prev" : $min);
            }
            $res = $tot{$am}{$bo}{$st}{$ch}{$loc[$x]};
            $min = $loc[$x];
          }
          $prev = $loc[$x];
        }
        if (@loc) {
          printf OUT "$am\t$bo\t$st\t$ch\t%s\t$res\n",
            ($min != $prev ? "$min-$prev" : $min);
        }
      }
    }
  }
}
