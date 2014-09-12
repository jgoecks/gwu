# lengths of a fastq file
#   output by amplicon

use strict;
use warnings;

# usage: perl length2.pl  <input fastq>  <output txt>  <percent>  <distance>
#    where <percent> is the fraction of the total reads to consider (e.g. 0.01, which is 1%)
#    and <distance> is the difference in length from the maximum to consider (e.g. 5 bp)

die "Need file name\n" if (scalar @ARGV < 1);
open(IN, $ARGV[0]) || die "cannot open\n";
my $total = 0;
my $count = 0;
my %len = ();

while (my $line1 = <IN>) {
  next if (substr($line1, 0, 1) ne "@");
  chomp $line1;
  my @spl = split(" ", $line1);
  my $line = <IN>;
  chomp $line;
  $count++;
  my $ln = length $line;
  $total += $ln;
  $len{$spl[2]}{$ln}++;
  $line = <IN>;
  die "Something ain't right\n" if ($line ne "+\n");
  $line = <IN>;
}
close IN;

print "Reads: $count\n";
printf "Length: %.1f\n", $total / $count;

my $pct = 0.01;  # percent of total to include
my $dist = 5;    # distance (bp) from max to include
if (scalar @ARGV > 1) {
  if (scalar @ARGV > 3) {
    $pct = $ARGV[2];
    $dist = $ARGV[3];
  }
  open(OUT, ">$ARGV[1]");
  print OUT "Amplicon\tLength\tPercent\tMode Length\n";
  open(PR, "../Winship/primers/primers.txt") || die "Cannot open primers.txt\n";
  while (my $line = <PR>) {
    my @spl = split(",", $line);
    if (exists $len{$spl[0]}) {
      my $max = 0;
      my $maxlen = 0;
      my $total = 0;
      foreach my $key (keys %{$len{$spl[0]}}) {
        if ($len{$spl[0]}{$key} > $max) {
          $max = $len{$spl[0]}{$key};
          $maxlen = $key;
        }
        $total += $len{$spl[0]}{$key};
      }
      print OUT "$spl[0]\t$total\n";
      for (my $x = 0; $x < 150; $x++) {
        if (exists $len{$spl[0]}{$x}) {
          if ($x == $maxlen) {
            print OUT "\t$x\t$len{$spl[0]}{$x}\n";
          } elsif (($len{$spl[0]}{$x} > $pct * $total) &&
              (($x < $maxlen - $dist) || ($x > $maxlen + $dist))) {
            printf OUT "\t$x\t%.1f\t$maxlen\n", 100*$len{$spl[0]}{$x}/$total;
            #printf OUT "$spl[0]\t$x\t%.1f\t$maxlen\n", 100*$len{$spl[0]}{$x}/$total;
          }
        }
      }
    }
  }
  close PR;
  close OUT;
}
