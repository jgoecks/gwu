# converts a pileup to use actual coordinates

use strict;
use warnings;

die "Need file names\n" if (scalar @ARGV < 2);
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(OUT, ">$ARGV[1]");
while (my $line = <IN>) {
  my @spl = split("\t", $line);
  my @div = split("_", $spl[0]);
  $spl[0] = $div[1];
  $spl[1] += $div[2];
  print OUT join("\t", @spl);
}
close IN;
close OUT;
