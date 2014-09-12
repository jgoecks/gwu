# Converts SAM file to fastq

use strict;
use warnings;

die "Need file name\n" if (scalar @ARGV < 2);
open(IN, $ARGV[0]) || die "Cannot open\n";
open(OUT, ">$ARGV[1]");

my $count = 0;
my $yes = 0;
my $no = 0;
while (my $line = <IN>) {
  chomp $line;
  next if (substr($line, 0, 1) eq "@");
  my @spl = split("\t", $line);
  if ($spl[1] & 0x800) {
    next;
  }

  if ($spl[1] & 0x4) {
    $no++;
  } else {
    $count++;
    my $st = 0; 
    my $end = length $spl[9];
    if ($spl[5] =~ m/^(\d+)S/) {  # change S to [SI] to include insertions
      $st = $1;
    }
    if ($spl[5] =~ m/(\d+)S$/) {
      $end -= $1;
    }
#if ($spl[5] =~ m/S/) {
#  print "$spl[9]\n$spl[5]\n", substr($spl[9], $st, $end-$st), "\n";
#  my $wait = <STDIN>;
#}
    print OUT '@', "$spl[0] $spl[2] $spl[3] $spl[4]\n",
      substr($spl[9], $st, $end-$st), "\n+\n",
      substr($spl[10], $st, $end-$st), "\n";
  }

}
close IN;
close OUT;

print "Reads printed: $count\n";
print "Not mapped: $no",
  "\n";
