# retrieve reads without primer matches

use strict;
use warnings;

my %un;
my $count = 0;
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
while (my $line = <IN>) {
  chomp $line;
  next if (substr($line, 0, 1) ne "@");

  my @spl = split(" ", $line);
  $un{$spl[0]} = 1;

  for (my $x = 0; $x < 3; $x++) {
    my $waste = <IN>;
  }
  $count++;
}
close IN;

print "Reads: $count\n";

open(FQ, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(FQO, ">$ARGV[2]");
my $print = 0;
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne "@");
  my @spl = split(" ", $line);
  my $q = <FQ>;
  my $w = <FQ>;
  my $e = <FQ>;
  if (exists $un{$spl[0]}) {
    print FQO "$line$q$w$e";
    $print++;
  }
}
print "Printed: $print\n";
