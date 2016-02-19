# checks a fasta file for unique sequences

use strict;
use warnings;

open(IN, $ARGV[0]) || die "Cannot open\n";

$/ = ">";
my $waste = <IN>;
my $count = 0;
my $dup = 0;
my %fa;
while (my $chunk = <IN>) {
  chomp $chunk;
  $count++;
  my @spl = split("\n", $chunk);
  shift @spl;
  my $seq = join("", @spl);
  if (exists $fa{$seq}) {
    $dup++;
  } else {
    $fa{$seq} = 1;
  }
}
close IN;

print "Reads: $count\n",
  "Duplicates: $dup\n";

