# attaches primers to reads

use strict;
use warnings;

die "Usage: perl attachPrimers.pl  <fastq file>  <primer file>  <output file>\n"
  if (scalar @ARGV < 3);

open(IN, $ARGV[0]) || die "Cannot open fastq file $ARGV[0]\n";

# load primers;
my %fwd; my %rev; my %frc; my %rrc;
open(PR, $ARGV[1]) || die "Cannot open primer file $ARGV[1]\n";
while (my $line = <PR>) {
  chomp $line;
  my @spl = split(",", $line);
  $fwd{$spl[0]} = $spl[1];
  $rev{$spl[0]} = $spl[2];
  #my $fr = reverse $spl[1];  # for rc reads -- but everything should be fwd
  #$fr =~ tr/ACGT/TGCA/;
  #$frc{$spl[0]} = $fr;
  #my $rr = reverse $spl[2];
  #$rr =~ tr/ACGT/TGCA/;
  #$rrc{$spl[0]} = $rr;
}
close PR;

open(OUT, ">$ARGV[2]");
while (my $line = <IN>) {
  my @spl = split(" ", $line);
  if (! exists $fwd{$spl[2]}) {
    die "Unknown primer: $spl[2]\n";
  }
  my $seq = <IN>;
  chomp $seq;
  my $plus = <IN>;
  my $qual = <IN>;
  chomp $qual;
  print OUT "$line$fwd{$spl[2]}$seq$rev{$spl[2]}\n$plus",
    "I" x length $fwd{$spl[2]}, $qual,
    "I" x length $rev{$spl[2]}, "\n";
}
close IN;
close OUT;
