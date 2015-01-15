#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Find putative small structural variants (based on
#   length variants) in a fastq file.

use strict;
use warnings;

sub usage {
  print q(Usage: perl getLengthVars.pl  <infile1>  <infile2>  <outfile>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer)
    <infile2>  File containing length variants (produced by findLengthVars.pl)
    <outfile>  Output file for reads
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

# load length variants
my %var;  # variants
my @ord;  # order of amplicons
open(AMP, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
while (my $line = <AMP>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split("\t", $line);  # $spl[0] is name, $spl[2] is length
  die "$ARGV[1] is improperly formatted\n" if (scalar @spl < 3);
  $var{$spl[0]}{$spl[2]} = $spl[2] - $spl[1]; # $spl[2]-$spl[1] is difference
  push @ord, $spl[0];
}
close AMP;

# parse fastq file
my $count = 0;
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(OUT, ">$ARGV[2]");
while (my $line = <IN>) {
  next if (substr($line, 0, 1) ne "@");
  chomp $line;
  my @spl = split(" ", $line);  # $spl[2] has amplicon name
  my $line2 = <IN>;
  chomp $line2;
  my $ln = length $line2;
  my $flag = 0;
  if (exists $var{$spl[2]} && exists $var{$spl[2]}{$ln}) {
    print OUT "$line ", abs $var{$spl[2]}{$ln},
      $var{$spl[2]}{$ln} < 0 ? "D" : "I", "\n$line2\n";
    $count++;
    $flag = 1;
  }
  for (my $x = 0; $x < 2; $x++) {
    $line = <IN>;
    print OUT $line if ($flag);
  }
}
close IN;

print "Reads printed to $ARGV[2]: $count\n";
