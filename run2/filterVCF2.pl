#!/usr/bin/perl

# John M. Gaspar
# Dec. 2014

# Filter a VCF file, limiting variants to specified
#   target regions.
# Version 2: considering MNPs

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterVCF.pl  <infile1>  <infile2>  <outfile>
  Required:
    <infile1>  Input VCF file
    <infile2>  BED file listing locations of primers
    <outfile>  Output VCF file
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

# open files
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(BED, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]");

# load target regions (convert to 1-based)
my %pos;
my %reg;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print "Improperly formatted line in BED file: $line\n",
      "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      print "Warning: skipping amplicon $spl[3] -- ",
        "located at chromosomes $spl[0] and $div[0]!?\n";
    } else {
      my $st; my $end;
      if ($spl[1] < $div[1]) {
        $st = $spl[2] + 1;
        $end = $div[1] + 1;
      } else {
        $st = $div[2] + 1;
        $end = $spl[1] + 1;
      }
      for (my $x = $st; $x < $end; $x++) {
        $reg{"$spl[0]\t$x"} = 1;
      }
    }
    delete $pos{$spl[3]};
  } else {
    $pos{$spl[3]} = "$spl[0]\t$spl[1]\t$spl[2]";
  }
}
close BED;

# parse vcf file, produce output
my $var = 0; my $good = 0;
while (my $line = <IN>) {
  if (substr($line, 0, 1) eq '#') {
    print OUT $line;
    next;
  }
  $var++;
  my @spl = split("\t", $line);
  my @div = split(',', $spl[4]);  # different alt. alleles

  # load CIGARs
  my @cig;
  my @cut = split(';', $spl[7]);
  for (my $x = 0; $x < scalar @cut; $x++) {
    my @brk = split('=', $cut[$x]);
    if ($brk[0] eq "CIGAR") {
      @cig = split(',', $brk[1]);
      last;
    }
  }
  die "CIGARs don't match alleles for $spl[0], $spl[1]\n"
    if (scalar @cig != scalar @div);

  # analyze each alternative allele
  my @all = ();  # indices of alleles to keep
  for (my $x = 0; $x < scalar @div; $x++) {
    my @diff = ();  # positions that are different
    my $pos = 0;
    while ($cig[$x] =~ m/(\d+)([IDMX])/g) {
      if ($2 ne 'M') {
        for (my $y = 0; $y < $1; $y++) {
          push @diff, $pos+$y;
          last if ($2 eq 'I');
        }
      } 
      $pos += $1 if ($2 ne 'I');
    }

    # check all diff. positions
    for (my $y = 0; $y < scalar @diff; $y++) {
      my $loc = $spl[1] + $diff[$y];
      if (exists $reg{"$spl[0]\t$loc"}) {
        push @all, $x;
        last;
      }
    }
  }

  # keep everything if any is good
  if (@all) {
    print OUT $line;
    $good++;
  }

}
close IN;
close OUT;

print "Variants in $ARGV[0]: $var\n",
  "Filtered variants in $ARGV[2]: $good\n";
