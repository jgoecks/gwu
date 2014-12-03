#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Produce a BED file listing combined target regions
#   (sorted) from a BED file listing primer positions.

use strict;
use warnings;

sub usage {
  print q(Usage: perl getTargetBED.pl  <infile>  <outfile>
  Required:
    <infile>   Input BED file
    <outfile>  Output BED file
);
  exit;
}

usage() if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

# open files
open(BED, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(OUT, ">$ARGV[1]");

# load primer locations: chr# - 5'Loc - 5'PLen - targLen - 3'PLen
my %pos;  # position of amplicons
my %loc;  # for sorted amplicon names
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
      delete $pos{$spl[3]};
    }
    if ($spl[1] < $div[1]) {
      $pos{$spl[3]} = "$spl[0]\t$spl[1]\t".($spl[2]-$spl[1]).
        "\t".($div[1]-$spl[2])."\t$div[2]";
    } else {
      $pos{$spl[3]} .= "\t".($spl[1]-$div[1]-$div[2]).
        "\t".($spl[2]-$spl[1]);
    }

    # save amplicon name, sorted
    if ($spl[0] !~ m/chr((\d+|[MXY]))/i) {
      print "Invalid chromosome name $spl[0]\n";
      next;
    }
    $loc{$1}{$spl[1]} = $spl[3];
  } else {
    $pos{$spl[3]} = "$spl[0]\t$spl[1]\t".($spl[2]-$spl[1]);
  }
}
close BED;

# produce output BED, sorted
if (exists $loc{"chrM"}) {
  my $prev = 0;
  foreach my $pl (sort {$a <=> $b} keys %{$loc{"chrM"}}) {
    my $amp = $loc{"chrM"}{$pl};
    my @spl = split("\t", $pos{$amp});
    if (!$prev) {
      print OUT "$spl[0]\t", $spl[1]+$spl[2];
      $prev = $spl[1]+$spl[2]+$spl[3];
    } elsif ($spl[1]+$spl[2] <= $prev) {
      $prev = $spl[1]+$spl[2]+$spl[3];
    } else {
      print OUT "\t$prev\n$spl[0]\t", $spl[1]+$spl[2];
      $prev = $spl[1]+$spl[2]+$spl[3];
    }
  }
  print OUT "\t$prev\n";

  delete $loc{"chrM"};
}

foreach my $chr (sort {$a <=> $b} keys %loc) {
  my $prev = 0;
  foreach my $pl (sort {$a <=> $b} keys %{$loc{$chr}}) {
    my $amp = $loc{$chr}{$pl};
    my @spl = split("\t", $pos{$amp});
    if (!$prev) {
      print OUT "$spl[0]\t", $spl[1]+$spl[2];
      $prev = $spl[1]+$spl[2]+$spl[3];
    } elsif ($spl[1]+$spl[2] <= $prev) {
      $prev = $spl[1]+$spl[2]+$spl[3];
    } else {
      print OUT "\t$prev\n$spl[0]\t", $spl[1]+$spl[2];
      $prev = $spl[1]+$spl[2]+$spl[3];
    }
  }
  print OUT "\t$prev\n";
}
close OUT;
