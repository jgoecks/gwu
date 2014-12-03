#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Altering the SAM file using length-variant realignments.

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterSAMLengthVars.pl  <infile1>  <infile2>  <outfile1>  <outfile2>
  Required:
    <infile1>  File containing realignment info (produced by alignLengthVars4.pl)
    <infile2>  Input SAM file
    <outfile1> Output SAM file
    <outfile2> Verbose output file
);
  exit;
}

usage() if (scalar @ARGV < 4 || $ARGV[0] eq "-h");

open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(SAM, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]");
open(LOG, ">$ARGV[3]");

# load new alignments
my $maxExt = 0.80;  # maximum match of external in/del to consider --
                    #   above this value, the primer is a good enough
                    #   match to explain the length variant as an
                    #   artifact
my %tot; my %xExt;  # counting variables
my %rep;  # potential new alignments
my $total = 0;
while (my $line = <IN>) {
  chomp $line;
  next if (substr($line, 0, 1) eq '#');
  $total++;
  my @spl = split("\t", $line);
  $tot{$spl[1]}++;
  if (($spl[$#spl-1] eq "external") && ($spl[$#spl] > $maxExt)) {
    $xExt{$spl[1]}++;
    next;
  }
  my $read = shift @spl;
  $rep{$read} = join("\t", @spl);
}
close IN;
print "Total realignments: $total\n",
  "Valid realignments: ", scalar keys %rep, "\n";

# parse SAM file
my %xInv; my %xWorse; my %xFil; my %yBetter; my %yUnmap;  # counting variables
my $mapq = 255;  # mapping quality (255 = "unavailable")
my $line = <SAM>;
while ($line) {
  if (substr($line, 0, 1) eq '@') {
    print OUT $line;
    $line = <SAM>;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);

  # no new alignment: copy and move on
  if (!exists $rep{$spl[0]}) {
    print OUT "$line\n";
    $line = <SAM>;
    next;
  }

  # compare new alignment to previous (if any)
  my @div = split("\t", $rep{$spl[0]});
  delete $rep{$spl[0]};

  # determine if 5' or 3' end was clipped
  my $clip5 = 0; my $clip3 = 0;
  my $prev = ($spl[1] & 0x10 ? revComp($spl[9]) : $spl[9]);
  if ($div[5] ne $prev) {
    my $diff = (length $div[5]) - (length $prev);
    my $x;
    for ($x = 0; $x < $diff+1; $x++) {
      if (substr($div[5], $x, length $prev) eq $prev) {
        $clip5 = $x;
        $clip3 = $diff - $x;
        last;
      }
    }
    if ($x == $diff+1) {
      print "Error: Read $spl[0] does not match sequence in SAM record\n";
      print OUT "$line\n";
      $line = <SAM>;
      next;
    }
  }
  # adjust realignment accordingly  (skip)
  #if ($clip5 || $clip3) {
    #$div[3] += $clip5;  # position
    #$div[4] =~ m/^(\d+)([IMD])/g;
    #if ($2 eq 'M' && $1 > $clip5) {
    #  $div[4] = ($1-$clip5).$2.$';
    #}
  #}

  # construct new SAM record
  my $rec = "$spl[0]\t0\t$div[3]\t$div[4]\t$mapq\t$div[5]".
    "\t*\t0\t0\t$div[6]\t$div[7]";
  my $as = -5*$div[1] - 1 ; # for AS, count subs as -5 each, in/del as -1 only
  my $xn = 0;               # assume no Ns in ref
  my $xm = $div[1];         # number of subs
  my $xo = 1;               # 1 gap open, by definition
  $div[2] =~ m/^(\d+)([ID])/;
  my $xg = $1;              # length of gap
  my $nm = $xm + $xg;       # "edit distance"
  $rec .= "\tAS:i:$as\tXN:i:$xn\tXM:i:$xm\tXO:i:$xo\tXG:i:$xg\t".
    "NM:i:$nm\tMD:Z:$div[8]\tYT:Z:UU";
  my $extra = "\tOP:i:$spl[3]\tOC:Z:$spl[5]";

  # determine if new alignment is valid --
  #   scoring function based on bowtie2's default,
  #   but considering only subs ($xm)
  my $len = length $div[6];  # length of ref.
  $len -= $xg if ($2 eq 'I');
  if ($xm > (0.6*$len+0.6) / 5) {
    print OUT "$line\n";
    $xInv{$div[0]}++;
    while ($line = <SAM>) {
      my @cut = split("\t", $line);
      last if ($cut[0] ne $spl[0]);
      print OUT $line;
    }
    next;
  }

#print "$spl[0]\t", (0.6*$len+0.6)/5, "\t$xm\t$div[8]\n";
#my $wait = <STDIN>;

  # previously unmapped
  if ($spl[1] & 0x4) {
    print OUT "$rec$extra\n";
    $yUnmap{$div[0]}++;
    $line = <SAM>;
    next;
  }

  # get previous alignment score
  my $oldXM = getTag("XM", @spl);
  if ($oldXM == 1000) {
    die "Cannot find XM in $line\n";
    print OUT "$line\n";
    $line = <SAM>;
    next;
  }

  # compare alignment scores
  if ($xm < $oldXM || (($div[$#div-1] eq "external")
      && ($xm == $oldXM))) {

    # print new alignment
    print OUT "$rec$extra\n";
    $yBetter{$div[0]}++;

    # add previous as secondary alignment
    $spl[1] = ($spl[1] | 0x100);
    print OUT join("\t", @spl), "\n";
    while ($line = <SAM>) {
      my @cut = split("\t", $line);
      last if ($cut[0] ne $spl[0]);
      print OUT "$line\n";
    }

  } else {
    print OUT "$line\n";
    $xWorse{$div[0]}++;

    # add $rec as secondary alignment
    my @cut = split("\t", $rec);
    $cut[1] = 0x100;  # secondary alignment
    $rec = join("\t", @cut);
    my $pr = 0;
    while ($line = <SAM>) {
      my @try = split("\t", $line);
      last if ($try[0] ne $spl[0]);
      $oldXM = getTag("XM", @try);
      if ((!$pr) && $xm <= $oldXM) {
        print OUT "$rec\n";
        $pr = 1;
      }
      print OUT "$line\n";
    }
    if (!$pr) {
      print OUT "$rec\n";
    }

  }

}
close SAM;
close OUT;

# count reads not in SAM
foreach my $x (keys %rep) {
  my @div = split("\t", $rep{$x});
  $xFil{$div[0]}++;
}

# print log information
print LOG "Amplicon\tBetterMap\tNewMap\tInvalidExternalInDel\t",
  "InvalidNewMap\tWorseMap\tFiltered\tTotal\n";
foreach my $x (sort keys %tot) {
  print LOG "$x\t",
    (exists $yBetter{$x} ? $yBetter{$x} : 0), "\t",
    (exists $yUnmap{$x} ? $yUnmap{$x} : 0), "\t",
    (exists $xExt{$x} ? $xExt{$x} : 0), "\t",
    (exists $xInv{$x} ? $xInv{$x} : 0), "\t",
    (exists $xWorse{$x} ? $xWorse{$x} : 0), "\t",
    (exists $xFil{$x} ? $xFil{$x} : 0), "\t",
    (exists $tot{$x} ? $tot{$x} : 0), "\n";
}
close LOG;

# reverse-complement sequence
sub revComp {
  my $seq = $_[0];
  $seq =~ tr/ACGT/TGCA/;
  return reverse $seq;
}

# retrieve value from optional field of SAM
sub getTag {
  my $tag = shift @_;
  my @spl = @_;
  my $ret = 1000;
  my $x;
  for ($x = 11; $x < scalar @spl; $x++) {
    my @cut = split(':', $spl[$x]);
    if ($cut[0] eq $tag) {
      $ret = $cut[$#cut];
      last;
    }
  }
  if ($x == scalar @spl) {
    print "Error: cannot find $tag for $spl[0]\n";
  }
  return $ret;
}
