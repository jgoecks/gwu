#!/usr/bin/perl

# John M. Gaspar
# Jan. 2015

# Find variants that are adjacent to homopolymers.

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterVars.pl  <infile>  <genome>  <outfile>  <logfile>
  Required:
    <infile>   Input VCF file listing variants, freebayes-style
    <genome>   Fasta file of reference genome
    <outfile>  Output VCF file listing good variants
    <logfile>  Verbose log file
  Optional:
    <length>   Homopolymer length adjacent to variant (default: 5)
    <abunSub>  Minimum abundance to keep a substitution next to a homopolymer (default: 0.20)
    <abunI/D>  Minimum abundance to keep an in/del next to a homopolymer (default: 0.40)
);
  exit;
}

usage() if (scalar @ARGV < 4 || $ARGV[0] eq "-h");

# open files
open(VCF, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(GEN, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
$/ = '>';
my $waste = <GEN>;
$/ = "\n";
open(OUT, ">$ARGV[2]");
open(LOG, ">$ARGV[3]");
print LOG "#CHROM\tPOS\tREF\tALT\tCIGAR\tGenomeSeg\n";

# get optional parameters
my $hom = 4;
$hom = $ARGV[4]-1 if (scalar @ARGV > 4);
my $abSub = 0.2;
$abSub = $ARGV[5] if (scalar @ARGV > 5);
my $abID = 0.4;
$abID = $ARGV[6] if (scalar @ARGV > 6);

# start analyzing VCF file
my $chr = "";
my $seq = "";
my $rem = 0; my $kept = 0;
while (my $line = <VCF>) {
  if (substr($line, 0, 1) eq '#') {
    if (substr($line, 1, 1) ne '#') {
      print OUT "##FORMAT=<ID=HP,Number=1,Type=Integer,",
        "Description=\"Length of adjacent homopolymer\">\n";
    }
    print OUT $line;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);
  $spl[8] .= ":HP";

  # determine abundance
  $spl[7] =~ m/AB\=(.*?)\;/;
  my $ab = $1;
  if ($ab == 0) {
    $spl[7] =~ m/AO\=(.*?)\;/;
    my $ct = $1;
    $spl[7] =~ m/DP\=(.*?)\;/;
    $ab = $ct / $1;
  }

  # determine position of variant(s) using CIGAR
  $spl[7] =~ m/CIGAR\=(.*?)\;/;
  my $cig = $1;
  my @loc = ();
  my $pos = 0; my $sub = 0;
  my $alt = '*';
  my $prev = "";
  my $base = "";  # different nt
  my $flag = 0;
  while ($cig =~ m/(\d+)([IDMX])/g) {
    if ($2 ne 'M') {
      if ($prev) {
        $alt .= " $prev *";
        $prev = "";
      }
      my $diff = "";
      for (my $x = 0; $x < $1; $x++) {
        push @loc, $spl[1] + $pos + $x;
        if ($2 eq 'D') {
          if ($base) {
            $flag = 1 if ($base ne substr($spl[3], $pos+$x, 1));
          } else {
            $base = substr($spl[3], $pos+$x, 1);
          }
          $alt .= substr($spl[3], $pos+$x, 1);
        } else {
          if ($base) {
            $flag = 1 if ($base ne substr($spl[4], $sub+$x, 1));
          } else {
            $base = substr($spl[4], $sub+$x, 1);
          }
          $alt .= substr($spl[4], $sub+$x, 1) if ($2 ne 'I');
        }
      }
      $alt .= '*';
    } else {
      if ($pos) {
        $prev = substr($spl[3], $pos, $1);
      }
      if ($base && length $' != 0) {
        for (my $x = 0; $x < $1; $x++) {
          $flag = 1 if ($base ne substr($spl[4], $sub+$x, 1));
        }
      }
    }
    last if ($flag);
    $pos += $1 if ($2 ne 'I');
    $sub += $1 if ($2 ne 'D');
  }
  if (scalar @loc == 0) {
    die "No variant position in $line\n";
  }
  if ($flag) {
    print LOG "$spl[0]\t$spl[1]\t$spl[3]\t$spl[4]\t$ab\t$cig\tnot checked\n";
    $spl[9] .= ":NA";
    print OUT join("\t", @spl), "\n";
    next;
  }

  # if not already loaded, get chromosome from genome
  if ($spl[0] ne $chr) {
    $chr = $spl[0];
    $seq = "";
    for (my $x = 0; $x < 2; $x++) {
      local $/ = '>';
      while (my $chunk = <GEN>) {
        chomp $chunk;
        my @div = split("\n", $chunk);
        my $ch = shift @div;
        if ($spl[0] eq $ch) {
          $seq = join("", @div);
          last;
        }
      }
      if (! $seq) {
        # if not found, try again
        close GEN;
        open(GEN, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
        my $waste = <GEN>;
      } else {
        last;
      }
    }
  }
  if (! $seq) {
    die "Cannot find chromosome $chr in $ARGV[1]\n";
  }

  # judge chrom segments (on either side of variant)
  #my $seg1 = substr($seq, $loc[0]-$hom-2, $hom+1);
  #my $seg2 = substr($seq, $loc[$#loc], $hom+1);
  my $seg1 = substr($seq, $loc[0]-11, 10);
  my $seg2 = substr($seq, $loc[$#loc], 10);
  $seg1 =~ tr/a-z/A-Z/;
  $seg2 =~ tr/a-z/A-Z/;

  my $hit = 0;
  $seg1 =~ m/($base*)$/;
  $hit += length $1;
  $seg2 =~ m/^($base*)/;
  $hit += length $1;
  $spl[9] .= ":$hit";
  print OUT join("\t", @spl), "\n";

  print LOG "$spl[0]\t$spl[1]\t$spl[3]\t$spl[4]\t$ab\t$cig\t$seg1 $alt $seg2\t$hit";
#  if ($hit > $hom) {
#    my $dec = ($cig =~ m/[ID]/g ? ($ab > $abID ? 0 : 1) : ($ab > $abSub ? 0 : 1));
#    if ($dec) {
#      print LOG "\tremoved";
#      $rem++;
#    } else {
#      print OUT join("\t", @spl), "\n";
#      $kept++;
#    }
#  } else {
#    print OUT join("\t", @spl), "\n";
#    $kept++;
#  }
  print LOG "\n";

}
close GEN;
close VCF;
close OUT;
close LOG;

#print "Variants kept: $kept\nRemoved: $rem\n";
