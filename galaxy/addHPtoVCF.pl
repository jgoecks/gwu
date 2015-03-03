#!/usr/bin/perl

# John M. Gaspar (jmgaspar@gwu.edu)
# Jan. 2015

# Add "HP" tag to VCF that lists homopolymer length
#   adjacent to each variant. Also changes 'AB=0' to
#   'AB=x' where x is AO/DP.

use strict;
use warnings;

sub usage {
  print q(Usage: perl addHPtoVCF.pl  <infile>  <genome>  <outfile>  <logfile>
  Required:
    <infile>   Input VCF file listing variants, freebayes-style
    <genome>   Fasta file of reference genome
    <outfile>  Output VCF file
    <logfile>  Verbose log file
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

# start analyzing VCF file
my $chr = "";
my $seq = "";
my $rem = 0; my $kept = 0;
my $pr = 0;
while (my $line = <VCF>) {
  if (substr($line, 0, 1) eq '#') {
    if (substr($line, 0, 6) eq '##INFO') {
      $pr = 1;
    } elsif ($pr) {
      print OUT "##INFO=<ID=HP,Number=1,Type=Integer,",
        "Description=\"Length of adjacent homopolymer that matches variant\">\n";
      $pr = 0;
    }
    print OUT $line;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);

  # determine abundance
  if ($spl[7] !~ m/AB\=(.*?)\;/) {
    die "No AB found in $line\n";
  }
  my $ab = $1;
  if ($ab == 0) {
    $spl[7] =~ m/AO\=(.*?)\;/;
    my $ct = $1;
    $spl[7] =~ m/DP\=(.*?)\;/;
    $ab = $ct / $1;
    $ab = int(10000000*$ab+0.5)/10000000;
  }

  # determine position of variant(s) using CIGAR
  if ($spl[7] !~ m/CIGAR\=(.*?)\;/) {
    die "No CIGAR found in $line\n";
  }
  my $cig = $1;
  my @loc = ();
  my $pos = 0; my $sub = 0;
  my $alt = '*';
  my $prev = "";
  my $base = "";  # different nt
  my $flag = 0;   # 1 if complex variant
  my $del = "";   # deleted bases if variant is a deletion
  my $hit = 0;    # length of homopolymer
  while ($cig =~ m/(\d+)([IDMX])/g) {
    if ($2 ne 'M') {
      if ($prev) {
        $alt .= " $prev *";
        $prev = "";
      }
      for (my $x = 0; $x < $1; $x++) {
        push @loc, $spl[1] + $pos + $x;
        if ($2 eq 'D') {
          if ($base) {
            $flag = 1 if ($base ne substr($spl[3], $pos+$x, 1));
          } else {
            $del .= substr($spl[3], $pos+$x, 1);
          }
          $alt .= substr($spl[3], $pos+$x, 1);
        } else {
          if ($base) {
            $flag = 1 if ($base ne substr($spl[4], $sub+$x, 1));
          } else {
            for (my $y = 0; $y < length $del; $y++) {
              if (substr($spl[4], $sub+$x, 1) ne substr($del, $y, 1)) {
                $flag = 1;  # complex variant, base does not match del
                last;
              }
            }
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
        # matching bases in complex variant must match variant bases
        for (my $x = 0; $x < $1; $x++) {
          $flag = 1 if ($base ne substr($spl[4], $sub+$x, 1));
          $hit++;  # add to homopolymer length
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
    # complex variant
    print LOG "$spl[0]\t$spl[1]\t$spl[3]\t$spl[4]",
      "\t$ab\t$cig\tnot checked\n";
    $hit = 0;
  } else {

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

    # judge chrom segments (10bp on either side of variant)
    my $seg1 = substr($seq, $loc[0]-11, 10);
    my $seg2 = substr($seq, $loc[$#loc], 10);
    $seg1 =~ tr/a-z/A-Z/;
    $seg2 =~ tr/a-z/A-Z/;

    # add homopolymer lengths
    if ($del) {
      # for deletion, do not need to match $base
      $seg1 =~ m/((.)\2*)$/;
      my $hit1 = length $1;
      $seg2 =~ m/^($2*)/;
      $hit1 += length $1;

      # check other side
      $seg2 =~ m/^((.)\2*)/;
      my $hit2 = length $1;
      $seg1 =~ m/($2*)$/;
      $hit2 += length $1;

      $hit += ($hit2 > $hit1 ? $hit2 : $hit1);  # save maximum

    } else {
      $seg1 =~ m/($base*)$/;
      $hit += length $1;
      $seg2 =~ m/^($base*)/;
      $hit += length $1;
    }
    print LOG "$spl[0]\t$spl[1]\t$spl[3]\t$spl[4]\t$ab",
      "\t$cig\t$seg1 $alt $seg2\t$hit\n";
  }

  $spl[7] =~ s/AB\=(.*?)\;/AB\=$ab\;/;
  $spl[7] .= ";HP=$hit";
  print OUT join("\t", @spl), "\n";

}
close GEN;
close VCF;
close OUT;
close LOG;
