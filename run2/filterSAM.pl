#!/usr/bin/perl

# John M. Gaspar
# Sept. 2014

# Filtering a SAM file for multi-mapping reads,
#   selecting reads that map to target regions.
# (formerly called samparsemult3.pl)

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterSAM.pl  <infile1>  <infile2>  <infile3>  <infile4>  <outfile>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer)
    <infile2>  SAM file containing mapping information
    <infile3>  BED file listing locations of amplicons
    <infile4>  File listing for each amplicon whether primers are exact/close (1)
                 or way off (0).  For example:
                   PTEN_t9_5       1
                   PTEN_t9_6       1
                   BRAF_t11_1      0
                   BRAF_t15_5      0
    <outfile>  Output file (another SAM)
);
  exit;
}

usage() if (scalar @ARGV < 5 || $ARGV[0] eq "-h");

open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(SAM, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(BED, $ARGV[2]) || die "Cannot open $ARGV[2]\n"; # regions+pr.bed
open(CL, $ARGV[3]) || die "Cannot open $ARGV[3]\n";  # closePrimers.txt

# load location information for amplicons
my %loc;
my @am;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  die "$ARGV[2] is not properly formatted\n" if (scalar @spl < 4);
  $loc{$spl[3]} = "$spl[0]\t" . ($spl[1]-50) . "\t$spl[2]"; # 50bp of wiggle room at 5' end
  push @am, $spl[3];
}
close BED;

# load primers that are exact/close or way off
my %pr;  # value is 1 if exact/close, 0 if way off
while (my $line = <CL>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split("\t", $line);
  $pr{$spl[0]} = $spl[1];
}
close CL;

# load amplicon info for reads from fastq
my %amp;
my $dup = 0;
my $count = 0;
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne '@');
  chomp $line;
  $count++;
  my @spl = split(" ", $line);
  die "$ARGV[0] is not properly formatted\n" if (scalar @spl < 3);
  my $x;
  for ($x = 0; $x < scalar @am; $x++) {
    if ($spl[2] eq $am[$x]) {
      my $que = substr($spl[0], 1);
      if (exists $amp{$que}) {
        die "Primers don't match\n" if ($spl[2] ne $amp{$que});
        $dup++;
      } else {
        $amp{substr($spl[0], 1)} = $am[$x];
      }
      last;
    }
  }
  die "Cannot find amplicon $spl[2] for read $spl[0]\n" if ($x == scalar @am);
  for ($x = 0; $x < 3; $x++) {
    $line = <FQ>;
  }
}
close FQ;
print "Reads: $count\n";
print "Should be the same: ", $dup + scalar keys %amp,
  "\n";


# parse SAM, produce filtered file
open(OUT, ">$ARGV[4]");
$count = 0; my $mapped = 0;
my $ct = 0; my $un = 0;
my $line = <SAM>;
while ($line) {
  if (substr($line, 0, 1) eq '@') {
    print OUT $line;
    $line = <SAM>;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);
  $count++;
  die "Supp. alignment: $line\n" if ($spl[1] & 0x800);
  if ($spl[1] & 0x4) {
    $un++;
    print OUT "$line\n";  # print unmapped as well
    $line = <SAM>;
    next;
  }

  # load amplicon info for read
  if (! exists $amp{$spl[0]}) {
    die "Cannot find correct amplicon for $spl[0]\n";
  }
  my $ramp = $amp{$spl[0]};
  my @cut = split("\t", $loc{$ramp});  # expected location: chr $cut[0], bp $cut[1]-$cut[2]

  my $as = getAS(@spl); # get alignment score
  my @res = ();  # array for saving primary maps and secondary correct map only (at this point)
  my $idx = -1;  # index for correct map

  # check mapping(s)
  my $map = 0;  # primary map to correct location
  my $sec = 0;  # secondary map to correct location
  my $mpri = 0; # multiple primary maps?
  my $mult = 0; # multiple maps?
  push @res, $line;  # save primary map


  if (($spl[2] eq $cut[0]) && ($spl[3] > $cut[1]) && ($spl[3] < $cut[2])) {
    $map = 1;
    $idx = 0;
  }
  my $diff = 100;  # minimum difference between secondary mapping and primary
  while ($line = <SAM>) {
    chomp $line;
    my @div = split("\t", $line);
    last if ($div[0] ne $spl[0]);

    my $xs = getAS(@div);
    $mult = 1;
    if ($xs == $as) {
      push @res, $line;  # save all primary maps
      $mpri = 1;
    }
    if ($as - $xs < $diff) {
      $diff = $as - $xs;
    }
    if (($div[2] eq $cut[0]) && ($div[3] > $cut[1]) && ($div[3] < $cut[2])) {
      if ($xs == $as) {
        $map = 1;
      } elsif (!$map) {
        push @res, $line;
        $sec = 1;
      }
      $idx = $#res;
    }
  }

  # check primers for multi-mapping possibility
  my $clo = 0;
  if (($map && $mult) || $sec) {
    if (exists $pr{$ramp}) {
      $clo = $pr{$ramp};
    } elsif ($diff < 11) {
      # print warning if secondary map is close (within 10)
      print "No multi-mapping considered for $ramp and diff is $diff (read $spl[0])\n";
    }
  }

  # consider only reads that mapped to right location
  if ($map) {
    die "Unknown correct location for $spl[0]\n" if ($idx == -1);

    if ($mult) {
      if ($mpri) {
        if ($clo) {
          # multiple primary, and primers close
          # if correct read not listed 1st, move to 1st
          if ($idx) {
            my $temp = $res[0];
            my @cut = split("\t", $res[$idx]);
            $cut[1] = ($cut[1] & 0xEFF);
            $res[0] = join("\t", @cut);
            @cut = split("\t", $temp);
            $cut[1] = ($cut[1] | 0x100);
            $res[$idx] = join("\t", @cut);
            $idx = 0;
          }

          # print all primary results
          for (my $x = 0; $x < scalar @res; $x++) {
            print OUT "$res[$x]\n";
          }
        } else {
          # multiple primary, but primers not close
          my @cut = split("\t", $res[$idx]);
          $cut[1] = ($cut[1] & 0xEFF);
          print OUT join("\t", @cut), "\n";
        }
        #$mmap{$ramp}++;
      } else {
        print OUT "$res[$idx]\n"; # 1 primary read is correct: print it
        #$primap{$ramp}++;
      }
    } else {
      print OUT "$res[$idx]\n"; # 1 read is correct: print it
      #$sinmap{$ramp}++;
    }
    $mapped++;
  } elsif ($sec) {
    # if primers aren't close, promote the read
    if (!$clo) {
      my @cut = split("\t", $res[$idx]);
      $cut[1] = ($cut[1] & 0xEFF);
      print OUT join("\t", @cut), "\n";
      $mapped++;
    } else {
    # else read probably came from other location: list as unmapped
      my @cut = split("\t", $res[0]);
      $cut[1] = 0x4;
      print OUT join("\t", @cut), "\n";
      $ct++;
    }
    #$smap{$ramp}++;  # secondary mapping
  } else {
    my @cut = split("\t", $res[0]);
    $cut[1] = 0x4;
    print OUT join("\t", @cut), "\n";
    #$umap{$ramp}++;  # unmapped to correct location
    $ct++;
  }

}
close SAM;
close OUT;
print "Reads analyzed: $count\n",
  "Unmapped: $un\n",
  "Kept as mapped: $mapped\n",
  "Changed to unmapped: $ct\n";


# get alignment score from sam line
sub getAS {
  my @spl = @_;
  my $as = 1;
  for (my $x = 11; $x < scalar @spl; $x++) {
    my @div = split(":", $spl[$x]);
    if ($div[0] eq "AS") {
      $as = pop @div;
      last;
    }
  }
  die "Cannot find alignment score\n" if ($as == 1);
  return $as;
}
