#!/usr/bin/perl

# John M. Gaspar
# July 2014

# Get reads in common from two fastq files. Make sure the
#   primers removed were from the same amplicon.

use strict;
use warnings;

sub usage {
  print q(Usage: perl getCommon.pl  <in1>  <in2>  <out1>  <out2>  <out3>  <out4>
  Required:
    <in1>   Input FASTQ file #1, with primers removed and amplicon identification
              in header (produced by removePrimer)
    <in2>   Input FASTQ file #2, with primers removed and amplicon identification
              in header (produced by removePrimer)
    <out1>  Output FASTQ file with common reads from #1
    <out2>  Output FASTQ file with common reads from #2
  Optional:
    <out3>  Output FASTQ file with singleton reads (from either #1 or #2)
    <out4>  Output FASTQ file with common reads that had different primers removed
              (like they were derived from a PCR chimera)
);
  exit;
}

usage() if (scalar @ARGV < 4 || $ARGV[0] eq "-h");

# open files
open(FQ1, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(FQ2, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT1, ">$ARGV[2]");
open(OUT2, ">$ARGV[3]");
open(SING, ">$ARGV[4]") if (scalar @ARGV > 4);
open(CHI, ">$ARGV[5]") if (scalar @ARGV > 5);

# load info from first file
my %pr;  # which amplicon the read came from
my %seq; # header, sequence, and qual scores
my $count = 0;
while (my $q = <FQ1>) {
  next if (substr($q, 0, 1) ne "@");
  my $w = <FQ1>;
  my $e = <FQ1>;
  my $r = <FQ1>;
  my @spl = split(" ", $q);

  # determine amplicon ID
  my $id = "";
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      last;
    }
  }
  die "$ARGV[0] is improperly formatted:\n  ",
    "no amplicon ID in $q\n" if (!$id);

  $pr{$spl[0]} = $id;
  $seq{$spl[0]} = $q . $w . $e . $r;
  $count++;
}
close FQ1;
print "Reads in $ARGV[0]: $count\n";

# parse second file
my $match = 0;
my $nonmatch = 0;
my $not1 = 0;
$count = 0;
while (my $q = <FQ2>) {
  next if (substr($q, 0, 1) ne "@");
  my $w = <FQ2>;
  my $e = <FQ2>;
  my $r = <FQ2>;

  # determine amplicon ID
  my @spl = split(" ", $q);
  my $id = "";
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      last;
    }
  }
  die "$ARGV[1] is improperly formatted:\n  ",
    "no amplicon ID in $q\n" if (!$id);

  if (exists $pr{$spl[0]}) {
    if ($pr{$spl[0]} eq $id) {
      # match -- print reads to outputs
      print OUT1 $seq{$spl[0]};
      print OUT2 "$q$w$e$r";
      $match++;
    } else {
      # chimera!?
      $nonmatch++;
      print CHI "$seq{$spl[0]}$q$w$e$r" if (scalar @ARGV > 5);
    }
    delete $pr{$spl[0]};
  } else {
    # singleton
    print SING "$q$w$e$r" if (scalar @ARGV > 4);
    $not1++;
  }
  $count++;
}
close FQ2;
close OUT1;
close OUT2;
close CHI if (scalar @ARGV > 5);

# singletons in $ARGV[0]
if (scalar @ARGV > 4) {
  foreach my $re (sort keys %pr) {
    print SING $seq{$re};
  }
  close SING;
}

print "Reads in $ARGV[1]: $count\n",
  "Matches: $match\n",
  "Singletons in $ARGV[0], not in $ARGV[1]: ", scalar keys %pr, "\n",
  "Singletons in $ARGV[1], not in $ARGV[0]: $not1\n",
  "Chimeras: $nonmatch\n";
