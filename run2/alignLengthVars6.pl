#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Align length variants to genomic segments.

use strict;
use warnings;

sub usage {
  print q(Usage: perl alignLengthVars.pl  <infile1>  <infile2>  <infile3>  <outfile>  <logfile>
  Required:
    <infile1>  File containing length-variant reads without primers attached
                 (produced by getLengthVars.pl)
    <infile2>  File listing primer and target sequences (produced by getPrimers.pl)
    <infile3>  BED file listing locations of primers
    <outfile>  Output file containing SAM mapping info
    <logfile>  Verbose output file listing possible CIGARs and scores for each read
);
  exit;
}

usage() if (scalar @ARGV < 5 || $ARGV[0] eq "-h");

# open files
open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(PR, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(BED, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(OUT, ">$ARGV[3]");
open(LOG, ">$ARGV[4]");

# load primer and target sequences
my %fwd; my %rev; my %targ;
while (my $line = <PR>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split(',', $line);
  next if (scalar @spl < 4);
  $fwd{$spl[0]} = $spl[1];
  $rev{$spl[0]} = $spl[2];
  $targ{$spl[0]} = $spl[3];
}
close PR;

# load expected amplicon location
my %loc;
while (my $line = <BED>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print "Improperly formatted line in BED file: $line\n",
      "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $loc{$spl[3]}) {
    my @div = split("\t", $loc{$spl[3]});
    $loc{$spl[3]} = "$spl[0]\t".($spl[2]+1)
      if ($spl[2]+1 < $div[1]);
  } else {
    $loc{$spl[3]} = "$spl[0]\t".($spl[2]+1);
  }
}
close BED;

# load (unique) reads from fastq
my $count = 0; my $uniq = 0;
my %seq;  # amplicons are keys, lengths are keys to %{$seq{$amp}},
          #   seqs are keys to %{$seq{$amp}{$len}},
          # values are lists of reads that share the seq (comma-separated)
my %qual; # for quality scores
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne '@');
  chomp $line;
  my @spl = split(" ", $line);  # $spl[0] is name, $spl[2] is amplicon
                                # $spl[$#spl] is cig I/D
  die "$ARGV[0] is improperly formatted\n"
    if (scalar @spl < 5 || $spl[$#spl] !~ m/[ID]/);
  my $read = substr($spl[0], 1);
  $line = <FQ>;
  chomp $line;
  my $len = length $line;
  if (exists $seq{$spl[2]}{$len}{$line}) {
    $seq{$spl[2]}{$len}{$line} .= ",$read";
  } else {
    $seq{$spl[2]}{$len}{$line} = "$spl[$#spl] $read";
    $uniq++;
  }
  $count++;
  $line = <FQ>;
  # save quality scores
  $line = <FQ>;
  chomp $line;
  $qual{substr($spl[0], 1)} = $line;
}
close FQ;
print "Reads loaded: $count\n",
  "Unique reads: $uniq\n";

# align reads to genomic segments
print LOG "Amplicon\tReads\tScore\tCIGAR(s)\tSequence\tReads";
print OUT "#Read\tAmplicon\tSubs\tIn/del\tChrom\tPos\tCIGAR",
  "\tSeq\tQual\tMD\tExternal\tScore\n";
foreach my $am (sort keys %seq) {

  # loop through length variants
  foreach my $len (sort keys %{$seq{$am}}) {

    if (! exists $loc{$am}) {
      print "Error! No location info for $am\n";
      last;
    } elsif (! exists $targ{$am}) {
      print "Error! No sequences loaded for $am\n";
      last;
    }

    print LOG "\n$am len=$len\n";
    my %rep = ();  # for I/D and list of reads
    my %cig = ();  # saves counts of reads that match each cigar
    my $tot = 0;   # total number of reads

    # loop through unique reads
    while (my $max = findMax(%{$seq{$am}{$len}})) {
      my @que = split(" ", $max); # $que[0] is sequence, $que[1] is read count
      $tot += $que[1];
      my @spl = split(" ", $seq{$am}{$len}{$que[0]}); # $spl[0] is I/D, $spl[1] lists reads

      # get possible CIGARs
      my $aln = alignSeq($que[0], $targ{$am}, $spl[0]);
      my @cut = split(" ", $aln);
      my $score = shift @cut;
      for (my $x = 0; $x < scalar @cut; $x++) {
        $cig{$cut[$x]} += $que[1];
      }

      # log info
      print LOG "\t$que[1]\t$score\t", join(",", @cut),
        "\t$que[0]\t$spl[1]\n";

      $rep{$que[0]} = "$spl[0] $spl[1]"; # save I/D and list of reads
      delete $seq{$am}{$len}{$que[0]};
    }

    # find consensus CIGAR
    my $max = 0;
    my $firstM = 1000;  # first M of cigar -- 0 indicates external in/del
    my $best = "";      # best cigar
    foreach my $key (keys %cig) {
      if ($cig{$key} > $max) {
        $max = $cig{$key};
        $best = $key;
        if ($key =~ m/\D0M$/) {
          $firstM = 0;
        } else {
          $key =~ m/^(\d+)M/;
          $firstM = $1;
        }
      } elsif ($cig{$key} == $max) {
        # preference for external in/dels
        if ($key =~ m/\D0M$/) {
          $best = $key;
          $firstM = 0;
        } else {
          $key =~ m/^(\d+)M/;
          if ($1 < $firstM) {
            $firstM = $1;
            $best = $key;
          }
        }
      }
    }
    printf LOG "consensus\t$max\t%.1f%%\t$best\n", 100*$max/$tot;
    print "Warning for $am, len=$len: $best not perfect\n" if ($max != $tot);

    # evaluate external deletions for artifacts (i.e. if primer is internal match)
    #   (skip external insertions -- how to evaluate them?)
    my $score = -1;
    if ((!$firstM) && $best =~ m/D/) {
      my $prim; my $del;  # primer and genomic sequences
      if ($best =~ m/\D0M$/) {
        # deletion at 3' end
        $prim = reverse $rev{$am};
        $del = reverse (substr($targ{$am}.$rev{$am}, $len, length $prim));
      } else {
        # deletion at 5' end
        $prim = $fwd{$am};
        $del = substr($fwd{$am}.$targ{$am}, - $len - length $prim, length $prim);
      }

      # score match of primer to genomic segment
      $score = scoreAlign($del, $prim);
      printf LOG "\texternal\t%.3f\n", $score;
    }

    # produce output -- info for a SAM record
    while (my $max = findMax(%rep)) {
      my @que = split(" ", $max); # $que[0] is sequence, $que[1] is read count
      my @cut = split(" ", $rep{$que[0]}); # $cut[0] has in/del length
      my @spl = split(",", $cut[1]); # list of reads
      my ($md0, $md1) = getMD($que[0], $targ{$am}, $best);
      for (my $x = 0; $x < scalar @spl; $x++) {
        print OUT "$spl[$x]\t$am\t$md1\t$cut[0]\t$loc{$am}\t$best\t$que[0]",
          "\t$qual{$spl[$x]}\t$md0";
        printf OUT "\texternal\t%.3f", $score if ($score != -1);
        print OUT "\n";
      }
      delete $rep{$que[0]};
    }
  }
}
close LOG;
close OUT;

# produce MD flag
sub getMD {
  my $que = $_[0];
  my $ref = $_[1];
  my $cig = $_[2];

  my $md = "";
  my $qpos = 0; my $rpos = 0;
  my $match = 0;
  my $sub = 0;
  while ($cig =~ m/(\d+)([IMD])/g) {
    my $len = $1;
    if ($2 eq 'M') {
      for (my $x = 0; $x < $len; $x++) {
        if (substr($que, $qpos + $x, 1) ne
            substr($ref, $rpos + $x, 1)) {
          $md .= $match . substr($ref, $rpos + $x, 1);
          $match = 0;
          $sub++;
        } else {
          $match++;
        }
      }
      $rpos += $len;
      $qpos += $len;
    } elsif ($2 eq 'D') {
      if ($match) {
        $md .= $match;
      } elsif (!$md) {
        $md = "0";
      }
      $match = 0;
      $md .= '^' . substr($ref, $rpos, $len);
      $rpos += $len;
    } elsif ($2 eq 'I') {
      $qpos += $len;
    }
  }
  $md .= $match if ($match);
  return $md, $sub;
}

# align sequences -- return cigar(s)
sub alignSeq {
  my $que = $_[0];
  my $ref = $_[1];

  # get cigar info
  $_[2] =~ m/(\d+)([ID])/;
  my $gap = $1;
  my $ins = 0;
  # rearrange so $que is shorter sequence
  if ($2 eq 'I') {
    $ins = 1;
    my $temp = $que;
    $que = $ref;
    $ref = $temp;
  }

  my $maxscore = -1000;  # best alignment score (actual max. is 0)
  my $cig = "";   # best alignment cigar

  # check each position -- $x is 5' end of gap
  for (my $x = 0; $x < 1 + length $que; $x++) {
    my $score = 0;
    my $rpos = 0;
    my $flag = 0;
    for (my $y = 0; $y < length $que; $y++) {
      $rpos += $gap if ($y == $x);
      if (substr($que, $y, 1) ne substr($ref, $rpos + $y, 1)) {
        if (--$score < $maxscore) {
          $flag = 1;
          last;
        }
      }
    }

    # save cigar
    if ($flag) {
      next;
    } elsif ($score == $maxscore) {
      $cig .= " ${x}M$gap" . ($ins ? "I" : "D") . ((length $que) - $x) . "M";
    } else {
      $maxscore = $score;
      $cig = "${x}M$gap" . ($ins ? "I" : "D") . ((length $que) - $x) . "M";
    }
  }

  return "$maxscore $cig";
}

# find sequence in hash with the most reads
sub findMax {
  my %hash = @_;
  my $max = 0;
  my $seq = "";
  foreach my $key (keys %hash) {
    my @spl = split(",", $hash{$key});
    if (scalar @spl > $max) {
      $max = scalar @spl;
      $seq = $key;
    }
  }
  return $max ? "$seq $max" : "";
}

# score alignment -- no in/dels allowed
# weighting function (from 5' end):
#   bases before last 20: 1
#   bases 1-10:           2*pos
#   bases 11-19:          3*pos
#   base 20:              5*pos
sub scoreAlign {
  my $que = $_[0];
  my $ref = $_[1];
  my $len = length $que;

  my $score = 0; my $total = 0;
  for (my $x = 0; $x < $len; $x++) {
    my $val = 21 - $len + $x;  # value of match at this position
    if ($val > 19) {
      $val *= 5;
    } elsif ($val > 10) {
      $val *= 3;
    } elsif ($val > 0) {
      $val *= 2;
    } else {
      $val = 1;
    }
    $score += $val if (substr($ref, $x, 1) eq substr($que, $x, 1));
    $total += $val;
  }
  return $score / $total;
}
