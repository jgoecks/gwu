# converts a pileup back to a fastq

use strict;
use warnings;

die "Need files\n" if (scalar @ARGV < 2);

open(IN, $ARGV[0]) || die "Cannot open\n";
open(OUT, ">$ARGV[1]");
my @read;
my @st;
my @seq;
my @qual;
my @mqu;
my $y = 0;
while (my $line = <IN>) {
  chomp $line;
  #print "$line\n";
  my @spl = split("\t", $line);
  my $pos = 0;
  my $done = 0;
  for (my $x = 0; $x < $spl[3]; $x++) {
    my $ch = substr($spl[4], $pos, 1);
    my $qu = substr($spl[5], $x, 1);
    if ($ch eq '^') {
      die "$line\nAdding things out of order: pos $x is not ", scalar @read,
        "\n" if (scalar @read != $x);
      $mqu[$x] = ord(substr($spl[4], ++$pos, 1)) - 33;  # save mapping quality
      $read[$x] = $spl[0];
      $st[$x] = $spl[1];
      $ch = substr($spl[4], ++$pos, 1);
    }

    if (($ch eq '.') || ($ch eq ',')) {
      $ch = $spl[2];
    } elsif ($ch eq '*') {
      $ch = "";
      $qu = "";
    }

    my $next = substr($spl[4], $pos+1, 1);
    if ($next eq '-') {
      my $num = substr($spl[4], $pos+2, 1);
      $pos += $num + 2;
    } elsif ($next eq '+') {
      my $num = substr($spl[4], $pos+2, 1);
      $ch .= substr($spl[4], $pos+3, $num);
      $pos += $num + 2;
      #$qu .= substr($  # no quality scores for inserted bases!? check.pl
      #print "$line\n$ch\n";
      #my $wait = <STDIN>;
    }
    $ch =~ tr/a-z/A-Z/;
    #print "Read $x: char is $ch, qual is $qu\n";
    #my $wait = <STDIN>;
    $seq[$x] .= $ch;
    $qual[$x] .= $qu;

    if (substr($spl[4], ++$pos, 1) eq '$') {
      $pos++;
      print OUT '@Read', "$y $read[$x] $st[$x] $mqu[$x]\n$seq[$x]\n+\n$qual[$x]\n";
      $y++;
      $st[$x] = -1;
      $done++;
    }
  }

  # remove finished sequences
  if ($done) {
    my $skip = 0;
    for (my $z = 0; $z < scalar @st - $done; $z++) {
      #while ($st[$z] == -1) {
      #  $read[$z] = $read[$z + $skip];
      #  $st[$z] = $st[$z + $skip];
      #  $mqu[$z] = $mqu[$z + $skip];
      #  $seq[$z] = $seq[$z + $skip];
      #  $qual[$z] = $qual[$z + $skip];
        #$z--;
      #}
      my $a;
      for ($a = 0; $a < scalar @st - $skip - $z; $a++) {
        last if ($st[$z + $skip + $a] != -1);
      }
      $skip += $a;
      $read[$z] = $read[$z + $skip];
      $st[$z] = $st[$z + $skip];
      $mqu[$z] = $mqu[$z + $skip];
      $seq[$z] = $seq[$z + $skip];
      $qual[$z] = $qual[$z + $skip];
    }
    for (my $z = 0; $z < $done; $z++) {
      pop @read;
      pop @st;
      pop @mqu;
      pop @seq;
      pop @qual;
    }
  #print "Finished: $done\nStill processing: ", scalar @read,
  #  "\n$spl[0]\t$spl[1]\t$spl[3]\n";
#  print "Reads parsed: $pos\n";
#  for (my $z = 0; $z < $pos; $z++) {
#    print "read$z  $seq[$pos]  $qual[$pos]\n";
#  }
  #my $wait = <STDIN>;

  }

}
close IN;
close OUT;

print "Reads printed: $y\n";
