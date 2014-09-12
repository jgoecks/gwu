# get only reads whose pair also had primer removed

use strict;
use warnings;

die "Need files on command line\n" if (scalar @ARGV < 5);
open(FQ, $ARGV[0]) || die "Cannot open\n";
open(F2, $ARGV[1]) || die "Cannot open\n";
my %pr;
my %seq;
my $count = 0;
while (my $q = <FQ>) {
  next if (substr($q, 0, 1) ne "@");
  my $w = <FQ>;
  my $e = <FQ>;
  my $r = <FQ>;
  my @spl = split(" ", $q);
  die "No primer match!?\n$q\n" if (scalar @spl < 3);
  $pr{$spl[0]} = $spl[2];
  $seq{$spl[0]} = $q . $w . $e . $r;
  $count++;
}
close FQ;
print "Reads loaded: $count\n";

my $match = 0;
my $nonmatch = 0;
my $not1 = 0;
$count = 0;
open(OUT, ">$ARGV[2]");
open(OU2, ">$ARGV[3]");
open(SIN, ">$ARGV[4]");
open(MIS, ">$ARGV[5]");
while (my $q = <F2>) {
  next if (substr($q, 0, 1) ne "@");
  my $w = <F2>;
  my $e = <F2>;
  my $r = <F2>;
  my @spl = split(" ", $q);
  die "No primer match!?\n$q\n" if (scalar @spl < 3);
  if (exists $pr{$spl[0]}) {
    if ($pr{$spl[0]} eq $spl[2]) {
      print OUT $seq{$spl[0]};
      print OU2 "$q$w$e$r";
      $match++;
    } else {
      $nonmatch++;
      print MIS "$seq{$spl[0]}$q$w$e$r";
      #print "$spl[0] $pr{$spl[0]}\t$spl[2]\n";
      #my $wait = <STDIN>;
    }
    delete $pr{$spl[0]};
  } else {
    print SIN "$q$w$e$r";
    $not1++;
    #print "$q\n";
    #my $wait = <STDIN>;
  }
  $count++;
}
close F2;
close OUT;
close OU2;

foreach my $re (keys %pr) {
  print SIN $seq{$re};
}
close SIN;
close MIS;

print "Query reads: $count\n",
  "Matches: $match\n",
  "Nonmatches: $nonmatch\n",
  "in 2, not in 1: $not1\n",
  "in 1, not in 2: ", scalar keys %pr, "\n";
