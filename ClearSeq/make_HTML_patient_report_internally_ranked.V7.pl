#!/usr/bin/perl

# Produces an HTML report for RNAseq.
# (edited JMG, 2/2016:
#   - using ANNOVAR output, not reformatted
#   - more robust parsing of Cufflinks output
#   - not loading all input files to memory
# )

use strict;
use warnings;

die "Error! Need five files on command-line\n" if (scalar @ARGV < 5);

open FILE1, "<$ARGV[0]" || die "Cannot open $ARGV[0]\n"; # list of genes in which to report variants
open FILE2, "<$ARGV[1]" || die "Cannot open $ARGV[1]\n"; # list of genes of which to report expression
open FILE3, "<$ARGV[2]" || die "Cannot open $ARGV[2]\n"; # variants annotated by ANNOVAR
open FILE4, "<$ARGV[3]" || die "Cannot open $ARGV[3]\n"; # Cufflinks FPKM
open FILE5, "<$ARGV[4]" || die "Cannot open $ARGV[4]\n"; # Text output from tophat_fusion_post

my$samplename="Sample";

# load lists of genes
my @variants_to_report = <FILE1>;
close FILE1;
chomp @variants_to_report;

my @all_expression_genes = <FILE2>;
close FILE2;
chomp @all_expression_genes;


#########################
# Pull out FPKM Values  #
#########################

# load header info from Cufflinks
my @info = qw(gene_id FPKM);
my @idx;
for (my $x = 0; $x < scalar @info; $x++) { $idx[$x] = -1; }
my $maxIdx = -1;
my $line = <FILE4>;
chomp $line;
my @head = split("\t", $line);
for (my $x = 0; $x < scalar @head; $x++) {
  for (my $y = 0; $y < scalar @info; $y++) {
    if ($head[$x] eq $info[$y]) {
      die "Error! Repeat of header item $info[$y]\n" if ($idx[$y] != -1);
      $idx[$y] = $x;
      $maxIdx = $x if ($x > $maxIdx);
      last;
    }
  }
}
for (my $x = 0; $x < scalar @idx; $x++) {
  die "Error! No Cufflinks info for $info[$x]\n" if ($idx[$x] == -1);
}

# load FPKM values
my @all_FPKM_store;
my @gene_FPKMs;
while (my $line = <FILE4>) {
  chomp $line;
  my @spl = split("\t", $line);
  die "Error! Not enough info in Cufflinks record:\n$line\n" if ($#spl < $maxIdx);

  push(@all_FPKM_store, $spl[$idx[1]]); # save all FPKM values

  # load FPKM values for genes of interest
  for (my $x = 0; $x < scalar @all_expression_genes; $x++) {
    if ($spl[3] eq $all_expression_genes[$x]) {
      if ($spl[9] == 0) {
        push(@gene_FPKMs, "$spl[3],$spl[10]");
      } else {
        push(@gene_FPKMs, "$spl[3],$spl[9]");
      }
      last;
    }
  }

}
close FILE4;
@all_FPKM_store = sort { $a <=> $b } @all_FPKM_store;


my@rows5=();
foreach(@gene_FPKMs){

	my ( $gene , $fpkm ) = split(",",$_);
	my$position_counter=1;
	my$logfpkm=log10($fpkm + 1);


	foreach(@all_FPKM_store){

		if($_ < $fpkm){
		
		$position_counter++;
		}
	}
	
	my$percentile=int(100*($position_counter/scalar(@all_FPKM_store)));
	my$row="\t\t\[ \'$gene\' , $logfpkm , $percentile \],\n";
	push(@rows5,$row);

}


######################################################
# Pull out coding mutations from annotated variants  #
######################################################

my@rows=();

# load header info from ANNOVAR
@info = qw(Func.refGene Gene.refGene AAChange.refGene 1000g esp snp cosmic);
@idx = ();
for (my $x = 0; $x < scalar @info; $x++) { $idx[$x] = -1; }
$maxIdx = -1;
$line = <FILE3>;
chomp $line;
@head = split("\t", $line);
for (my $x = 0; $x < scalar @head; $x++) {
  for (my $y = 0; $y < scalar @info; $y++) {
    if ($head[$x] =~ m/$info[$y]/) {
      next if ($y < 3 && $head[$x] ne $info[$y]);
      die "Error! Repeat of header item $info[$y]\n" if ($idx[$y] != -1);
      $idx[$y] = $x;
      $maxIdx = $x if ($x > $maxIdx);
      last;
    }
  }
}
for (my $x = 0; $x < scalar @idx; $x++) {
  die "Error! No ANNOVAR info for $info[$x]\n" if ($idx[$x] == -1);
}

while (my $line = <FILE3>) {
  chomp $line;
  my @var = split("\t", $line);
  die "Error! Not enough info for variant:\n$line\n" if ($#var < $maxIdx);

  # skip variants
  next if ($var[$idx[0]] !~ m/(exonic|splicing)/);
  next if ($var[$idx[3]] ne '.' && $var[$idx[3]] > 0.01); # 1000genomes
  next if ($var[$idx[4]] ne '.' && $var[$idx[4]] > 0.01); # esp
  next if ($var[$idx[5]] ne '.' && $var[$idx[6]] eq '.'); # in dbsnp, not cosmic

  $var[$idx[5]] = "NA" if ($var[$idx[5]] eq '.');

  # save 1st cosmic match
  if ($var[$idx[6]] eq '.') {
    $var[$idx[6]] = "NA";
  } else {
    my @arr = split(';', $var[$idx[6]]);
    $var[$idx[6]] = $arr[0];
  }

  my @aachange = split(',', $var[$idx[2]]);
  my @aainfo = split(':', $aachange[0]);

  foreach (@variants_to_report) {
    chomp $_;
    if ($_ =~ m/^$var[$idx[1]]/) {
      my $row="\t\t\[\'$var[$idx[1]]\', \'$aainfo[1]\',\'$aainfo[2]\',\'$aainfo[3]\',\'$aainfo[4]\',\'$var[$idx[5]]\',\'$var[$idx[6]]\' ],\n";
      push(@rows,$row);
    }
  }
}
close FILE3;

#################################
# Heatmap of missense and indel #
#################################

my@tmp_rows=@rows;
my@rows6=();
push(@rows6,"		data.addColumn('number', 'test');\n");
foreach(@variants_to_report){
	
	chomp$_;
	my $gene_to_report=$_;	
	push(@rows6,"		data.addColumn('number', '$gene_to_report')\n");

	}

push(@rows6,"		data.addRows(2);\n");
push(@rows6,"		data.setCell(1, 0, 'Nonsense');\n");
push(@rows6,"		data.setCell(0, 0, 'Missense');\n");
push(@rows6,"		data.setCell(0, 1, 0);\n");
push(@rows6,"		data.setCell(1, 1, 5);\n");

my$genecounter=2;
foreach(@variants_to_report){
	
	chomp$_;
	my $gene_to_report=$_;	
	my $mutation_counter=0;	
	my $indel_counter=0;

	foreach(@tmp_rows){

		s/\s//g;
		s/'//g;
		s/[[]//g;

		my ( $gene , $aachange1 , $aachange2 , $aachange3, $aachange4 , $rest ) = split(",",$_);


		if(($gene_to_report =~ m/^$gene$/)&&($aachange4 !~ m/fs|X/)){
			$mutation_counter++;
		}
		
		if(($gene_to_report =~ m/^$gene$/)&&($aachange4 =~ m/fs|X/)){
			$indel_counter++;			
			}
		}

	if($mutation_counter > 0){
		
		push(@rows6,"		data.setCell(0, $genecounter, 0);\n");
		}

	else{
		push(@rows6,"		data.setCell(0, $genecounter, null);\n");
		}

	if($indel_counter > 0){
		
		push(@rows6,"		data.setCell(1, $genecounter, 5);\n");
		}

	else{
		push(@rows6,"		data.setCell(1, $genecounter, null);\n");
		}

	$genecounter++;

	}

#################
# Gene Fusions  #
#################

my@rows4=();
while (my $line = <FILE5>) {
  chomp $line;
  my @fusion_genes = split("\t", $line);
	my $geneA=$fusion_genes[1];
	my $geneB=$fusion_genes[4];
	my $chrA=$fusion_genes[2];
	my $chrB=$fusion_genes[5];
	my $chrs="$chrA-$chrB";
	my $LHS=$fusion_genes[3];
	my $RHS=$fusion_genes[6];
	my $strands="na";
	my $support=$fusion_genes[7] + $fusion_genes[8];

	foreach(@all_expression_genes){	
		chomp$_;
		if (($geneA =~ m/^$_$/)||($geneB =~ m/^$_$/)){
			my$row="\t\t\[\'$geneA-$geneB\',\' FUSION \',\' $chrs \',\' $LHS \',\' $RHS \',\' $strands \', $support \],\n";
			push(@rows4,$row);
		}
	}
}
close FILE5;

###############
# HTML Output #
###############

print"
<html>

  <head>
  <h1>$samplename RNASeq Report</h1>
  <h2>Coding Variants in Recurrently Mutated Genes</h2>

<script type='text/javascript' src='http://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {});
      google.load('prototype', '1.6');
    </script> 
    <script type='text/javascript' src='http://systemsbiology-visualizations.googlecode.com/svn/trunk/src/main/js/load.js'></script>
    <script type='text/javascript'>
        systemsbiology.load('visualization', '1.0', {packages:['bioheatmap']});
    </script>
    <script type='text/javascript'>
      google.setOnLoadCallback(drawHeatMap);
      function drawHeatMap() {
          var data = new google.visualization.DataTable();
          data.addColumn('string', 'Gene Name');
          ";
	print @rows6;
	print"
          heatmap = new org.systemsbiology.visualization.BioHeatMap(document.getElementById('heatmapContainer'));
          heatmap.draw(data, {});
      }
    </script>
  </head>

  <body>
    <div id='heatmapContainer'></div>
  </body>
	<p>Mutations detected in ROI genes. Mutations reported here have at least 8-fold coverage and at least 2 non-reference reads. Variants found at a frequency of > 1% in the 1000 genomes project or > 1% in the Exome Variant Server, or those that appear in dbSNP, are not listed.
</p>
    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {packages:['table']});
      google.setOnLoadCallback(drawTable);
      function drawTable() {
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Gene');
 	data.addColumn('string', 'RefSeq ID');
 	data.addColumn('string', 'Exon');
 	data.addColumn('string', 'Nucleotide Change');
 	data.addColumn('string', 'Amino Acid Change');
 	data.addColumn('string', 'DBSNP ID (build 137)');
 	data.addColumn('string', 'COSMIC ID');
        data.addColumn('number', 'Number of Reference Reads');
	data.addColumn('number', 'Number of Variant Reads');

        data.addRows([";
 
	print @rows;
	print"

        ]);

        var table = new google.visualization.Table(document.getElementById('table_div'));
        table.draw(data, {showRowNumber: true});
      }
    </script>
  </head>

  <body>
    <div id='table_div' style='width: 900px;'></div>
  </body>

  <h2>Expression of Clinically Relevant and Commonly Translocated Genes</h2>

<p>Gene expression estimates result from cufflinks quanitfication of RefSeq genes. Gene FPKM values are displayed on a log scale.
FPKM = Fragments per kilobase of exon per million mapped reads.</p>

  <head>
    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {packages: ['corechart']});
    </script>
    <script type='text/javascript'>
      function drawVisualization() {
        // Some raw data (not necessarily accurate)
        var data = google.visualization.arrayToDataTable([
          ['Gene', 'Log10(FPKM+1)', 'Expression Percentile Rank'],
	";

	print "\n";
	print @rows5;
	print"

	]);

        var options = {
         title : 'Gene Expression',
        vAxes: {
		0: {title: 'Log10(FPKM+1)'}, 
		1: {title: 'Percentile Rank'}},
        hAxis: {title: 'Gene'},
        series: {
            0:{ type: 'bars', targetAxisIndex: 0 },
            1: { type: 'line', targetAxisIndex: 1}
        	}
        };

        var chart = new google.visualization.ComboChart(document.getElementById('chart_combo'));
        chart.draw(data, options);
      }
      google.setOnLoadCallback(drawVisualization);
    </script>
  </head>
  <body>
    <div id='chart_combo' style='width: 900px; height: 500px;'></div>
  </body>


  <h2>Fusion Transcript Predictions</h2>
<p>Fusion transcript predictions result from a modifed version on Tophat Fusion. Predictions listed here are limited to the above 
therapeutically relevant and known translocation genes.</p>

    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {packages:['table']});
      google.setOnLoadCallback(drawTable);
      function drawTable() {
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Potential Fusion / Disruption');
 	data.addColumn('string', 'Type');
 	data.addColumn('string', 'Chromosomes');
 	data.addColumn('string', 'Mapping1');
 	data.addColumn('string', 'Mapping2');
 	data.addColumn('string', 'Strands');
        data.addColumn('number', 'Supporting Reads');

        data.addRows([";
 
	print "\n";

	my@unique_rows4=uniq(@rows4);
	print @unique_rows4;

	print"
        ]);

        var table = new google.visualization.Table(document.getElementById('table_div2'));
        table.draw(data, {showRowNumber: true});
      }

    </script>
</head>
  <body>
    <div id='table_div2' style='width: 900px;'></div>
  </body>
</html>";

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

sub log10 {
my $n = shift;
return log($n)/log(10);
}

