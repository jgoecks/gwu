#!/usr/bin/perl -w

#use HTML::Table;
#use CGI;
#use CGI::Pretty;
#my $cgi = new CGI;

open FILE1, "<$ARGV[0]";
open FILE2, "<$ARGV[1]";
open FILE3, "<$ARGV[2]";
open FILE4, "<$ARGV[3]";
open FILE5, "<$ARGV[4]";

my$samplename="Sample";
my@variants_to_report= <FILE1> ;
my@gene_expression_to_report= <FILE2> ;
my@sample_annotated_variants= <FILE3> ;
my@cufflinks_FPKM= <FILE4> ;
my@fusions= <FILE5> ;

close FILE1;
close FILE2;
close FILE3;
close FILE4;
close FILE5;

#########################
# Pull out FPKM Values  #
#########################

my@all_FPKM_store=();
foreach(@cufflinks_FPKM){

	# Skip header
	next if ($_ =~ m/FPKM/);
	chomp;
	my @gene_FPKM = split("\t",$_);
	push(@all_FPKM_store,$gene_FPKM[9]);

}

@all_FPKM_store=sort { $a <=> $b } @all_FPKM_store;

my@all_expression_genes=@gene_expression_to_report;

my@gene_FPKMs=();
foreach(@all_expression_genes){

	chomp$_;
	my$gene_to_look_for=$_;

	foreach(@cufflinks_FPKM){

		my@cufflinks_gene_FPKM = split("\t",$_);
		if($cufflinks_gene_FPKM[3] =~ m/^$gene_to_look_for$/){
		
			if($gene_FPKM[9]==0){
				push(@gene_FPKMs,"$cufflinks_gene_FPKM[3],$cufflinks_gene_FPKM[10]");
			}
			else{
				push(@gene_FPKMs,"$cufflinks_gene_FPKM[3],$cufflinks_gene_FPKM[9]");
			}

		}
	}
}

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

my@row=();

foreach(@sample_annotated_variants){

	chomp;
	next if ($_ =~ m/Chromosome/);
	my@annotated_variant=split("\t",$_);
	my$gene=$annotated_variant[6];
	my@aachange=split(":",$annotated_variant[8]);
	my$intron_exon_other=$annotated_variant[5];
	my$exonic_function=$annotated_variant[7];
	my$EVS=$annotated_variant[9];
	my$onethousandg=$annotated_variant[10];
	my$dbsnp=$annotated_variant[11];
	my$cosmic=$annotated_variant[13];
	my$ref_allele_count=$annotated_variant[22];
	my$variant_allele_count=$annotated_variant[23];
	
	next if ($EVS ne '' && $EVS > 0.01);
	next if ($onethousandg ne '' && $onethousandg > 0.01);
	next if ($intron_exon_other !~ m/^exonic|spicing/);
	next if ($dbsnp ne '' && $cosmic eq '');

	if($dbsnp eq''){
		$dbsnp="NA";
	}

	if($cosmic eq''){
		$cosmic="NA";
	}

	my ( $cosmicID , $therest )= split(";",$cosmic);

	foreach(@variants_to_report){
	
		chomp$_;	
		if ($_ =~ m/^$gene/){
			
			my ($nucleotide_change , $therest )=split(",",$aachange[4]);
			my$row="		\[\'$gene\', \'$aachange[1]\',\'$aachange[2]\',\'$aachange[3]\',\'$nucleotide_change\',\'$dbsnp\',\'$cosmicID\', $ref_allele_count , $variant_allele_count ],\n";
			push(@rows,$row);
			}
		}
	}

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
foreach(@fusions){

	chomp;
	my @fusion_genes = split("\t",$_);
	my $geneA=$fusion_genes[1];
	my $geneB=$fusion_genes[4];
	my $chrA=$fusion_genes[2];
	my $chrB=$fusion_genes[5];
	my $chrs=$fusion_genes="$chrA-$chrB";
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
	<p>Mutations detected in ROI genes. Muations reported here have at least 8-fold coverage and at least 2 non-reference reads. Variants found at a frequency of > 1% in the 1000 genomes project or > 1% in the Exome Variant Server, or those that appear in dbSNP are not listed.
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

