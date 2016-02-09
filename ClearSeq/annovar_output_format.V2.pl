#!/usr/bin/perl
# Formats the output from the table_annovar.pl script within Annovar.
# Gets rid of intronic an intergenic varaints.
# Gets rid of some of the functional outputs as I don't know what they are.

open (ANNOVAR, "<$ARGV[0]") or die "cannot open $ARGV[0]: $!";

my @ANNOVAR=<ANNOVAR>;

print "Chromosome	Start Position	End Position	Reference Allele	Variant Allele	Intron/Exon/Other	Gene:RefGene	ExonicFunction:RefGene	Amino Acid Change:RefGene	 Occurrences in Exome Variant Server 6500 samples (http://evs.gs.washington.edu/EVS/)	 Occurrences in 1000 genomes 2012apr all ethnicities	DBSNP 137 ID	Occurrences in COSMIC build 67	COSMIC Details	SIFT	PolyPhen2_HDIV	 Average per-sample depth of bases with Phred score >= 15	Somatic Status	Genotype	Genotype Quality	Raw Read Depth as reported by SAMtools	Quality Read Depth of bases with Phred score >= 15	Depth of reference-supporting bases (reads1)	Depth of variant-supporting bases (reads2)	Variant allele frequency	P-value from Fisher's Exact Test (germline or unknown variants)	P-value from Fisher's Exact Test (Somatic or LOH variants)	Average quality of reference-supporting bases (qual1)	Average quality of variant-supporting bases (qual2)	Depth of reference-supporting bases on forward strand (reads1plus)	Depth of reference-supporting bases on reverse strand (reads1minus)	Depth of variant-supporting bases on forward strand (reads2plus)	Depth of variant-supporting bases on reverse strand (reads2minus)\n";

foreach(@ANNOVAR){

	chomp$_;
	my@annovar=split(/\t/,$_);			
	my$intron_exon=$annovar[5];

	my$tmp=$annovar[-3];
#	my$tmp=$annovar[-4];
	$tmp=~s/SOMATIC;//g;
	my@info=split(/:/,$annovar[-1]);
	my@zygosity=split(/;/,$tmp);

	my$Chr=$annovar[0];	
	my$Start=$annovar[1];	
	my$End=$annovar[2];	
	my$Ref=$annovar[3];	
	my$Alt=$annovar[4];	
	my$Func_refGene=$annovar[5];	
	my$Gene_refGene=$annovar[6];	
	my$ExonicFunc_refGene=$annovar[7];	
	my$AAChange_refGene=$annovar[8];	
	my$esp6500si_all=$annovar[9];	
	my$o1000g2012apr_all=$annovar[10];	
	my$snp137=$annovar[11];	
	my$LJB2_SIFT=$annovar[12];	
	my$LJB2_PolyPhen2_HDIV=$annovar[13];	
	my$LJB2_PP2_HDIV_Pred=$annovar[14];	
	my$LJB2_PolyPhen2_HVAR=$annovar[15];	
	my$LJB2_PolyPhen2_HVAR_Pred=$annovar[16];	
	my$LJB2_LRT=$annovar[17];	
	my$LJB2_LRT_Pred=$annovar[18];	
	my$LJB2_MutationTaster=$annovar[19];	
	my$LJB2_MutationTaster_Pred=$annovar[20];	
	my$LJB_MutationAssessor=$annovar[21];	
	my$LJB_MutationAssessor_Pred=$annovar[22];	
	my$LJB2_FATHMM=$annovar[23];	
	my$LJB2_GERP=$annovar[24];	
	my$LJB2_PhyloP=$annovar[25];	
	my$LJB2_SiPhy=$annovar[26];	
	my$cosmic67=$annovar[27];	
	my$GT;
	my$GQ;
	my$SDP;
	my$RD;
	my$AD;
	my$FREQ;
	my$PVAL;
	my$PVAL2;
	my$RBQ;
	my$ABQ;
	my$RDF;
	my$RDR;
	my$ADF;
	my$ADR;
	my$ADP;
	my$WT;

	#Number of COSMIC samples
	my( $IDequals , $IDs , $numbers)=split("=",$annovar[27]);
	# Get rid of anything hat isn't a number or a comma	
	$numbers =~ tr/0-9,//cd;
	my@cosmic_counts=split(",",$numbers);

	my $total_COSMIC_mutations = 0;
	for ( @cosmic_counts ) {
 		   $total_COSMIC_mutations += $_;
	}

	#A sample without a control
	if ($_  =~ m/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR/){

		$GT = $info[0]; #"Genotype"
		$GQ = $info[1]; #"Genotype Quality"
		$SDP = $info[2]; #"Raw Read Depth as reported by SAMtools"
		$DP = $info[3]; #"Quality Read Depth of bases with Phred score >= 15"
		$RD = $info[4]; #"Depth of reference-supporting bases (reads1)"
		$AD = $info[5]; #"Depth of variant-supporting bases (reads2)"
		$FREQ = $info[6]; #"Variant allele frequency"
		$FREQ=~s/%//g;
		$FREQ=$FREQ/100;

		$PVAL = $info[7]; #"P-value from Fisher's Exact Test"
		$PVAL2 = "NA"; #"Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls" 
		$RBQ = $info[8]; #"Average quality of reference-supporting bases (qual1)"
		$ABQ = $info[9]; #"Average quality of variant-supporting bases (qual2)"
		$RDF = $info[10]; #"Depth of reference-supporting bases on forward strand (reads1plus)"
		$RDR = $info[11]; #"Depth of reference-supporting bases on reverse strand (reads1minus)"
		$ADF = $info[12]; #"Depth of variant-supporting bases on forward strand (reads2plus)"
		$ADR = $info[13]; #"Depth of variant-supporting bases on reverse strand (reads2minus)"

		$ADP = $zygosity[0]; #"Average per-sample depth of bases with Phred score >= 15";
		$WT = $zygosity[1]; #"Number of samples called reference (wild-type)";

		if($WT=~ m/WT=0/){
			$WT="Unknown";
		}

	}

	# A tumor sample with a matched control
	else{
		$GT = $info[0]; #"Genotype"
		$GQ = $info[1]; #"Genotype Quality"
		$SDP = "NA"; #"Raw Read Depth as reported by SAMtools"
		$DP = $info[2]; #"Quality Read Depth of bases with Phred score >= 15"
		$RD = $info[3]; #"Depth of reference-supporting bases (reads1)"
		$AD = $info[4]; #"Depth of variant-supporting bases (reads2)"
		$FREQ = $info[5]; #"Variant allele frequency"
		$FREQ=~s/%//g;
		$FREQ=$FREQ/100;
		
		# P value from @zygosity
		$PVAL = $zygosity[3]; #"Fisher's Exact Test P-value of tumor+normal versus no variant for Germline calls"
		$PVAL2 = $zygosity[4]; #"Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls" 
		$PVAL=~s/GPV=//g;
		$PVAL2=~s/SPV=//g;
		$RBQ = "NA"; #"Average quality of reference-supporting bases (qual1)"
		$ABQ = "NA"; #"Average quality of variant-supporting bases (qual2)"
		@fr_read_counts=split(",",$info[6]);
		$RDF = $fr_read_counts[0]; #"Depth of reference-supporting bases on forward strand (reads1plus)"
		$RDR = $fr_read_counts[1]; #"Depth of reference-supporting bases on reverse strand (reads1minus)"
		$ADF = $fr_read_counts[2]; #"Depth of variant-supporting bases on forward strand (reads2plus)"
		$ADR = $fr_read_counts[3]; #"Depth of variant-supporting bases on reverse strand (reads2minus)"
		$ADP = $zygosity[0]; #"Total depth of quality bases";
		$WT = $zygosity[1]; #"Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)";

		if($WT=~ m/SS=0/){
			$WT="Reference";
		}
		elsif($WT=~ m/SS=1/){
			$WT="Germline";
		}
		elsif($WT=~ m/SS=2/){
			$WT="Somatic";
		}
		elsif($WT=~ m/SS=3/){
			$WT="LOH";
		}
		elsif($WT=~ m/SS=4/){
			$WT="Unknown";
		}
		else{
			$WT="Unknown";
		}
	}


	my@output=join("	",$Chr,$Start,	
 	$End,	
 	$Ref,	
	$Alt,	
 	$Func_refGene,	
	$Gene_refGene,	
 	$ExonicFunc_refGene,	
 	$AAChange_refGene,	
	$esp6500si_all,	
	$o1000g2012apr_all,	
	$snp137,
	$total_COSMIC_mutations,
	$cosmic67,
	$LJB2_SIFT,	
 	$LJB2_PolyPhen2_HDIV,		

	$ADP,
 	$WT,
	$GT,
	$GQ,
 	$SDP,
 	$DP,
 	$RD,
 	$AD,
 	$FREQ,
 	$PVAL,
 	$PVAL2,
 	$RBQ,
 	$ABQ,
 	$RDF,
 	$RDR,
 	$ADF,
 	$ADR);	

	if(($intron_exon =~ /intronic/)||($intron_exon =~ /intergenic/)||($Chr =~ /Chr/)){
	}
	else{
	print "@output\n";
		}

	}


