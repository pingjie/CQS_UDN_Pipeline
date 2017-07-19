#!/usr/bin/perl
# extract_fields_from_vcf.pl LB_UDN813101_step1.vcf ENS,RFG,UCG
### UPDATE 6/9/2016: changing varID coding in order to be compatible with ANNOVAR-derived varID.

$vcfFile=$ARGV[0];
$fields=$ARGV[1];
@several=split(',',$fields);
$nFields=scalar(@several);

if ($vcfFile=~/(UDN\d+)/) {
        $outFile=$1."_".$fields.".tab"; # as the result file is largely distant to VCF, the output file is not named to vcf2txt,but to *.tab.
};
open vcfFI,$vcfFile || die "Cannot open vcf\n";
open fieldFO,">$outFile" || die "Cannot write outfile";
$tabFields=join("\t",@several);
print fieldFO "chr\tchrStart\tchrEnd\tref\talt\tqual\tfilter\t$tabFields\tvarID\n";
while (<vcfFI>) {
	if (/^#/) {
		next;
	};
	if (/^(\S+)\t  (\S+)\t  (\S+)\t  (\S+)\t (\S+)\t (\S+)\t (\S+)\t  (\S+)\t/x) {
		$annoPrefix="$1\t$2\t$2\t$4\t$5\t$6\t$7";
		$chr=$1;
		$oStart=$2;
		$oRef=$4;
		$oAlt=$5;
		if ( length($oRef)==length($oAlt)& length($oRef)==1 ) {
			$start=$oStart;
			$end=$start;
			$ref=$oRef;
			$alt=$oAlt;
		} elsif ( length($oRef)<length($oAlt)& length($oRef)==1 ) { #pure insertion
			$start=$oStart;
			$end=$oStart;
			$ref='-';
			$alt=substr($oAlt,1,length($oAlt)-1);
		} elsif ( length($oRef)>length($oAlt)& length($oAlt)==1 ) { #pure deletion
			$start=$oStart+1;
			$end=$oStart+length($oRef)-1;
			$ref=substr($oRef,1,length($oRef)-1);
			$alt='-';
		} else {
			die 'Peculiar substitution happens!\n';
		};
		$varID="var_$chr\_$start\_$end\_$ref\_$alt";
		$INFO=$8;
		undef($string);
		for ($i=0;$i<=($nFields-1);$i++) {
			$field=$several[$i];
			undef($value);
			if ($INFO=~/;$field\=([^;]+)/) { #$INFO=~/$field\=(\S{10})/     worked for ENSGN
				$value=$1;
			} elsif ($INFO=~/^$field\=([^;]+)/) {
				$value=$1;
			};
			if ($i==0) {
				$string=$value;
			} else {
				$string=$string."\t".$value;
			};
		};
		print fieldFO "$annoPrefix\t$string\t$varID\n";
	};
};

close vcfFI;
close fieldFO;


