#!/usr/bin/perl
### ./strandRatio_bl.pl UDN829699_step1.txt [0.3]
### this script must be used after step 1 and work as a step2 filter. Input file must be of a VCF format.
### INPUT & OUTPUT naming: input file name should contain *step1* and the output will change accordingly to *step2*
### it reads in an input file generated with '--withzyg --includeinfo' argument, and writes an output file with ANNOVAR-converted GT info, depth, qual, coverage, and GQ.

  
use List::Util qw(min max);
open FI,$ARGV[0] || die "Cannot open input file!!!\n";
$outFile=$ARGV[0];
$outFile=~s/step1/step2/; ### this script must be used after step1 and to lead to step2 result.
open FO,">$outFile" || die "Cannot write outfile!!!\n";

if (defined($ARGV[1])) {
	$th=$ARGV[1];
} else {
	$th=0.15;
};
while (<FI>) {
	if (/\t(\S+)\t (\S+)$/x) { # 9th & 10th columns
		#print "inside each line format judgement loop\n";
        @fmt=split /:/,$1;
        @fmtData=split /:/,$2;
        undef(%GT);
        my %GT;
        for ($j=0;$j<=(scalar(@fmt)-1);$j++) {
        	$GT{$fmt[$j]} = $fmtData[$j];
        };
        $varFrac=0; # or called 'strand ratio'
     	$AD=0;
	    if (  exists($GT{'VR'}) and exists($GT{'DP'}) and ($GT{'DP'}>0) ) { # VR is equivalent to AD of Hudson Alpha.
			#print "inside varFrac judgement loop\n";
			$AD=$GT{'VR'}; ## VR field has only 1 integer. No multiple Alleles, so addressing way is simpler than HA.
            $varFrac = $AD/$GT{'DP'};
        };
		my $coverage;
		my $GQ;
        if ($varFrac>=$th) {
			if ( exists($GT{'DP'}) and exists($GT{'RR'}) and exists($GT{'VR'}) ) {
				$coverage=$GT{'DP'}.':'.$GT{'RR'}.':'.$GT{'VR'};
			} else {
				$coverage='NA';
			}
			if ( exists($GT{'GQ'}) ) {
				$GQ=$GT{'GQ'};
			} else {
				$GQ='NA';
			}
			chomp($_);
			if ($_ =~ /^((\S+\s){8})(.*)$/) { 
				$former=$1;
				$latter=$3;
                print FO "$former$coverage\t$GQ\t$latter\n";
			} 
        };
	};
};
close FI;
close FO;
