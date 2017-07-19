#!/usr/bin/perl
### ./srrandRatio.pl UDN829699_step1.txt [0.3]
### this script must be used after step 1 and work as a step2 filter. ideally, the input file name should be UDN###_step1.txt
### it reads in an input file generated with '--withzyg --includeinfo' argument, and writes an output file of only 6 beginning columns (the 6th col for GT info: het/hom).

  
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
	my $refReads;
	if (/\t(\S+)\t(\S+)$/x) {
                @fmt=split /:/,$1;
                @fmtData=split /:/,$2;
                undef(%GT);
                my %GT;
                for ($j=0;$j<=(scalar(@fmt)-1);$j++) {
                        $GT{$fmt[$j]} = $fmtData[$j];
                };
                $varFrac=0; # or called 'strand ratio'
		my @ads;
                if (  exists($GT{'AD'}) and exists($GT{'DP'}) and ($GT{'DP'}>0) ) {
                        @ads = split /,/,$GT{'AD'};
			$refReads = $ads[0];
			shift @ads;
			if (scalar(@ads)>1) {
				$AD=max(@ads); #when multi-allele exists, assume the allele with largest AD is the concerned allele.
			} else {
				$AD=$ads[0];
			}; 
                        $varFrac = $AD/$GT{'DP'};
                };
                my $coverage; 
                my $GQ;
                if ($varFrac>=$th) {
                        if ( exists($GT{'DP'}) and defined($refReads) and defined(@ads) ) {
				$strAds = join(',',@ads);
                                $coverage=$GT{'DP'}.':'.$refReads.':'.$strAds;
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
                                print FO "$former$coverage\t$GQ\n";
                        }

                };
	};
};
close FI;
close FO;
