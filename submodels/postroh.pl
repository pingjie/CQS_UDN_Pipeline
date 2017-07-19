#!/usr/bin/perl

use strict;
use warnings;
use Carp;

print "Chrom\tStart\tEnd\tState\tMeanQual\trunLen_Ms\n";

my (%beg,%prev,%cur,$qual,$nqual);
while (my $line=<STDIN>)
{
    if ($line=~/^#/) {
	next;
    };
    @cur{qw(chr pos state qual)} = split(/\s+/,$line);
    if ( !%beg )
    {
        %beg   = %cur; 
        %prev  = %cur; 
        $qual += $cur{qual};
        $nqual++;
        next; 
    }
    if ( $beg{chr} eq $cur{chr} && $beg{state} eq $cur{state}) 
    { 
        %prev  = %cur;
        $qual += $cur{qual};
        $nqual++;
        next; 
    }
    $qual /= $nqual;
    my $runLen=($prev{pos}-$beg{pos})/1000000;	
    printf "$beg{chr}\t$beg{pos}\t$prev{pos}\t$beg{state}\t%.1f\t%.1f\n",$qual,$runLen;
    %beg   = %cur;
    $qual  = $cur{qual};
    $nqual = 1;
}
if ( %beg && %cur )
{
    $qual /= $nqual;
    my $runLen = ($prev{pos}-$beg{pos})/1000000; 
    printf "$beg{chr}\t$beg{pos}\t$prev{pos}\t$beg{state}\t%.1f\t%.1f\n",$qual,$runLen;
}

