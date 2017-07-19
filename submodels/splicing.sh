#!/bin/bash
varID=$1
prime=$2 # 2nd argument must be either 3 or 5, mapping to score3.pl and score5.pl, respectively.
if [ ! -d "splicing" ]; then
	mkdir splicing
fi
if [ -f *.fa ]; then
	rm *.fa
fi
chr=`echo $varID | cut -d _ -f 2`
site=`echo $varID | cut -d _ -f 3`
site_to=`echo $varID | cut -d _ -f 4`
ref=`echo $varID | cut -d _ -f 5`
alt=`echo $varID | cut -d _ -f 6` 
if [ "$prime" -eq "3" ]; then
	let maxEntLen=23
else
	let maxEntLen=9
fi
let iEnd=maxEntLen-1
#### my $i indicates the in-sequence index for mutation interruption.
if [ "$ref" != "-" ]; then ### for SNP & deletion ###  i: mutation happen after i nucleotides, i.e., residing at i+1 position).
	for i in `seq 0 $iEnd`; do #totally 23 different sequences will be generated for SNP; 22 sequences for deletion - see continue line.
		let chrStart=site-i-1 # -1 is because twoBitToFa is 0-based instead of 1 based. 
		let chrEnd=chrStart+iEnd+1 # +1 is because twoBitToFa requires a non-inclusive end site.
		if [ "$alt" = "-" ]; then
			if [ "$i" -eq "0" ]; then 
				continue 
			fi
			let chrEnd=chrEnd+site_to-site+1
		fi
		twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit -udcDir=. fwd0.fa -seq=chr$chr -start=$chrStart -end=$chrEnd
		longWdt=`tail -1 fwd0.fa` # possibly longer wildtype string 
		echo ${longWdt:0:$maxEntLen}>> fwd.fa
		############## Consider sequence pair after mutation ####################
		let altRelSite=i+1
		if [ "$alt" != "-" ]; then # snp
			sed "2s/[A-Za-z]/$alt/$altRelSite" fwd0.fa | tail -1 >> mut_fwd.fa
		else # deletion
			let j=maxEntLen-i
			Str=`tail -1 fwd0.fa`
			leftStr=${Str:0:$i} # extract string left to to-be-deleted substr
			rightStr=${Str:(-$j)} # extract j-nt string right to to-be-deleted substr
			echo "$leftStr$rightStr" >> mut_fwd.fa			
		fi
	done
else ### for insertion ###
	insLen=${#alt}
	let iStart=-insLen+1 # e.g., inserting 4 nucleotides, iStart will be -3.
	for i in `seq $iStart $iEnd`; do
		if [ "$i" -le "0" ]; then
			let chrStart=site+1-1 # +1 because we are to extract sequence on the right to site, so starting 1 nt more and on;
					      # -1 because because twoBitToFa is 0-based instead of 1 based. 	
			let j=insLen+i # j: how long a substring,extracted from the right of $alt, will show up in the current sequence window.
			insStr=${alt:(-$j)}
			let k=0 # k denotes after how many bits of the current refGenome window all or part of the $alt will be inserted.		
		else
			let chrStart=site-i # -1 is because twoBitToFa is 0-based instead of 1 based.
			insStr=$alt
			let k=i
		fi
		let chrEnd=chrStart+iEnd+1 # Although we only need at most 22nt from ref Genome, we need 23nt for WDT scoring.
					 # +1 because twoBitToFa requires a non-inclusive end site. 
		twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit -udcDir=. fwd0.fa -seq=chr$chr -start=$chrStart -end=$chrEnd
		tail -1 fwd0.fa >> fwd.fa
		############# Consider sequence pair after mutation ####################
		Str=`tail -1 fwd0.fa`
		leftStr=${Str:0:$k} # the left k-nt string.
		rightStr=${Str:$k} # the right string starting from k+1 position.
		tapedStr="$leftStr$insStr$rightStr"
		truncatedStr=${tapedStr:0:$maxEntLen}
		echo $truncatedStr >> mut_fwd.fa
	done
fi
nRows=`cat fwd.fa | wc -l`
for k in `seq 1 $nRows`; do 
	row=`head -$k fwd.fa | tail -1`
	/scratch/cqs/udn/pipeline/submodels/revcompStr.sh $row >> rev.fa
    row=`head -$k mut_fwd.fa | tail -1`
    /scratch/cqs/udn/pipeline/submodels/revcompStr.sh $row >> mut_rev.fa
done
`/scratch/cqs/udn/pipeline/submodels/score$prime.pl fwd.fa >fwd.score`
`/scratch/cqs/udn/pipeline/submodels/score$prime.pl mut_fwd.fa >mut_fwd.score`
`/scratch/cqs/udn/pipeline/submodels/score$prime.pl rev.fa >rev.score`
`/scratch/cqs/udn/pipeline/submodels/score$prime.pl mut_rev.fa >mut_rev.score`
fwd=`awk -v max=-100 '{cnt++; if($2>max){max=$2;idx=cnt}} END {print max "_" idx}' fwd.score`
mut_fwd=`awk -v max=-100 '{cnt++; if($2>max){max=$2;idx=cnt}} END {print max "_" idx}' mut_fwd.score`
rev=`awk -v max=-100 '{cnt++; if($2>max){max=$2;idx=cnt}} END {print max "_" idx}' rev.score`
mut_rev=`awk -v max=-100 '{cnt++; if($2>max){max=$2;idx=cnt}} END {print max "_" idx}' mut_rev.score`
echo -e "fwd_seq.$prime\tfwd_score.$prime\tmutFwd_seq.$prime\tmutFwd_score.$prime\trev_seq.$prime\trev_score.$prime\tmutRev_seq.$prime\tmutRev_score.$prime" > splicing/$varID.score$prime
paste fwd.score mut_fwd.score rev.score mut_rev.score >> splicing/$varID.score$prime
rm *.fa
rm *.score
echo -e "$fwd\t$mut_fwd\t$rev\t$mut_rev"  ## four scores: wildtype forward, mutated forward, wildtype reverse-complement, mutated reverse-complement 
