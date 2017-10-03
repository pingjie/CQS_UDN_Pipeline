#!/bin/bash -i

##### Beginning of Configuration

SeqPlatform=	### 1 Baylor, 2 Hudson Alpha

ProBand_ID=
ProBand_Gender=	### Male 1, Female 2

ParentIDs=() ### use space to separate different IDs ### If a parent is a nofilter ID, do not include it here.
ParentGenders=()  ### use space to separate different Genders
ParentAffStat=()   ### use space to separate different affected status

BroSisIDs=()   ### use space to separate different IDs
BroSisGenders=()    ### use space to separate different Genders
BroSisAffStat=()    ### use space to separate different affected status

nofilterIDs=()  ### for testing purpose only ### allow one or more nofilter members.
nofilterGender=()
nofilterAffStat=() ### since the affect status is unambiguous, let's set it as 1.

RAW_VCF_DIR=	### Raw VCF directory, you'd better put all $UDN**.vcf under the same folder 
WORK_DIR=	### Working directory
PhenoFile=

### Parameters

StrandRatioCutOff=0.3
th=0.01
n_Homo=5
thStrict=1E-4

##### End of Configuration

### make several directories

PREPROCESS_DIR=$WORK_DIR/0_preprocessed_steps
FILTER_DIR=$WORK_DIR/1_filtered_steps
IH_DIR=$WORK_DIR/2_ihFilter

mkdir $PREPROCESS_DIR $FILTER_DIR $IH_DIR

mkdir $IH_DIR\/prioritization

### Software Configurations & auxiliary files 

setpkgs -a R_3.2.2_gcc
setpkgs -e java_1.8
export R_LIBS="/home/pingj1/R/x86_64-pc-linux-gnu-library/3.2:$R_LIBS"
export PATH=/home/pingj1/local/bin:$PATH
CODES_DIR=/workspace/pingj1/soft/CQS_UDN_Pipeline/submodels # To be finalized
ANNOVARDBDIR=/scratch/yuh9/software/annovar/humandb/
AFfile=/scratch/cqs/udn/hg19_g1k2015_roh.tab.gz
GNOMADdata=/scratch/cqs/udn/gnomad_freq.rdata
EXAC_metrics_file=/scratch/cqs/udn/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
APPRISfile=/scratch/cqs/udn/apprisP.refSeq.txt

### merge ids ####

FamilyIDs=(${ParentIDs[@]} ${BroSisIDs[@]})
FamilyGenders=(${ParentGenders[@]} ${BroSisGenders[@]})
FamilyAffStatus=(${ParentAffStat[@]} ${BroSisAffStat[@]})
if [[ -n ${nofilterIDs[@]} ]]; then
	FamilyIDs=(${FamilyIDs[@]} ${nofilterIDs[@]})
	FamilyGenders=(${FamilyGenders[@]} ${nofilterGender[@]})
	FamilyAffStatus=(${FamilyAffStatus[@]} ${nofilterAffStat[@]})
fi
### preprocessing VCF files ###
for IDs in $ProBand_ID "${FamilyIDs[@]}"; do
        if [[ -f $RAW_VCF_DIR/$IDs\.vcf ]]; then
                ln -s $RAW_VCF_DIR/$IDs.vcf $PREPROCESS_DIR/$IDs\.vcf
        else
                if [ $SeqPlatform -eq 1 ]; then
                        bgzip -c "$(find $RAW_VCF_DIR -regex ".*\-$IDs\-.*SNPs.*vcf")" > $PREPROCESS_DIR/$IDs\.SNPs.vcf.gz
                        bgzip -c "$(find $RAW_VCF_DIR -regex ".*\-$IDs\-.*INDELs.*vcf")" > $PREPROCESS_DIR/$IDs\.INDELs.vcf.gz
                        tabix -p vcf $PREPROCESS_DIR/$IDs\.SNPs.vcf.gz
                        tabix -p vcf $PREPROCESS_DIR/$IDs\.INDELs.vcf.gz
                        bcftools concat -a $PREPROCESS_DIR/$IDs\.SNPs.vcf.gz $PREPROCESS_DIR/$IDs\.INDELs.vcf.gz -o $PREPROCESS_DIR/$IDs\.vcf
                else
                        ln -s "$(find $RAW_VCF_DIR -regex ".*$IDs.*vcf")" $PREPROCESS_DIR/$IDs\.vcf
                fi
        fi
done

### Filter Probands

cd $FILTER_DIR

convert2annovar.pl --format vcf4 --filter PASS --withzyg --includeinfo $PREPROCESS_DIR/$ProBand_ID\.vcf > $FILTER_DIR/$ProBand_ID\_step1.vcf

if [ $SeqPlatform -eq 1 ]; then
	perl $CODES_DIR/strandRatio_bl.pl $FILTER_DIR/$ProBand_ID\_step1.vcf $StrandRatioCutOff
else
	perl $CODES_DIR/strandRatio_ha.pl $FILTER_DIR/$ProBand_ID\_step1.vcf $StrandRatioCutOff
fi

table_annovar.pl $FILTER_DIR/$ProBand_ID\_step2.vcf $ANNOVARDBDIR -out $FILTER_DIR/$ProBand_ID -build hg19 --otherinfo -protocol refGene,snp138,esp6500siv2_all,exac03,1000g2015aug_all,gnomad_exome,gnomad_genome -operation g,f,f,f,f,f,f -arg '-splicing 8,,,,,,' 2> $FILTER_DIR/annovar.out

nCol=`head -1 $FILTER_DIR/$ProBand_ID\.hg19_multianno.txt | wc -w`
let nCol=nCol+4
cut -f 1-$nCol $FILTER_DIR/$ProBand_ID\.hg19_multianno.txt > $FILTER_DIR/temp.annovar
header=`head -1 $FILTER_DIR/temp.annovar`
echo -e "$header\tqual\tdepth\tcoverage\tGQ" > $FILTER_DIR/temp.header
nRow=`cat $FILTER_DIR/temp.annovar | wc -l`
let nRow=nRow-1
tail -$nRow $FILTER_DIR/temp.annovar > $FILTER_DIR/temp.annovar1
cat $FILTER_DIR/temp.header $FILTER_DIR/temp.annovar1 > $FILTER_DIR/$ProBand_ID\.annovar.txt
rm $FILTER_DIR/temp*

if [ $SeqPlatform -eq 1 ]; then
	awk -F "\t" '($6~/splicing/ || ($9!="" && $9!="synonymous SNV" && $9!="unknown")) {print $0}' $FILTER_DIR/$ProBand_ID\.annovar.txt > $FILTER_DIR/$ProBand_ID\_step3.0.txt
else
	awk -F "\t" '($6=="upstream" || $6=="downstream" || $6=="upstream;downstream" || $6=="UTR5" || $6=="UTR3" || $6=="UTR5;UTR3" || $6=="intronic" || $6~/splicing/ || $6=="exonic;splicing" || ($9!="" && $9!="synonymous SNV" && $9!="unknown") ) {print $0}' $FILTER_DIR/$ProBand_ID\.annovar.txt > $FILTER_DIR/$ProBand_ID\_step3.0.txt 
fi

mv $FILTER_DIR\/$ProBand_ID.annovar.txt $WORK_DIR

Rscript $CODES_DIR/collateExacMetrics.R $FILTER_DIR/$ProBand_ID\_step3.0.txt $EXAC_metrics_file
if [ $SeqPlatform -eq 2 ]; then
	Rscript $CODES_DIR/polymorphFilter_ha.R $ProBand_ID\_step3.txt $th $n_Homo
else
	perl $CODES_DIR/extract_fields_from_vcf.pl $PREPROCESS_DIR/$ProBand_ID\.vcf ARIC_AA,ARIC_EA
	mv $FILTER_DIR/$ProBand_ID\_ARIC_AA,ARIC_EA.tab $FILTER_DIR/$ProBand_ID\_aric.txt

	Rscript $CODES_DIR/aric_freqFilter.R $FILTER_DIR/$ProBand_ID\_step3.txt $FILTER_DIR/$ProBand_ID\_aric.txt $th

	perl $CODES_DIR/extract_fields_from_vcf.pl $PREPROCESS_DIR/$ProBand_ID\.vcf IET,IEO
	mv $FILTER_DIR/$ProBand_ID\_IET,IEO.tab $FILTER_DIR/$ProBand_ID\_ESEindel.txt

	Rscript $CODES_DIR/eseIndel_freqFilter.R $ProBand_ID\_step4.1.txt $ProBand_ID\_ESEindel.txt $th

	Rscript $CODES_DIR/gnomad_freqFilter.R $ProBand_ID\_step4.2.txt $th $n_Homo $GNOMADdata

	Rscript $CODES_DIR/bulk4_freqFilter.R $ProBand_ID\_step4.3.txt $th
fi

### ROH for Proband

bcftools roh -G30 -I --AF-file $AFfile $PREPROCESS_DIR/$ProBand_ID\.vcf | perl $CODES_DIR/postroh.pl > $FILTER_DIR/$ProBand_ID\.roh
bcftools query -f'%CHROM\t%POS\t[%GT]\n' $PREPROCESS_DIR/$ProBand_ID\.vcf | awk '$3~/^[01\/]+$/' | sed 's,/,\t,' | awk -F "\t" '{print $1,$2,$3+$4}' > $FILTER_DIR/$ProBand_ID\.dos 
Rscript $CODES_DIR/roh.R $ProBand_ID 1

### Filter Family VCFs

echo "Preparing step1 & step2 for all family members"
for (( i=0; i<${#FamilyIDs[@]}; i++ ))
do
	bcftools roh -G30 -I --AF-file $AFfile $PREPROCESS_DIR/${FamilyIDs[i]}\.vcf | perl $CODES_DIR/postroh.pl > $FILTER_DIR/${FamilyIDs[i]}\.roh
	bcftools query -f'%CHROM\t%POS\t[%GT]\n' $PREPROCESS_DIR/${FamilyIDs[i]}\.vcf | awk '$3~/^[01\/]+$/' | sed 's,/,\t,' | awk -F "\t" '{print $1,$2,$3+$4}' > $FILTER_DIR/${FamilyIDs[i]}\.dos 
	Rscript $CODES_DIR/roh.R ${FamilyIDs[i]} ${FamilyAffStatus[i]}

	if [ "${FamilyAffStatus[i]}" -eq 1 ]; then
		convert2annovar.pl --format vcf4 --filter PASS --withzyg --includeinfo $PREPROCESS_DIR/${FamilyIDs[i]}\.vcf > $FILTER_DIR/${FamilyIDs[i]}\_step1.txt

		if [ $SeqPlatform -eq 1 ]; then
			perl $CODES_DIR/strandRatio_bl.pl $FILTER_DIR/${FamilyIDs[i]}\_step1.txt $StrandRatioCutOff
		else
			perl $CODES_DIR/strandRatio_ha.pl $FILTER_DIR/${FamilyIDs[i]}\_step1.txt $StrandRatioCutOff
		fi
	else
		convert2annovar.pl --format vcf4 --withzyg $PREPROCESS_DIR/${FamilyIDs[i]}\.vcf > $FILTER_DIR/${FamilyIDs[i]}\_step0.txt
	fi
done

Rscript $CODES_DIR/roh_final.R $ProBand_ID

### Step 5 ####

cd $IH_DIR

find $FILTER_DIR -iregex ".*step.\.txt" -exec cp {} $IH_DIR \;
Rscript $CODES_DIR/appris.R $ProBand_ID\_step4.txt $APPRISfile
if [ "$SeqPlatform" = "1" ]; then
	ihFilter=ihFilter_bl.R
else
	ihFilter=ihFilter_ha.R
fi
if [ ${#ParentAffStat[@]} -eq 2 ] && [ $((${ParentAffStat[0]} + ${ParentAffStat[1]})) -eq 0 ]; then
	Rscript $CODES_DIR/$ihFilter parTrio $ProBand_ID\_step4.txt ${ParentIDs[0]} ${ParentIDs[1]} thStrict=$thStrict Gender=$ProBand_Gender
	cp $ProBand_ID\_parFiltered.txt $ProBand_ID\_tmpFiltered.txt
	mv $ProBand_ID\_parFiltered.txt $ProBand_ID\_0_${ParentIDs[0]}\_${ParentIDs[1]}_parFiltered.txt
	mv parFiltered.txt 0_${ParentIDs[0]}\_${ParentIDs[1]}_parFiltered.txt
	unset FamilyIDs[0]
	unset FamilyIDs[1]
	unset FamilyAffStatus[0]
	unset FamilyAffStatus[1]
	FamilyIDs=("${FamilyIDs[@]}")
	FamilyAffStatus=("${FamilyAffStatus[@]}")
else
	Rscript $CODES_DIR/$ihFilter impute $ProBand_ID\_step4.txt $thStrict
	cp $ProBand_ID\_parFiltered.txt $ProBand_ID\_tmpFiltered.txt
	mv $ProBand_ID\_parFiltered.txt $ProBand_ID\_imputeFiltered.txt
fi

for (( i=0; i<${#FamilyIDs[@]}; i++ ))
do
	let j=i+1
	if [[ ${FamilyAffStatus[i]} -eq 1 ]]; then
		Rscript $CODES_DIR/$ihFilter common $ProBand_ID\_tmpFiltered.txt ${FamilyIDs[i]}
		cp $ProBand_ID\_commonFiltered.txt $ProBand_ID\_tmpFiltered.txt
		cp commonFiltered.txt tmpFiltered.txt
		mv $ProBand_ID\_commonFiltered.txt $ProBand_ID\_$j\_${FamilyIDs[i]}_commonFiltered.txt
		mv commonFiltered.txt $j\_${FamilyIDs[i]}_commonFiltered.txt
	else
		Rscript $CODES_DIR/$ihFilter unique $ProBand_ID\_tmpFiltered.txt ${FamilyIDs[i]}
		cp $ProBand_ID\_uniqFiltered.txt $ProBand_ID\_tmpFiltered.txt
		cp uniqFiltered.txt tmpFiltered.txt
		mv $ProBand_ID\_uniqFiltered.txt $ProBand_ID\_$j\_${FamilyIDs[i]}_uniqFiltered.txt
		mv uniqFiltered.txt $j\_${FamilyIDs[i]}_uniqFiltered.txt
	fi
	
	if [[ ${#nofilterIDs[@]} -gt 0 ]] && [[ " ${nofilterIDs[@]} " == *" ${FamilyIDs[i]} "* ]]; then
		tail -n +2 $IH_DIR\/tmpFiltered.txt > $IH_DIR\/salvaged.txt
		cp $IH_DIR\/$ProBand_ID\_tmpFiltered.txt $IH_DIR\/$ProBand_ID\_tmpFiltered1.txt
		cat $IH_DIR\/$ProBand_ID\_tmpFiltered1.txt $IH_DIR\/salvaged.txt > $IH_DIR\/$ProBand_ID\_tmpFiltered.txt
		rm $IH_DIR\/$ProBand_ID\_tmpFiltered1.txt *${FamilyIDs[i]}*Filtered.txt $IH_DIR\/salvaged.txt
	fi
	
done
mv $IH_DIR\/$ProBand_ID\_tmpFiltered.txt $IH_DIR\/$ProBand_ID\_finalFiltered.txt

############ Phenolyzer prioritization ########

cd $IH_DIR

tail -n +2 $IH_DIR\/$ProBand_ID\_finalFiltered.txt | cut -f 3 | sort | uniq > $IH_DIR\/prioritization\/$ProBand_ID\_candGenes.txt
sed -i '/,/s/,.*//' $IH_DIR\/prioritization\/$ProBand_ID\_candGenes.txt
cp $PhenoFile $IH_DIR\/prioritization\/$ProBand_ID\_terms.txt

/scratch/yuh9/software/phenolyzer/disease_annotation.pl prioritization\/$ProBand_ID\_terms.txt -f -p -ph -logistic --gene prioritization\/$ProBand_ID\_candGenes.txt -out prioritization\/$ProBand_ID 2> prioritization\/phenolyzer.log
Rscript $CODES_DIR/prio_cmm.R $IH_DIR\/prioritization\/$ProBand_ID $IH_DIR\/$ProBand_ID

### splicing scores ####
while IFS='' read line || [[ -n "$line" ]]; do
        # Each line is read into variable $line, and cells in $line are separated by spaces.
        varID=`echo $line | cut -d " " -f 2`
        genomicFunc=`echo $line | cut -d " " -f 3`
	if [ "$varID" = "varID" ]; then
		echo -e "$line\tfwd.ss3\tmut_fwd.ss3\trev.ss3\tmut_rev.ss3\tfwd.ss5\tmut_fwd.ss5\trev.ss5\tmut_rev.ss5" > $WORK_DIR/$ProBand_ID\_spliceScored.txt
		continue
	fi
	if [ "$genomicFunc" != "exonic" ]; then
		score5=`$CODES_DIR/splicing.sh $varID 5`
		score3=`$CODES_DIR/splicing.sh $varID 3`
		score8="$score3\t$score5"
		paste splicing\/$varID.score3 splicing\/$varID.score5 > splicing\/$varID.score
	else
		score8="\t\t\t\t\t\t\t"
	fi
	echo -e "$line\t$score8" >> $WORK_DIR/$ProBand_ID\_spliceScored.txt
done < $IH_DIR\/$ProBand_ID\_prioCmmed.txt

### Generate xlxs table

Rscript $CODES_DIR\/convertFormat.R $WORK_DIR\/$ProBand_ID\_spliceScored.txt $SeqPlatform

echo "JOB DONE!"

