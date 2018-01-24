options(java.parameters = "-Xmx8g")
#Sys.setlocale(locale="C")
if (!require('xlsx')) {
	stop('The package xlsx was not installed')
}

is.empty <- function(x, mode=NULL){
  if (is.null(mode)) mode <- class(x)
  identical(vector(mode,1),c(x,vector(class(x),1)))
}

library(data.table)
library(stringr)

geneomim <- fread("/scratch/cqs/udn/genemap2.txt")

processline<-function(x){
	loc<-unlist(strsplit(x$varID,"_"))
	chr<-paste("chr",loc[2],sep="")
	pos<-loc[3]
	change<-paste(loc[5],loc[6],sep=" → ")
	rsID<-x$snp138
	gnomad <- x$gnomadAF #max(x$gnomADexomeFreq, x$gnomADgenomeFreq)
	gnomadHomo <- max(c(x$gnomadHom,x$gnomadHemi),na.rm=T)
	if (gnomadHomo>0) {
		gnomad=paste(gnomad, paste(gnomadHomo,'homoz'))
		gnomad=gsub('NA ','',gnomad)
	}
	genename<-x$geneName
	proGT<-ifelse(x$GT=="het","●○","●●")
	udnid<-grep("UDN",colnames(x))
	udnGT<-sapply(x[udnid],function(x){switch(as.character(x),wdt="○○",het="●○",hom="●●")}) #homo to hom: UPDATED by Hui 05/28/2017 ##some don't have genotype

	udnGT[unlist(lapply(udnGT,is.null))]<-NA
	comment <- x$comment
	for (n in 1:length(comment)) {
  		if( !is.empty(which(geneomim$Approved_Symbol == genename[n])) ) {
    		if ( geneomim$Phenotypes[which(geneomim$Approved_Symbol == genename[n])] != "") {
      			phenotypeOMIM <- geneomim$Phenotypes[which(geneomim$Approved_Symbol == genename[n])]
      			tmp <- strsplit(phenotypeOMIM, ", ")[[1]]
      			for (t in 1:length(tmp)) {
        			if (str_detect(tmp[t], "\\d+ \\(\\d\\)")) {
          				tmp[t] <- paste0("MIM:", strsplit(tmp[t], " ")[[1]][1])
        			}
      			}
      			phenotypeOMIM <- paste(tmp, collapse = ", ")
      
      			comment[n] <- paste0(phenotypeOMIM, " || ", comment[n])
    		}
  		}
	}

	rank <- x$Rank
	if (x$Func.refGene =="exonic") {
		effect<-x$ExonicFunc.refGene 
		effect[effect=="nonsynonymous SNV"]<-"missense"
		tmp <- unlist(strsplit(x$AAChange.refGene,","))[1]
		detail<-unlist(strsplit(tmp,':',fixed=T))
		dnachange<-detail[grep("c.",detail,fixed=T)][1]
		pchange<-detail[grep("p.",detail)][1]
		refseq<-detail[2]
	} else if (x$Func.refGene =="splicing" | x$Func.refGene =="UTR5" | x$Func.refGene =="UTR3" | x$Func.refGene =="ncRNA_splicing") {
		effect<-x$Func.refGene
		tmp <- unlist(strsplit(x$GeneDetail.refGene,","))[1]
		detail<-unlist(strsplit(tmp,':',fixed=T))
		refseq<-detail[1]
		dnachange<-detail[grep("c.",detail,fixed=T)][1]
		pchange<-NA
	} else {
		effect<-x$Func.refGene
		refseq<-""
		dnachange<-NA
		pchange<-NA
	}
	row1<-c(paste(genename,refseq,sep="\r\n"),chr,change,as.character(effect),proGT,udnGT,x$qual,round(x$g1kFreq,digits=4),round(x$mis_z,digits=2),rep(NA,3),rank,comment)
	row2<-c(NA,pos,dnachange,NA,NA,rep(NA,length(udnGT)),x$GQ,x$espFreq,round(x$pLI,digits=2),rep(NA,3),'✓',NA)
	row3<-c(NA,rsID,pchange,NA,NA,rep(NA,length(udnGT)),x$coverage,x$exacFreq,gnomad,rep(NA,5))
	combinerow<-rbind(row1,row2,row3)
	return(combinerow)
}


processMode<-function(data,ihMode){
	modedata<-data[data$ihMode==ihMode,]
	result<-NULL
	if (dim(modedata)[1]>0){
		for (i in 1:dim(modedata)[1]){
			result<-rbind(result,processline(modedata[i,]))
		}
	}
	return(result)
}

mergecells<-function(sheet,rowvec,colvec){
	for (j in 1:(length(rowvec)-1)){
		for (i in 1:length(colvec)){
			addMergedRegion(sheet,rowvec[j],rowvec[j]+2,colvec[i],colvec[i])
		}
	}
}

#ihType=c('De Novo Variants','Compound Heterozygous Variants','Homozygous Variants','X-Linked Variants')
writeModeBlock <- function(ihType,variants,sheet,startrow,template,mergecol,mergecol2,cs1,cs2,cs3) {
	len=ncol(template)
	ihtype<-matrix(c(ihType, rep(NA,len-1)),nrow=1)
	cb1<-CellBlock(sheet,startrow,1,1,len)
	CB.setMatrixData(cb1,ihtype,1,1,cellStyle=cs1)
	addMergedRegion(sheet,startrow,startrow,1,len)
	startrow<-startrow+1
	cb2<-CellBlock(sheet,startrow,1,3,len)
	CB.setMatrixData(cb2,template,1,1,cellStyle=cs2)
	mergecells(sheet,startrow:(startrow+1),mergecol)
	startrow<-startrow+3
	#variants<-processMode(data,"ARs") 
	if (!is.null(variants)){
		cb3<-CellBlock(sheet,startrow,1,dim(variants)[1],len)
		CB.setMatrixData(cb3,variants,1,1,cellStyle=cs3)
		mergecells(sheet,seq(startrow,startrow+dim(variants)[1],3),mergecol2)
		startrow<-startrow+dim(variants)[1]
	}
	list(sheet = sheet, startrow = startrow)
}

#INPUT platform: 1 for Baylor, 2 for Hudson Alpha
convertFormat <- function(filename, platform) {
	platform <- as.character(platform)
	data<-read.table(filename,header=T,as.is=T,fill=T,sep="\t", quote="")
	data <- data[order(data$geneName),] ## UPDATED by Hui 06/02/2017
	outputfile<-gsub(".txt",".formatted.xlsx",filename)
	##remove inROH#####
	only.inROH <- which(data$ihMode=="inROH ")
	if (length(only.inROH)>0) {
		for (i in only.inROH) {
			data$ihMode[i] <- paste('inROH',ifelse(data$GT[i]=='het','deNovo','ARs'))
		}
	}
	data$ihMode<-gsub("inROH ","",data$ihMode)
	#generate header
	UDNid<-colnames(data)[grep("UDN",colnames(data))]
	platform <- switch(platform, '1'='Baylor', '2'='Hudson Alpha')
	template1<-c("Gene","Chr","Change","Effect","Proband",UDNid,"Quality","1KG AF","Missense Z",platform,"Pheno DB","Omicia","Yu Shyr","Comments")
	template2<-c(NA,"Position","cLevel",NA,NA,rep(NA,length(UDNid)),"GQ","EVS AF","LoF pLI",rep(NA,5))
	template3<-c(NA,"rs#","pLevel",NA,NA,rep(NA,length(UDNid)),"Coverage","ExAC AF","GnomAD",rep(NA,5))
	len<-length(template1)
	template<-rbind(template1,template2,template3)
	mergecol<-c(1:len)[is.na(template2)]
	mergecol2<-mergecol[-(length(mergecol)-1)]
##define xlsx format
	wb<-createWorkbook()
	sheet<-createSheet(wb)
##ihMode type
	cs1<-CellStyle(wb)+Fill(backgroundColor="#d9d9d9",foregroundColor="#d9d9d9")+Font(wb,heightInPoints=20,isBold=TRUE, name="Calibri")
##header ####
	cs2<-CellStyle(wb)+Fill(backgroundColor="#dce6f1",foregroundColor="#dce6f1")+Alignment(h="ALIGN_CENTER",vertical= "VERTICAL_CENTER")+Font(wb,heightInPoints=12,isBold=TRUE,name="Calibri")+ Border(color="#d9d9d9", position=c("LEFT","RIGHT"))
##variants####
	cs3<-CellStyle(wb)+Alignment(h="ALIGN_CENTER",wrapText=T)+Font(wb,heightInPoints=12,isBold=FALSE,name="Calibri")+DataFormat("#,##0.00") ## UPDATED by Hui 06/02/2017 for wrapText
	startrow<-1
	variants <- processMode(data,"deNovo")
	res<-writeModeBlock('De Novo Variants',variants,sheet,startrow,template,mergecol,mergecol2,cs1,cs2,cs3)
	variants1 <- processMode(data,"ARch")
	variants2 <- processMode(data,'ARch.tier1')
	variants3 <- processMode(data,'ARch.tier2')
	variants <- rbind(variants1,variants2,variants3)
	res<-writeModeBlock('Compound Heterozygous Variants',variants,res$sheet,res$startrow,template,mergecol,mergecol2,cs1,cs2,cs3)
	variants <- processMode(data,'ARs')
	res <- writeModeBlock('Homozygous Variants',variants,res$sheet,res$startrow,template,mergecol,mergecol2,cs1,cs2,cs3)
	variants <- processMode(data,'Xlinked')
	res <- writeModeBlock('X-Linked Variants',variants,res$sheet,res$startrow,template,mergecol,mergecol2,cs1,cs2,cs3)
	sheet=res$sheet
	saveWorkbook(wb,outputfile)
}

###Main##########################
args <- commandArgs(TRUE)
convertFormat(args[1], args[2])
