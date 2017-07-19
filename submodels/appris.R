appris <- function(step4file,apprisList='/gpfs21/scratch/cqs/udn/apprisP.refSeq.txt') {
	step4 <- read.delim(step4file,head=T,as.is=T)
	apprisList <- scan(apprisList,'')
	GeneDetail <- step4$GeneDetail.refGene
	AAChange <- step4$AAChange.refGene
	step4$GeneDetail.refGene <- sapply(GeneDetail,getAppris,apprisList,'gene')
	step4$AAChange.refGene <- sapply(AAChange,getAppris,apprisList,'AAChange')	
        cat(sum(nchar(step4$GeneDetail.refGene)==0 & nchar(step4$AAChange.refGene)==0),'out of',nrow(step4),'rows were deleted due to non-APPRIS principals.\n')
	step4 <- step4[nchar(step4$GeneDetail.refGene)!=0 | nchar(step4$AAChange.refGene)!=0,]
	write.table(step4,step4file,row.names=F,sep='\t',quote=F)
}
getAppris <- function(multiItems,apprisList,type=c('gene','AAChange')[2]) {
	multiVec <- unlist(strsplit(multiItems,','))
	if (length(multiVec)==1)
		multiVec <- unlist(strsplit(multiItems,';'))
	ix <- switch(type,
		gene=1,
		AAChange=2
	)
	nmIDs <- sapply(multiVec,function(x) {unlist(strsplit(x,':'))[ix]} )
	inList.ix <- which(nmIDs%in%apprisList)
	if (length(inList.ix)>0) {
		multiVec <- multiVec[inList.ix]
		newString <- paste(multiVec,collapse=',')
	} else {
		newString <- ''
	}
	newString
}
args <- commandArgs(TRUE)
appris(args[1],args[2])
