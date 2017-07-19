options(echo=F)
### collateExacMetrics(previous,knowFile): from knowledge file knowFile we extract missense Z score (mis_z) and probability of loss-of-function intolorence (pLI) and collate them to the variants included in previous file (step3 file). 
collateExacMetrics <- function(previous,knowFile,step='3') {
        proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
        proID <- substr(previous,proID.start,proID.start+8)
	prev <- read.table(previous,header=T,as.is=T,sep='\t',quote='')
	knowTbl <- read.delim(knowFile,as.is=T)
	knowTbl <- knowTbl[,c('gene','mis_z','pLI')]
	curr <- merge(prev,knowTbl,by.x='Gene.refGene',by.y='gene',all.x=T)
	curr <- data.frame(curr[,c(2:7,1)],curr[,-(1:7)])
	write.table(curr,paste0(proID,'_step',step,'.txt'),row.names=F,sep='\t',quote=F)
	
}
args=(commandArgs(TRUE))
previous=args[[1]]
knowFile=args[[2]]
cmdTxt <- paste0('collateExacMetrics("',previous,'","',knowFile,'")')
cat('COMMAND: ',cmdTxt,'\n')
eval(parse(text=cmdTxt))
