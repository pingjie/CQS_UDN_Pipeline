options(echo=F)
### UPDATED 6/10/2016 for accomodating ANNOVAR-formated step3 previous file and new variant ID format.
eseIndel_freqFilter <- function(previous,freqFile,cntTH=0,step=4.2) {
	cntTH <- as.numeric(cntTH)
        proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
        proID <- substr(previous,proID.start,proID.start+8)
        prev <- read.table(previous,header=T,as.is=T,sep='\t',quote='')
        freq <- read.table(freqFile,header=T,as.is=T,sep='\t',quote='')
        freq <- freq[,c('varID','IET','IEO')]
	freq$IET <- as.numeric(freq$IET)
	freq$IEO <- as.numeric(freq$IEO)
        freq[is.na(freq)] <- 0
        curr <- merge(prev,freq,by.x='varID',by.y='varID',all.x=T)
        curr <- curr[curr$IET<=cntTH & curr$IEO<=cntTH,] # these two columns seem to be always filled with a number, say 0.
        write.table(curr,paste0(proID,'_step',step,'.txt'),row.names=F,sep='\t',quote=F)

}

args=(commandArgs(TRUE))
previous=args[[1]]
freqFile=args[[2]]
cntTH=args[[3]]
cmdTxt <- paste0('eseIndel_freqFilter("',previous,'","',freqFile,'",',cntTH,')')
cat('COMMAND: ',cmdTxt,'\n')
eval(parse(text=cmdTxt))
