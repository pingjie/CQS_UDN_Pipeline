options(echo=F)
### UPDATE 5/24/2017: changed argument order. freqFile was put as the 4th (last) argument.
### UPDATE 4/20/2017: set n_Homo threshold to 5 (Filter out variants with n_Homo >=5)
### exac_homo_check(): Referring to downloaded ExAC file, exclude variants according to ExAC freq & num_of_homo.
exac_freqFilter <- function(previous,th=0.001,n_Homo=5,freqFile,step='4.3') { # Only refer to Num_of_homo in the exac file. But still call it a freqFile.
	th <- as.numeric(th)
	n_Homo <- as.numeric(n_Homo)
        proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
        proID <- substr(previous,proID.start,proID.start+8)
        prev <- read.table(previous,header=T,as.is=T,sep='\t',quote='')
        freq <- read.table(freqFile,header=T,as.is=T,sep='\t',quote='')
	freq <- freq[,c('varID','exac_homo','exac_freq')] # will only apply rule on exac_homo, but keep associated exac_freq as a side reference.
	curr <- merge(prev,freq,by.x='varID',by.y='varID',all.x=T)	
	curr <- curr[curr$exac_homo<n_Homo|is.na(curr$exac_homo),]
	curr <- curr[curr$exac_freq<=th|is.na(curr$exac_freq),]
        write.table(curr,paste0(proID,'_step',step,'.txt'),row.names=F,sep='\t',quote=F)
}

args=(commandArgs(TRUE))
previous=args[[1]]
th=args[[2]]
n_Homo=args[[3]]
freqFile=args[[4]]


cmdTxt <- paste0('exac_freqFilter("',previous,'",',th,',',n_Homo,',"',freqFile,'")')
cat('COMMAND: ',cmdTxt,'\n')
eval(parse(text=cmdTxt))

