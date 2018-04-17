options(echo=F)
#args <- c('TW_UDN525928_SNP_step3.2.txt','TW_UDN525928_SNP_ET,EV,AT,AV.txt')
### UPDATED 6/9/2016 for accomodating ANNOVAR-formated step3 previous file and new variant ID format.
aric_freqFilter <- function(previous,freqFile,th=0.0001,step=4.1) {
	th <- as.numeric(th)
	proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
	proID <- substr(previous,proID.start,proID.start+8)
<<<<<<< HEAD
	freq <- read.table(freqFile,header=T,as.is=T,sep='\t',quote='')
	prev <- read.delim(previous,header=T,as.is=T,colClasses=c(Ref='character',Alt='character'),sep='\t',quote='') # previous file is ANNOVAR-streamlined step3 file.
	varID <- paste('var',prev[,1],prev[,2],prev[,3],prev[,4],prev[,5],sep='_') 
	prev <- data.frame(varID=varID,prev[,6:ncol(prev)])
	var2freq <- data.frame(varID=freq$varID,freqAA=freq$ARIC_AA,freqEA=freq$ARIC_EA)
	curr <- merge(prev,var2freq,by.x='varID',by.y='varID',all.x=T)
	judgeRow <- function(x,th) {
		x=as.numeric(x)
		y=sapply(x,function(x,th) is.na(x)|x<=th,th)
		any(y)
	}
	idx <- apply(curr[,(ncol(curr)-1):ncol(curr)],1,judgeRow,th=th)
	curr <- curr[idx,]
	colnames(curr)[(ncol(curr)-1):ncol(curr)] <- c('ARIC_freqAA','ARIC_freqEA')
	write.table(curr,paste0(proID,'_step',step,'.txt'),row.names=F,sep='\t',quote=F)
=======
        freq <- read.table(freqFile,header=T,as.is=T,sep='\t',quote='')
        prev <- read.table(previous,header=T,as.is=T,sep='\t',quote='') # previous file is ANNOVAR-streamlined step3 file.
        varID <- paste('var',prev[,1],prev[,2],prev[,3],prev[,4],prev[,5],sep='_') 
        prev <- data.frame(varID=varID,prev[,6:ncol(prev)])
        var2freq <- data.frame(varID=freq$varID,freqAA=freq$ARIC_AA,freqEA=freq$ARIC_EA)
        curr <- merge(prev,var2freq,by.x='varID',by.y='varID',all.x=T)
        judgeRow <- function(x,th) {
                x=as.numeric(x)
                y=sapply(x,function(x,th) is.na(x)|x<=th,th)
                any(y)
        }
        idx <- apply(curr[,(ncol(curr)-1):ncol(curr)],1,judgeRow,th=th)
        curr <- curr[idx,]
        colnames(curr)[(ncol(curr)-1):ncol(curr)] <- c('ARIC_freqAA','ARIC_freqEA')
        write.table(curr,paste0(proID,'_step',step,'.txt'),row.names=F,sep='\t',quote=F)
>>>>>>> origin/master
}

args=(commandArgs(TRUE))
previous=args[[1]]
freqFile=args[[2]]
th=args[[3]]
cmdTxt <- paste0('aric_freqFilter("',previous,'","',freqFile,'",',th,')')
cat('COMMAND: ',cmdTxt,'\n')
eval(parse(text=cmdTxt))

