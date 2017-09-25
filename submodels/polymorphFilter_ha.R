require(data.table)

load("/scratch/cqs/udn/gnomad_freq.rdata")

options(echo=F)
### previous file also includes ANNOVAR annotation columns.
polymorphFilter <- function(previous, th = 0.01, n_Homo = 5) {
	#### extract udnID ##########
	th <- as.numeric(th)
	udnID.match <- gregexpr('UDN\\d+', previous, perl = T, ignore.case = T)
	match.start <- unlist(udnID.match)
	match.len <- attr(udnID.match[[1]], 'match.length')
	udnID <- substr(previous, match.start, match.len)
	#### end extracting ############
	prev <- read.delim(previous, as.is = T, colClasses = c(Ref = 'character', Alt = 'character'), sep = '\t', quote = '') ### assume header is present.
	varID <- paste('var', prev[, 1], prev[, 2], prev[, 3], prev[, 4], prev[, 5], sep = '_')
	varID.fct <- factor(varID, levels = unique(varID), ordered = T)
	#prev <- do.call(rbind,by(prev,varID.fct,function(x) x[1,])) ### for weird case where redundant varID occur (transpassing exonic & rna). #mutated 3/23/2017 as ncRNA part has been muted.
	row.names(prev) <- unique(varID)
	annovar <- prev
	exacCol <- as.numeric(annovar[, 'ExAC_ALL'])
	espCol <- as.numeric(annovar[, 'esp6500siv2_all'])
	g1kCol <- as.numeric(annovar[, 'X1000g2015aug_all'])
	gnomAD.gCol <- as.numeric(annovar[, 'gnomAD_genome_ALL'])
	gnomAD.eCol <- as.numeric(annovar[, 'gnomAD_exome_ALL'])
	exac <- is.na(exacCol) | exacCol <= th
	esp <- is.na(espCol) | espCol <= th
	g1k <- is.na(g1kCol) | g1kCol <= th
	gnomAD <- (is.na(gnomAD.eCol) | gnomAD.eCol <= th) & (is.na(gnomAD.gCol) | gnomAD.gCol <= th)
	annovar <- annovar[exac&esp&g1k&gnomAD, c('qual','depth','coverage','GQ','snp138','ExAC_ALL','esp6500siv2_all','X1000g2015aug_all','gnomAD_exome_ALL','gnomAD_genome_ALL')]
	colnames(annovar)[(ncol(annovar)-4):ncol(annovar)] <- c('exacFreq','espFreq','g1kFreq','gnomADexomeFreq','gnomADgenomeFreq') 
	exacMetrics <- prev[exac&esp&g1k&gnomAD,c('mis_z','pLI','Otherinfo')]
	colnames(exacMetrics)[ncol(exacMetrics)] <- 'GT' 
	curr <- cbind(varID=rownames(prev)[exac&esp&g1k&gnomAD],prev[exac&esp&g1k&gnomAD,6:10],annovar,exacMetrics)
	write.table(curr,paste0(udnID,'_step4.1.txt'),row.names=F,sep='\t',quote=F)
	step4.2 <- gnomad_freqFilter(paste0(udnID,'_step4.1.txt'),th=th,n_Homo=n_Homo,step='4.2') 
	snp138.ix <- which(colnames(step4.2)=='snp138')
	curr1 <- step4.2[,1:snp138.ix]
	curr2 <- step4.2[,(snp138.ix+1):(ncol(step4.2)-3)]
	curr <- cbind(curr1,step4.2[,c('gnomadAF','gnomadHom','gnomadHemi')],curr2) #step4.2[,1:(ncol(step4.2)-2)]#:ncol(step4.2)]
	write.table(curr,paste0(udnID,'_step4.txt'),row.names=F,sep='\t',quote=F)
}

gnomad_freqFilter <- function(previous,th=0.001,n_Homo=5,step='4.3') { # Only refer to Num_of_homo in the exac file. But still call it a freqFile.
	th <- as.numeric(th)
	n_Homo <- as.numeric(n_Homo)
	proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
	proID <- substr(previous,proID.start,proID.start+8)
	prev <- read.table(previous,header=T,as.is=T,sep='\t',quote='')
	curr <- merge(prev, freq, by.x = 'varID', by.y = 'varID', all.x = T)
	curr <- curr[(curr$gnomadHom<n_Homo|is.na(curr$gnomadHom))&(curr$gnomadHemi<n_Homo|is.na(curr$gnomadHemi)),]
	curr <- curr[curr$gnomadAF<=th|is.na(curr$gnomadAF),]
	write.table(curr,paste0(proID,'_step',step,'.txt'),row.names=F,sep='\t',quote=F)
	curr
}

args <- commandArgs(TRUE)
for (i in 1:length(args)) {
        cat(args[[i]],'\n')
}
if (length(args)==2) {
        th <- args[[2]]
} else {
        th <- 0.01
}
if (length(args)==3) {
        n_Homo <- args[[3]]
} else {
        n_Homo <- 5
}

polymorphFilter(args[[1]], th, n_Homo)

