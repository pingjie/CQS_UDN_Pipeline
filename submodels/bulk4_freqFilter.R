options(echo=F)
### previous file also includes ANNOVAR annotation columns.
### UPGRADED 3/30/2017 for taking into account of gnomad_exome and gnomad_genome
### UPDATED 4/20/2017: function name changed to bulk4_freqFilter to reflect 4 filtering references: 1000Genome, ESP, ExAC, and gnomAD.
### UPDATED 6/4/2017: keep additional columns "snp138",'gnomadAF','gnomadHom','gnomadHemi' (last 3 coming from gnomad_freqFilter.
bulk4_freqFilter <- function(previous,th=0.01) {
	#### extract udnID ##########
	th=as.numeric(th)
	udnID.match <- gregexpr('UDN\\d+',previous,perl=T,ignore.case=T)
	match.start <- unlist(udnID.match)
	match.len <- attr(udnID.match[[1]],'match.length')
	udnID <- substr(previous,match.start,match.len)
	#### end extracting ############
	prev <- read.delim(previous,as.is=T,sep='\t',quote='') ### assume header is present.
        annovar <- prev
        exacCol <- as.numeric(annovar[,'ExAC_ALL'])
        espCol <- as.numeric(annovar[,'esp6500siv2_all'])
        g1kCol <- as.numeric(annovar[,'X1000g2015aug_all'])
        gnomAD.gCol <- as.numeric(annovar[,'gnomAD_genome_ALL'])
        gnomAD.eCol <- as.numeric(annovar[,'gnomAD_exome_ALL'])
        exac <- is.na(exacCol) | exacCol<=th
        esp <- is.na(espCol) | espCol <=th
        g1k <- is.na(g1kCol) | g1kCol <=th
        gnomAD <- (is.na(gnomAD.eCol) | gnomAD.eCol<=th) & (is.na(gnomAD.gCol) | gnomAD.gCol<=th)
        annovar <- annovar[exac&esp&g1k&gnomAD,c('qual','depth','coverage','GQ','snp138','gnomadAF','gnomadHom','gnomadHemi','ExAC_ALL','esp6500siv2_all','X1000g2015aug_all','gnomAD_exome_ALL','gnomAD_genome_ALL')] #modified 2/21/2017
	colnames(annovar)[(ncol(annovar)-4):ncol(annovar)] <- c('exacFreq','espFreq','g1kFreq','gnomADexomeFreq','gnomADgenomeFreq') # modified 2/21/2017
	exacMetrics <- prev[exac&esp&g1k&gnomAD,c('mis_z','pLI','Otherinfo')] # added 2/21/2017
	colnames(exacMetrics)[ncol(exacMetrics)] <- 'GT' # added 2/21/2017
        curr <- cbind(prev[exac&esp&g1k&gnomAD,1:6],annovar,exacMetrics) #modified 2/21/2017
        write.table(curr,paste0(udnID,'_step4.txt'),row.names=F,sep='\t',quote=F)
}

args=(commandArgs(TRUE))
for (i in 1:length(args)) {
        cat(args[[i]],'\n')
}


if (length(args)==2) {
        th=args[[2]]
} else {
        th=0.01
}

bulk4_freqFilter(args[[1]],th)

