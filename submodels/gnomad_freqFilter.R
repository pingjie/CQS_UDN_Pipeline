require(data.table)

#load("/scratch/cqs/udn/gnomad_freq.rdata")

options(echo=F)
### UPDATE 6/4/2017: modified based on exac_freqFilter().
### Now work on freqFile=gnomad_UDN000000_exomes.txt; but can be working on gnomad_UDN000000_combined.txt if the file is obtained with the same headers.
gnomad_freqFilter <- function(previous, th=0.001, n_Homo=5, gnomaddata,step='4.3') { # Only refer to Num_of_homo in the exac file. But still call it a freqFile.
	th <- as.numeric(th)
	n_Homo <- as.numeric(n_Homo)
	proID.start <- regexpr('UDN\\d+', previous, perl = T, ignore.case = T)
	proID <- substr(previous, proID.start, proID.start + 8)
	prev <- read.table(previous, header = T, as.is = T, sep = '\t', quote = '')
	load(gnomaddata) # this cause an object "freq" to occur in working space.
	colnames(freq) <- c('varID','gnomadAF','gnomadHom','gnomadHemi')
	curr <- merge(prev, freq, by.x = 'varID', by.y = 'varID', all.x = T)	
	curr$gnomadHom <- as.numeric(curr$gnomadHom)
	curr$gnomadHemi <- as.numeric(curr$gnomadHemi)
	curr$gnomadAF <- as.numeric(curr$gnomadAF)
	curr <- curr[(curr$gnomadHom < n_Homo | is.na(curr$gnomadHom)) & (curr$gnomadHemi < n_Homo | is.na(curr$gnomadHemi)), ]
	curr <- curr[curr$gnomadAF <= th | is.na(curr$gnomadAF), ]
	write.table(curr, paste0(proID, '_step', step, '.txt'), row.names = F, sep = '\t', quote = F)
}

args <- commandArgs(TRUE)
previous <- args[[1]]
th <- args[[2]]
n_Homo <- args[[3]]
gnomaddata <- args[[4]]

gnomad_freqFilter(previous, th, n_Homo,gnomaddata)
