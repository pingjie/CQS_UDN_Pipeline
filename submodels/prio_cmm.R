### INPUT allRows (the same for 4 following functions): rows under the name of one particular candidate gene, from full file '*.annotated_gene_scores'
### INPUT seedListFile='*.seed_gene_list'
sum.biosystem <- function(allRows,seedListFile) {
	colnames(allRows) <- c('sourceID','sourceName','seed','relScore')
	interstRows <- allRows[grepl('(BIOSYSTEM)',allRows$sourceID,fixed=T),]
	if (nrow(interstRows)>0) {
		seeds <- sapply(strsplit(interstRows$seed,' ',fixed=T), function(x) x[length(x)])
		SEEDS <- read.delim(seedListFile,as.is=T)[,2]
		matched <- match(SEEDS,seeds)
		seedHit <- SEEDS[which(!is.na(matched))][1]
		top_seed_state <- paste(interstRows[matched[!is.na(matched)][1],c(2:3)],collapse=' ')
		top_seed_state <- gsub('In the same','In the same pathways',top_seed_state)
		### identify the most frequent pathway #####
		pwNames <- gsub('In the same (','',interstRows$sourceName,fixed=T); pwNames <- gsub(')','',pwNames,fixed=T)
		pwNames <- strsplit(pwNames,'; ')
		tmp <- table(unlist(pwNames))
		top1pw <- names(tmp)[which.max(tmp)[1]]
		top1pw.rows <- grep(top1pw,interstRows$sourceName,fixed=T)
		### extract biosystem ID for the most frequent pathway ###
		order_in_row <- which(pwNames[[top1pw.rows[1]]]==top1pw)[1]
		sysIDs <- gsub( 'BIOSYSTEM:','',interstRows$sourceID[ top1pw.rows[1] ],fixed=T )	
		sysIDs <- gsub( ' (BIOSYSTEM)','',sysIDs,fixed=T )
		sysIDs <- unlist(strsplit(sysIDs,' ',fixed=T))
		sysID <- sysIDs[order_in_row]
		### summarize how many seeds appear in this pathway ###	
		seeds.samePW <- seeds[top1pw.rows]
		matched <- match(SEEDS,seeds.samePW)
		top1pw.topseed <- seeds.samePW[matched[!is.na(matched)][1]] ## the highest seed belonging to the most frequent pathway.
		seedFrac <- sum(!is.na(matched))/length(SEEDS)
		seedN <- sum(!is.na(matched))
		top_pw_state <- paste0('Most important pathway (',top1pw,', BIOSYSTEM:',sysID,') covers ',seedN,' or ', round(seedFrac,2),' seedGenes, say, ',top1pw.topseed)
		biosystem.states <- c(top_seed_state,top_pw_state)
		biosystem.states	
		#seedHit <- SEEDS[which(!is.na(matched))][1]
		list(state=paste(biosystem.states,collapse='; '),seed=seedHit)
	} else {
		list(state=NA,seed=NA)
	}
}

sum.HTRI <- function(allRows,seedListFile) {
        colnames(allRows) <- c('sourceID','sourceName','seed','relScore')
	interstRows <- allRows[grepl('(HTRI)',allRows$sourceID,fixed=T),]
	if (nrow(interstRows)>0) {
		seeds <- sapply(strsplit(interstRows$seed,' ',fixed=T), function(x) x[length(x)])
		SEEDS <- read.delim(seedListFile,as.is=T)[,2]
		matched <- match(SEEDS,seeds)
		hitRow <- interstRows[matched[!is.na(matched)][1],]	
		state <- paste0(hitRow[,3],' in experiment (',hitRow[,2],')',', ',gsub(' (HTRI)','',hitRow[,1],fixed=T),collapse='')
		seedHit <- SEEDS[which(!is.na(matched))][1]
		list(state=state,seed=seedHit)
	} else {
                list(state=NA,seed=NA)
	}
}

sum.family <- function(allRows,seedListFile) {
        colnames(allRows) <- c('sourceID','sourceName','seed','relScore')
        interstRows <- allRows[grepl('(GENE_FAMILY)',allRows$sourceID,fixed=T),]
	if (nrow(interstRows)>0) {
		seeds <- sapply(strsplit(interstRows$seed,' ',fixed=T), function(x) x[length(x)])
		SEEDS <- read.delim(seedListFile,as.is=T)[,2]
		matched <- match(SEEDS,seeds)
		hitRow <- interstRows[matched[!is.na(matched)][1],]
		seedHit <- SEEDS[which(!is.na(matched))][1]
		state <- paste0(hitRow[,2:3],collapse=' ')
		list(state=state,seed=seedHit)
        } else {
                list(state=NA,seed=NA)
        }

}

sum.HPRD <- function(allRows,seedListFile) {
        colnames(allRows) <- c('sourceID','sourceName','seed','relScore')
        interstRows <- allRows[grepl('(HPRD)',allRows$sourceID,fixed=T),]
        if (nrow(interstRows)>0) {
		seeds <- sapply(strsplit(interstRows$seed,' ',fixed=T), function(x) x[length(x)])
		SEEDS <- read.delim(seedListFile,as.is=T)[,2]
		matched <- match(SEEDS,seeds) # pick the highest-ranking seed gene as evidence.
		hitRow <- interstRows[matched[!is.na(matched)][1],]
		seedHit <- SEEDS[which(!is.na(matched))][1]
		state <- paste0(c('In protein-protein network interacts',hitRow[,c(2,3)],',',hitRow[,1]),collapse=' ')
		state=gsub(' (HPRD)','',state,fixed=T)
		list(state=state,seed=seedHit)
        } else {
                list(state=NA,seed=NA)
        }

}

### seedDetailFile: outTerms.merge_gene_scores
sumSeed <- function(seed,seedDetailFile) {
	seedTbl <- read.delim(seedDetailFile,skip=1,as.is=T,row.names=NULL,head=F)
	seedRows <- which(is.na(seedTbl[,4])) # Row Numbers for each gene overview
	names(seedRows) <- seedTbl[seedRows,1]
	startRow <- seedRows[seed]+1
	if (names(seedRows)[length(seedRows)]!=seed)
		stopRow <- seedRows[which(names(seedRows)==seed)+1]-1
	else
		stopRow <- nrow(seedTbl)
	seedTbl <- seedTbl[startRow:stopRow,] ## subrows of full table, for the gene 'seed' only.
	if (sum(grepl('OMIM:',seedTbl[,1],fixed=T))>0) seedTbl <- seedTbl[grepl('OMIM:',seedTbl[,1]),] ## just report OMIM-related evidence, if available.
	diseases <- paste(unique(seedTbl[,2]),collapse=', ') ## unique diseases
	phenos <- paste(unique(seedTbl[,3]),collapse=', ') ## unique features
	refs <- gsub(' \\(.*\\)','',seedTbl[,1],perl=T)
	omims <- unique(gsub('OMIM:','',refs[grepl('OMIM:',refs)]))
	orphanet <- unique(gsub('ORPHANET:','',refs[grepl('ORPHANET:',refs)]))
        pubmed <- unique(gsub('PUBMED:','',refs[grepl('PUBMED:',refs)]))	
	concat <- function(x,str) {
		if (length(x)>0)
			paste(c(str,x),collapse=':')
		else
			x
	}
	sources <- paste(c(concat(omims,'OMIM'),concat(orphanet,'ORPHANET'),concat(pubmed,'PUBMED')),collapse=',') ## unique source IDs
	state <- paste0(seed,' is associated with diseases (',diseases,'), covering features (',phenos,'); ',sources,collapse='')
	state	
}
### fullGeneFile='outTerms.annotated_gene_scores'
get.geneFile <- function(gene,fullGeneFile) {
	fullTbl <- read.delim(fullGeneFile,skip=1,as.is=T,row.names=NULL,head=F)	
	geneRows <- grep('Normalized score:',fullTbl[,4])
	names(geneRows) <- fullTbl[geneRows,1]
        startRow <- geneRows[gene]+1
	if (gene!=names(geneRows)[length(geneRows)])
        	stopRow <- geneRows[which(names(geneRows)==gene)+1]-1
	else
		stopRow <- nrow(fullTbl)
        geneTbl <- fullTbl[startRow:stopRow,]
	geneTbl
}

#####################################################
################ Running scripts below ##############
#####################################################
args=commandArgs(TRUE)
## 1st arg: Phenolyzer output prefix
## 2nd arg: path/udnID. A file named args[2]_finalFiltered.txt must be existent in the current wd.
prefix=args[1]
udnID=args[2]
gene_list_file=paste0(prefix,'.annotated_gene_list')#args[1] #gene_list_file
gene_detail_file=paste0(prefix,'.annotated_gene_scores')#args[2] #gene_detail_file
seed_list_file=paste0(prefix,'.seed_gene_list')#args[3] #seed_list_file
seed_detail_file=paste0(prefix,'.merge_gene_scores')#args[4] #seed_detail_file
#source('sum.geneDetail.R')
candTbl <- read.delim(gene_list_file,as.is=T)
genes <- candTbl$Gene
statements <- character(length(genes))
for (i in 1:length(genes)) {
	gene <- genes[i]
	#cat(i,'th gene',gene,'...\n')
	if (candTbl$Status[i]=='Predicted') {
		geneRows <- get.geneFile(gene,gene_detail_file)
		systems <- sum.biosystem(geneRows,seed_list_file)
		htri <- sum.HTRI(geneRows,seed_list_file)
		family <- sum.family(geneRows,seed_list_file)
		hprd <- sum.HPRD(geneRows,seed_list_file)
		seeds <- sort(unique(c(systems$seed,htri$seed,family$seed,hprd$seed)))
		states <- c(family$state,htri$state,hprd$state,systems$state)
		SEEDS <- list(family$seed,htri$seed,hprd$seed,systems$seed)
		nonNA <- sapply(SEEDS,function(x) length(setdiff(x,NA))>0)
		states <- states[nonNA]
		geneState <- paste(states,collapse=' || ')
		seedState <- paste(sapply(seeds,sumSeed,seed_detail_file),collapse=' || ')
		statement <- paste(c(geneState,seedState),collapse=' ... ')
	} else {
		statement <- sumSeed(gene,seed_detail_file)
	}
	statements[i] <- statement
}

candTbl <- data.frame(candTbl,comment=statements)
seedGenes <- candTbl$Gene[candTbl$Status=='SeedGene']
filteredTbl <- read.delim(paste0(udnID,'_finalFiltered.txt'),as.is=T) ## read in udnID_finalFiltered.txt file.
if (is.element('ARch.tier2',filteredTbl$ihMode)) { ## Filter out useless ARch.tier2 variants, keeping only Phenolyzer seedGene variants.
	if (sum(candTbl$Status=='SeedGene')>0)
		filteredTbl <- filteredTbl[filteredTbl$ihMode!='ARch.tier2'|filteredTbl$geneName%in%seedGenes,]
	else
		filteredTbl <- filteredTbl[filteredTbl$ihMode!='ARch.tier2',]
} 
cmmedTbl <- merge(filteredTbl,candTbl,by.x='geneName',by.y='Gene',all.x=T)
colnames(cmmedTbl)[1] <- 'geneName'
if (  sum(grepl(',',cmmedTbl$geneName))>0  ) {
	indices <- grep(',',cmmedTbl$geneName)
	for (i in indices) {
		gene1st <- unlist(strsplit(cmmedTbl$geneName[i],','))[1]
		if (is.element(gene1st,candTbl$Gene)) {
			ix <- match(gene1st,candTbl$Gene)
			cmmedTbl[i,(ncol(cmmedTbl)-4):ncol(cmmedTbl)] <- candTbl[ix,c(1,3:6)]
		}
	}
}
write.table(cmmedTbl,paste0(udnID,'_prioCmmed.txt'),row.names=F,sep='\t',quote=F)
